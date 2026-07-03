import glob
import re
from pathlib import Path

import duckdb
import pandas as pd


# Taxonomy rank columns for the split OTU table
TAXONOMY_RANKS = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Details']


def natural_sort_key(s):
    match = re.search(r'(\d+)', str(s))
    if match:
        return (0, int(match.group(1)))
    return (1, str(s))


def _split_taxonomy_series(taxonomy_series: pd.Series, db_format: str = 'unite') -> pd.DataFrame:
    """Split a Taxonomy Series (semicolon-separated) into 8 rank columns.

    Rank positions by format:
      unite / eukariome / cbs:  Kingdom;Phylum;Class;Order;Family;Genus;Species
      silva:                    Domain(=Kingdom);Phylum;Class;Order;Family;Genus[;Species]
      none:                     best-effort positional assignment

    Species is always stored as binomial (Genus + epithet).
    Anything beyond 7 core ranks is collected into 'Details'.
    """

    def _split_one(tax_str):
        if not isinstance(tax_str, str) or not tax_str.strip():
            return [''] * 8

        # Replace underscores with spaces (databases like UNITE use
        # underscores as word separators, e.g. "Tuber_melanosporum")
        parts = [p.strip().replace('_', ' ') for p in tax_str.split(';') if p.strip()]

        result = [''] * 8  # Kingdom, Phylum, Class, Order, Family, Genus, Species, Details

        # Assign parts to rank positions (max 7 core ranks)
        for i in range(min(len(parts), 7)):
            result[i] = parts[i]

        # Anything beyond 7 parts → Details
        if len(parts) > 7:
            result[7] = ';'.join(parts[7:])

        # --- Ensure Species is binomial ---
        genus = result[5]
        species = result[6]
        if species:
            sp_words = species.split()
            if len(sp_words) == 1 and genus:
                # Single epithet → prepend genus for binomial
                result[6] = f"{genus} {species}"
            elif len(sp_words) > 2:
                # More than binomial → keep first two words, rest → Details
                result[6] = f"{sp_words[0]} {sp_words[1]}"
                extras = ' '.join(sp_words[2:])
                result[7] = f"{extras};{result[7]}" if result[7] else extras
            # len == 2: already binomial, keep as-is

        return result

    split_data = taxonomy_series.apply(_split_one)
    return pd.DataFrame(
        split_data.tolist(),
        columns=TAXONOMY_RANKS,
        index=taxonomy_series.index
    )


def _build_reformat_sql(db_format: str, flat_sh: bool = False) -> str:
    """Return a DuckDB SQL expression that extracts the canonical OTU_ID
    from the raw index column, matching the logic in reformat_tmp_indices.

    This replaces Step 5 (reformat_tmp_indices) entirely — the OTU ID
    reformatting now happens inside the DuckDB query, avoiding a full
    read-write cycle on every .tmp.txt file.

    flat_sh (UNITE only): collapse rows to the SH species hypothesis code,
    so multiple reference accessions sharing the same SH become a single
    OTU row.
    """
    if db_format == 'unite':
        if flat_sh:
            # "species|accession|SH|refs|tax" → "SH" (parts[2])
            return "split_part(otu_id, '|', 3)"
        # "species|accession|SH|refs|tax" → "accession|SH" (parts[1]|parts[2])
        return ("split_part(otu_id, '|', 2) || '|' || split_part(otu_id, '|', 3)")
    elif db_format == 'eukariome':
        # "accession;taxonomy..." → "accession"
        return "split_part(otu_id, ';', 1)"
    elif db_format in ('cbs', 'none'):
        # "accession|extra" → "accession"
        return "split_part(otu_id, '|', 1)"
    else:
        # silva and others: keep as-is
        return "otu_id"


def _build_taxonomy_split_sql() -> str:
    """Return DuckDB SQL columns that split a semicolon-separated Taxonomy
    string into 8 rank columns (Kingdom..Details).

    This replaces _split_taxonomy_series (Python row-by-row apply) with
    a single vectorized SQL expression evaluated by DuckDB's engine.

    Binomial species logic:
      - 1 word epithet → prepend Genus (e.g. "melanosporum" → "Genus melanosporum")
      - >2 words → keep first two, rest to Details
      - Underscores replaced with spaces (UNITE convention)
    """
    # Helper: shorthand for extracting and cleaning the Nth taxonomy part
    # DuckDB arrays are 1-indexed; COALESCE handles missing parts
    def _part(n):
        return f"trim(replace(COALESCE(string_split(Taxonomy, ';')[{n}], ''), '_', ' '))"

    # Count spaces in species = number of words - 1
    sp = _part(7)
    sp_nspaces = f"(length({sp}) - length(replace({sp}, ' ', '')))"

    return f"""
        {_part(1)} AS Kingdom,
        {_part(2)} AS Phylum,
        {_part(3)} AS "Class",
        {_part(4)} AS "Order",
        {_part(5)} AS Family,
        {_part(6)} AS Genus,
        -- Species: ensure binomial
        CASE
            WHEN {sp} = '' THEN ''
            WHEN {sp_nspaces} = 0 AND {_part(6)} != ''
                THEN {_part(6)} || ' ' || {sp}
            WHEN {sp_nspaces} >= 2
                THEN string_split({sp}, ' ')[1] || ' ' || string_split({sp}, ' ')[2]
            ELSE {sp}
        END AS Species,
        -- Details: extra species words + parts beyond rank 7
        CASE
            WHEN {sp_nspaces} >= 2 AND array_length(string_split(Taxonomy, ';')) > 7
                THEN array_to_string(string_split({sp}, ' ')[3:], ' ')
                     || ';' || array_to_string(string_split(Taxonomy, ';')[8:], ';')
            WHEN {sp_nspaces} >= 2
                THEN array_to_string(string_split({sp}, ' ')[3:], ' ')
            WHEN array_length(string_split(Taxonomy, ';')) > 7
                THEN array_to_string(string_split(Taxonomy, ';')[8:], ';')
            ELSE ''
        END AS Details
    """


def makeOtu_duckdb(results_dir: str, output_file: str, taxonomy_file: str,
                   db_format: str = 'unite', flat_sh: bool = False):
    results_path = Path(results_dir)
    output_path = Path(output_file)

    input_glob_pattern = str(results_path / '*.tmp.txt')
    all_files_raw = glob.glob(input_glob_pattern)
    # OPT-B4: normalise backslashes to forward slashes so that the
    # DuckDB `split_part(filename, '/', -1)` expression below extracts the
    # sample name correctly on Windows (native, non-Docker) runs.
    clean_files = [
        f.replace('\\', '/')
        for f in all_files_raw
        if not Path(f).name.startswith('._')
    ]

    if not clean_files:
        print(f"No valid .tmp.txt files found in {results_dir}")
        return

    print(f"Starting count aggregation with DuckDB ({len(clean_files)} valid files)...")

    col_types = {
        'otu_id': 'VARCHAR', 'length': 'BIGINT',
        'mapped_count': 'BIGINT', 'unmapped_count': 'BIGINT'
    }

    col_types_sql = "{" + ", ".join(f"'{k}': '{v}'" for k, v in col_types.items()) + "}"
    clean_files_sql = "[" + ", ".join(
        "'" + f.replace("'", "''") + "'" for f in clean_files
    ) + "]"

    # --- OPT-A: OTU ID reformatting inside DuckDB ---
    # Replaces Step 5 (reformat_tmp_indices) entirely
    # OPT-B7: optional UNITE Flat-SH aggregation (collapses ref accessions
    # that share the same Species Hypothesis into a single OTU row).
    use_flat_sh = flat_sh and db_format == 'unite'
    reformat_expr = _build_reformat_sql(db_format, flat_sh=use_flat_sh)
    if flat_sh and not use_flat_sh:
        print(f"[Flat-SH] Ignored: only applies to db_format='unite' (current: {db_format})")
    elif use_flat_sh:
        print("[Flat-SH] UNITE references grouped by SH species hypothesis.")

    # --- OPT-5: Taxonomy join + OPT-A split inside DuckDB ---
    tax_path = Path(taxonomy_file) if taxonomy_file else None
    has_tax = tax_path and tax_path.exists()

    # --- OPT-B: Taxonomy split inside DuckDB ---
    tax_split_sql = _build_taxonomy_split_sql()

    if has_tax:
        tax_path_sql = str(tax_path).replace("'", "''")

        # Taxonomy CTE — when Flat-SH is active, project OTU_ID to the SH code
        # (= second pipe-separated field of the cached tax-map key "acc|SH")
        # and dedupe, so multiple accessions sharing an SH collapse to one row.
        tax_read_csv = (
            "read_csv('" + tax_path_sql + "', delim='\\t', header=True, "
            "columns={'OTU_ID': 'VARCHAR', 'Taxonomy': 'VARCHAR'})"
        )
        if use_flat_sh:
            taxonomy_cte_sql = (
                "taxonomy AS (\n"
                "    SELECT split_part(OTU_ID, '|', 2) AS OTU_ID,\n"
                "           any_value(Taxonomy) AS Taxonomy\n"
                "    FROM " + tax_read_csv + "\n"
                "    WHERE split_part(OTU_ID, '|', 2) != ''\n"
                "    GROUP BY split_part(OTU_ID, '|', 2)\n"
                ")"
            )
        else:
            taxonomy_cte_sql = "taxonomy AS (SELECT * FROM " + tax_read_csv + ")"

        sql_query = f"""
        WITH counts AS (
            SELECT * FROM (
                PIVOT (
                    SELECT
                        replace(split_part(filename, '/', -1), '.tmp.txt', '') AS sample_name,
                        {reformat_expr} AS otu_id,
                        mapped_count
                    FROM read_csv(
                        {clean_files_sql}, delim='\\t', header=False,
                        columns={col_types_sql}, filename=True
                    )
                    WHERE otu_id != '*'
                )
                ON sample_name
                USING sum(mapped_count)
                GROUP BY otu_id
            )
        ),
        {taxonomy_cte_sql},
        joined AS (
            SELECT t.OTU_ID, t.Taxonomy, c.* EXCLUDE (otu_id)
            FROM taxonomy t
            INNER JOIN counts c ON t.OTU_ID = c.otu_id
        )
        SELECT OTU_ID, Taxonomy,
               {tax_split_sql},
               j.* EXCLUDE (OTU_ID, Taxonomy)
        FROM joined j;
        """
    else:
        sql_query = f"""
        PIVOT (
            SELECT
                replace(split_part(filename, '/', -1), '.tmp.txt', '') AS sample_name,
                {reformat_expr} AS otu_id,
                mapped_count
            FROM read_csv(
                {clean_files_sql}, delim='\\t', header=False,
                columns={col_types_sql}, filename=True
            )
            WHERE otu_id != '*'
        )
        ON sample_name
        USING sum(mapped_count)
        GROUP BY otu_id;
        """

    final_df = pd.DataFrame()
    try:
        con = duckdb.connect()
        final_df = con.execute(sql_query).df()
        con.close()

        if final_df.empty:
            print("WARNING: No valid OTU data found (all reads may have been "
                  "removed by the quality filter). Creating empty OTU table.")
            empty_df = pd.DataFrame(columns=['Taxonomy'])
            empty_df.index.name = 'OTU_ID'
            try:
                empty_df.to_excel(output_path, engine='openpyxl')
                print(f"Empty OTU table saved to: {output_path}")
            except Exception as e:
                print(f"Error saving empty OTU table: {e}")
            return

        if has_tax:
            # DuckDB join may produce duplicate otu_id variants
            otu_id_dupes = [c for c in final_df.columns
                           if c.lower().startswith('otu_id') and c != 'OTU_ID']
            if otu_id_dupes:
                final_df = final_df.drop(columns=otu_id_dupes)
            final_df = final_df.set_index('OTU_ID')
        else:
            final_df = final_df.set_index('otu_id')
            final_df.index.name = 'OTU_ID'

        # fillna + astype in-place where possible
        final_df = final_df.fillna(0)
        num_cols = final_df.select_dtypes(include='number').columns
        final_df[num_cols] = final_df[num_cols].astype(int)

        # Strip whitespace from taxonomy string columns
        str_cols = final_df.select_dtypes(include='object').columns
        for c in str_cols:
            final_df[c] = final_df[c].str.strip()

        # --- Fallback: split taxonomy in Python if DuckDB didn't produce rank columns ---
        if 'Taxonomy' in final_df.columns and 'Kingdom' not in final_df.columns:
            print(f"Splitting taxonomy into rank columns (format: {db_format})...")
            rank_df = _split_taxonomy_series(final_df['Taxonomy'], db_format=db_format)
            for i, col in enumerate(TAXONOMY_RANKS):
                final_df.insert(
                    final_df.columns.get_loc('Taxonomy') + 1 + i,
                    col, rank_df[col]
                )

        # Natural-sort sample columns
        print("Applying natural sort to sample columns (barcode01, barcode02, ...)...")
        meta_cols = ['Taxonomy'] + TAXONOMY_RANKS
        tax_cols = [c for c in meta_cols if c in final_df.columns]
        sample_columns = [c for c in final_df.columns if c not in tax_cols]
        sorted_sample_columns = sorted(sample_columns, key=natural_sort_key)
        final_df = final_df[tax_cols + sorted_sample_columns]

    except Exception as e:
        print(f"CRITICAL ERROR during DuckDB count aggregation: {e}")
        import traceback
        traceback.print_exc()
        return

    # --- OPT-A5: Save TSV first (instant), then XLSX for user convenience ---
    tsv_output = output_path.with_suffix('.tsv')
    print(f"Saving final table (TSV) to {tsv_output}...")
    try:
        final_df.to_csv(tsv_output, sep='\t')
        print(f"Final OTU table (TSV) saved to: {tsv_output}")
    except Exception as e:
        print(f"WARNING: Error saving TSV file: {e}")

    # --- OPT-3: xlsxwriter is 5-15x faster than openpyxl for write-only ---
    print(f"Saving final table (XLSX) to {output_path}...")
    try:
        final_df.to_excel(output_path, engine='xlsxwriter')
        print(f"Final OTU table (XLSX) saved to: {output_path}")
    except ImportError:
        try:
            final_df.to_excel(output_path, engine='openpyxl')
            print(f"Final OTU table saved to: {output_path} (using openpyxl fallback)")
        except ImportError:
            print("Neither 'xlsxwriter' nor 'openpyxl' found. XLSX skipped (TSV already saved).")
    except Exception as e:
        print(f"WARNING: Error saving XLSX file: {e} (TSV already saved)")
