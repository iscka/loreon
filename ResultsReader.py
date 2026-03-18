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


def makeOtu_duckdb(results_dir: str, output_file: str, taxonomy_file: str,
                   db_format: str = 'unite'):
    results_path = Path(results_dir)
    output_path = Path(output_file)

    input_glob_pattern = str(results_path / '*.tmp.txt')
    all_files_raw = glob.glob(input_glob_pattern)
    clean_files = [f for f in all_files_raw if not Path(f).name.startswith('._')]

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

    # --- OPT-5: Taxonomy join inside DuckDB ---
    # When a taxonomy file is available, perform the join inside DuckDB
    # before materializing the Python DataFrame — avoids loading the
    # full taxonomy TSV (500k+ rows for SILVA) into pandas.
    tax_path = Path(taxonomy_file) if taxonomy_file else None
    has_tax = tax_path and tax_path.exists()

    if has_tax:
        tax_path_sql = str(tax_path).replace("'", "''")
        sql_query = f"""
        WITH counts AS (
            SELECT * FROM (
                PIVOT (
                    SELECT
                        replace(split_part(filename, '/', -1), '.tmp.txt', '') AS sample_name,
                        otu_id, mapped_count
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
        taxonomy AS (
            SELECT * FROM read_csv(
                '{tax_path_sql}', delim='\\t', header=True,
                columns={{'OTU_ID': 'VARCHAR', 'Taxonomy': 'VARCHAR'}}
            )
        )
        SELECT t.OTU_ID, t.Taxonomy, c.* EXCLUDE (otu_id)
        FROM taxonomy t
        INNER JOIN counts c ON t.OTU_ID = c.otu_id;
        """
    else:
        sql_query = f"""
        PIVOT (
            SELECT
                replace(split_part(filename, '/', -1), '.tmp.txt', '') AS sample_name,
                otu_id, mapped_count
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
            # DuckDB join may produce duplicate otu_id variants (otu_id, otu_id_1, etc.)
            # The EXCLUDE clause should prevent this, but drop any leftovers to be safe
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
        # Convert only numeric columns to int (skip 'Taxonomy' if present)
        num_cols = final_df.select_dtypes(include='number').columns
        final_df[num_cols] = final_df[num_cols].astype(int)

        # --- Split Taxonomy into 8 rank columns ---
        if 'Taxonomy' in final_df.columns:
            print(f"Splitting taxonomy into rank columns (format: {db_format})...")
            rank_df = _split_taxonomy_series(final_df['Taxonomy'], db_format=db_format)
            # Insert rank columns right after Taxonomy
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
        return

    # --- OPT-3: xlsxwriter is 5-15x faster than openpyxl for write-only ---
    print(f"Saving final table to {output_path}...")
    try:
        final_df.to_excel(output_path, engine='xlsxwriter')
        print(f"Final OTU table (with taxonomy) saved to: {output_path}")
    except ImportError:
        # Fall back to openpyxl if xlsxwriter not installed
        try:
            final_df.to_excel(output_path, engine='openpyxl')
            print(f"Final OTU table saved to: {output_path} (using openpyxl fallback)")
        except ImportError:
            print("Neither 'xlsxwriter' nor 'openpyxl' found. Saving as CSV.")
            csv_output = output_path.with_suffix('.csv')
            final_df.to_csv(csv_output)
            print(f"Final OTU table saved to: {csv_output}")
    except Exception as e:
        print(f"Error saving final file {output_path}: {e}")
