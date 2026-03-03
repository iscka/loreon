import glob
import re
from pathlib import Path

import duckdb
import pandas as pd


def natural_sort_key(s):
    match = re.search(r'(\d+)', str(s))
    if match:
        return (0, int(match.group(1)))
    return (1, str(s))


def makeOtu_duckdb(results_dir: str, output_file: str, taxonomy_file: str):
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

    sql_query = f"""
    WITH all_data AS (
        SELECT
            replace(split_part(filename, '/', -1), '.tmp.txt', '') AS sample_name,
            otu_id, mapped_count
        FROM read_csv(
            {clean_files_sql}, delim='\t', header=False,
            columns={col_types_sql}, filename=True
        )
        WHERE otu_id != '*'
    )
    PIVOT ( SELECT * FROM all_data )
    ON sample_name
    USING sum(mapped_count)
    GROUP BY otu_id;
    """

    counts_df = pd.DataFrame()
    try:
        con = duckdb.connect()
        counts_df = con.execute(sql_query).df()
        con.close()

        if counts_df.empty:
            print("No valid OTU data found (only unmapped?).")
            return

        counts_df = counts_df.set_index('otu_id').fillna(0).astype(int)

        print("Applying natural sort to sample columns (barcode01, barcode02, ...)...")
        sample_columns = counts_df.columns.tolist()
        sorted_sample_columns = sorted(sample_columns, key=natural_sort_key)
        counts_df = counts_df[sorted_sample_columns]

    except Exception as e:
        print(f"CRITICAL ERROR during DuckDB count aggregation: {e}")
        return

    tax_df = None
    if taxonomy_file:
        tax_path = Path(taxonomy_file)
        print(f"Loading taxonomy file: {tax_path.name}")
        try:
            if not tax_path.exists():
                print(f"WARNING: Taxonomy file not found: {tax_path}")
            else:
                tax_df = pd.read_csv(tax_path, sep='\t', index_col=0)
                tax_df.index.name = "OTU_ID"
                print("Taxonomy file loaded successfully.")
        except Exception as e:
            print(f"WARNING: Cannot read taxonomy file {tax_path}: {e}")
    else:
        print("No taxonomy file provided. Output will contain counts only.")

    if tax_df is not None:
        print("Joining taxonomy with counts (inner join)...")
        final_table_df = tax_df.join(counts_df, how='inner')
    else:
        final_table_df = counts_df
        final_table_df.index.name = "OTU_ID"

    print(f"Saving final table to {output_path}...")
    try:
        final_table_df.to_excel(output_path, engine='openpyxl')
        print(f"Final OTU table (with taxonomy) saved to: {output_path}")
    except ImportError:
        print("Module 'openpyxl' not found. Saving as CSV.")
        csv_output = output_path.with_suffix('.csv')
        final_table_df.to_csv(csv_output)
        print(f"Final OTU table saved to: {csv_output}")
    except Exception as e:
        print(f"Error saving final file {output_path}: {e}")
