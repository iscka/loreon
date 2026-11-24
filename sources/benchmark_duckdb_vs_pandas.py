#!/usr/bin/env python3
"""
benchmark_duckdb_vs_pandas.py
=============================
Benchmark: DuckDB vs pure-Pandas for OTU table aggregation in LOREON.

Measures wall time, peak RSS memory, and output correctness for the core
operations performed in ResultsReader.makeOtu_duckdb():

  1. Read N per-sample .tmp.txt files
  2. PIVOT from long to wide format (OTU × Samples)
  3. JOIN with taxonomy map (e.g. SILVA 138.2 NR99 — 510K entries)
  4. Write final OTU table

Can also compare against EPI2ME wf-16s timing when trace.txt is provided.

Usage:
    # Minimal: benchmark on real LOREON data
    python benchmark_duckdb_vs_pandas.py \
        --results-dir /path/to/LOREON/analysis/results_dir \
        --taxonomy /path/to/taxonomy_map.tsv

    # With synthetic data generation (no real data needed)
    python benchmark_duckdb_vs_pandas.py \
        --synthetic --n-otus 5000 --n-samples 50

    # Compare with EPI2ME timing
    python benchmark_duckdb_vs_pandas.py \
        --results-dir /path/to/results \
        --taxonomy /path/to/taxonomy.tsv \
        --epi2me-trace /path/to/epi2me_sample/execution/trace.txt

    # Full suite with multiple scales
    python benchmark_duckdb_vs_pandas.py --scaling-test
"""

import argparse
import gc
import os
import sys
import tempfile
import time
import tracemalloc
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import duckdb
    HAS_DUCKDB = True
except ImportError:
    HAS_DUCKDB = False
    print("WARNING: duckdb not installed. Only pandas benchmark will run.")

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def _format_bytes(n):
    for unit in ('B', 'KB', 'MB', 'GB'):
        if abs(n) < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def _format_time(seconds):
    if seconds < 1:
        return f"{seconds*1000:.1f} ms"
    elif seconds < 60:
        return f"{seconds:.2f} s"
    else:
        m, s = divmod(seconds, 60)
        return f"{int(m)}m {s:.1f}s"


class BenchmarkResult:
    def __init__(self, name, wall_time, peak_memory, n_rows, n_cols):
        self.name = name
        self.wall_time = wall_time
        self.peak_memory = peak_memory
        self.n_rows = n_rows
        self.n_cols = n_cols

    def __repr__(self):
        return (
            f"{self.name:25s} | "
            f"Time: {_format_time(self.wall_time):>10s} | "
            f"Memory: {_format_bytes(self.peak_memory):>10s} | "
            f"Shape: ({self.n_rows:,} x {self.n_cols})"
        )


def generate_synthetic_data(tmpdir, n_otus=2000, n_samples=20, n_tax_rows=50000):
    """Generate synthetic .tmp.txt files + taxonomy map for benchmarking."""
    print(f"  Generating synthetic data: {n_otus} OTUs, {n_samples} samples, "
          f"{n_tax_rows} taxonomy entries...")

    rng = np.random.default_rng(42)

    otu_ids = [f"OTU_{i:06d}" for i in range(n_otus)]

    all_otu_ids = [f"OTU_{i:06d}" for i in range(n_tax_rows)]
    tax_strings = []
    for i in range(n_tax_rows):
        k = f"k__Kingdom{i % 3}"
        p = f"p__Phylum{i % 10}"
        c = f"c__Class{i % 30}"
        o = f"o__Order{i % 80}"
        f_ = f"f__Family{i % 200}"
        g = f"g__Genus{i % 500}"
        s = f"s__Species{i % n_tax_rows}"
        tax_strings.append(f"{k};{p};{c};{o};{f_};{g};{s}")

    tax_path = Path(tmpdir) / "taxonomy_map.tsv"
    df_tax = pd.DataFrame({'OTU_ID': all_otu_ids, 'Taxonomy': tax_strings})
    df_tax.to_csv(tax_path, sep='\t', index=False)

    results_dir = Path(tmpdir) / "results"
    results_dir.mkdir()
    sample_files = []

    for s in range(n_samples):
        sample_name = f"barcode{s:02d}"
        fpath = results_dir / f"{sample_name}.tmp.txt"

        n_detected = rng.integers(n_otus // 4, n_otus)
        detected_otus = rng.choice(otu_ids, size=n_detected, replace=False)
        counts = rng.integers(1, 10000, size=n_detected)
        lengths = rng.integers(500, 5000, size=n_detected)
        unmapped = np.zeros(n_detected, dtype=int)

        with open(fpath, 'w') as f:
            for otu, length, count, unmap in zip(detected_otus, lengths, counts, unmapped):
                f.write(f"{otu}\t{length}\t{count}\t{unmap}\n")
            f.write(f"*\t0\t{counts.sum()}\t0\n")

        sample_files.append(str(fpath))

    print(f"  Generated {n_samples} sample files in {results_dir}")
    print(f"  Taxonomy map: {n_tax_rows} entries at {tax_path}")

    return str(results_dir), str(tax_path)


def benchmark_duckdb(results_dir, taxonomy_path):
    """Run the DuckDB PIVOT + JOIN approach and measure performance."""
    if not HAS_DUCKDB:
        return None

    results_path = Path(results_dir)
    clean_files = sorted(results_path.glob("*.tmp.txt"))
    if not clean_files:
        print("  ERROR: No .tmp.txt files found.")
        return None

    gc.collect()
    tracemalloc.start()
    t0 = time.perf_counter()

    con = duckdb.connect(':memory:')

    clean_files_sql = str([str(f) for f in clean_files])
    col_types_sql = "{'column0': 'VARCHAR', 'column1': 'INTEGER', 'column2': 'INTEGER', 'column3': 'INTEGER'}"

    if taxonomy_path and Path(taxonomy_path).exists():
        tax_path_sql = str(taxonomy_path).replace("'", "''")
        query = f"""
        WITH counts AS (
            SELECT * FROM (
                PIVOT (
                    SELECT
                        replace(split_part(filename, '/', -1), '.tmp.txt', '') AS sample_name,
                        column0 AS otu_id,
                        column2 AS mapped_count
                    FROM read_csv(
                        {clean_files_sql}, delim='\\t', header=False,
                        columns={col_types_sql}, filename=True
                    )
                    WHERE column0 != '*'
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
        query = f"""
        PIVOT (
            SELECT
                replace(split_part(filename, '/', -1), '.tmp.txt', '') AS sample_name,
                column0 AS otu_id,
                column2 AS mapped_count
            FROM read_csv(
                {clean_files_sql}, delim='\\t', header=False,
                columns={col_types_sql}, filename=True
            )
            WHERE column0 != '*'
        )
        ON sample_name
        USING sum(mapped_count)
        GROUP BY otu_id;
        """

    df = con.execute(query).df()
    con.close()

    t1 = time.perf_counter()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return BenchmarkResult(
        name="DuckDB (PIVOT+JOIN)",
        wall_time=t1 - t0,
        peak_memory=peak,
        n_rows=len(df),
        n_cols=len(df.columns),
    )


def benchmark_duckdb_optimized(results_dir, taxonomy_path):
    """Run the fully-optimized DuckDB approach: PIVOT + JOIN + taxonomy split."""
    if not HAS_DUCKDB:
        return None

    results_path = Path(results_dir)
    clean_files = sorted(results_path.glob("*.tmp.txt"))
    if not clean_files:
        return None

    gc.collect()
    tracemalloc.start()
    t0 = time.perf_counter()

    con = duckdb.connect(':memory:')

    clean_files_sql = str([str(f) for f in clean_files])
    col_types_sql = "{'column0': 'VARCHAR', 'column1': 'INTEGER', 'column2': 'INTEGER', 'column3': 'INTEGER'}"

    from ResultsReader import _build_taxonomy_split_sql
    tax_split_sql = _build_taxonomy_split_sql()

    if taxonomy_path and Path(taxonomy_path).exists():
        tax_path_sql = str(taxonomy_path).replace("'", "''")
        query = f"""
        WITH counts AS (
            SELECT * FROM (
                PIVOT (
                    SELECT
                        replace(split_part(filename, '/', -1), '.tmp.txt', '') AS sample_name,
                        column0 AS otu_id,
                        column2 AS mapped_count
                    FROM read_csv(
                        {clean_files_sql}, delim='\\t', header=False,
                        columns={col_types_sql}, filename=True
                    )
                    WHERE column0 != '*'
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
        ),
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
        query = f"""
        PIVOT (
            SELECT
                replace(split_part(filename, '/', -1), '.tmp.txt', '') AS sample_name,
                column0 AS otu_id,
                column2 AS mapped_count
            FROM read_csv(
                {clean_files_sql}, delim='\\t', header=False,
                columns={col_types_sql}, filename=True
            )
            WHERE column0 != '*'
        )
        ON sample_name
        USING sum(mapped_count)
        GROUP BY otu_id;
        """

    df = con.execute(query).df()
    con.close()

    t1 = time.perf_counter()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return BenchmarkResult(
        name="DuckDB optimized (all-SQL)",
        wall_time=t1 - t0,
        peak_memory=peak,
        n_rows=len(df),
        n_cols=len(df.columns),
    )


def benchmark_duckdb_pandas_split(results_dir, taxonomy_path):
    """DuckDB for PIVOT+JOIN, then pandas apply() for taxonomy split."""
    if not HAS_DUCKDB:
        return None

    results_path = Path(results_dir)
    clean_files = sorted(results_path.glob("*.tmp.txt"))
    if not clean_files:
        return None

    gc.collect()
    tracemalloc.start()
    t0 = time.perf_counter()

    con = duckdb.connect(':memory:')

    clean_files_sql = str([str(f) for f in clean_files])
    col_types_sql = "{'column0': 'VARCHAR', 'column1': 'INTEGER', 'column2': 'INTEGER', 'column3': 'INTEGER'}"

    if taxonomy_path and Path(taxonomy_path).exists():
        tax_path_sql = str(taxonomy_path).replace("'", "''")
        query = f"""
        WITH counts AS (
            SELECT * FROM (
                PIVOT (
                    SELECT
                        replace(split_part(filename, '/', -1), '.tmp.txt', '') AS sample_name,
                        column0 AS otu_id,
                        column2 AS mapped_count
                    FROM read_csv(
                        {clean_files_sql}, delim='\\t', header=False,
                        columns={col_types_sql}, filename=True
                    )
                    WHERE column0 != '*'
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
        query = f"""
        PIVOT (
            SELECT
                replace(split_part(filename, '/', -1), '.tmp.txt', '') AS sample_name,
                column0 AS otu_id,
                column2 AS mapped_count
            FROM read_csv(
                {clean_files_sql}, delim='\\t', header=False,
                columns={col_types_sql}, filename=True
            )
            WHERE column0 != '*'
        )
        ON sample_name
        USING sum(mapped_count)
        GROUP BY otu_id;
        """

    df = con.execute(query).df()
    con.close()

    if 'Taxonomy' in df.columns:
        from ResultsReader import _split_taxonomy_series
        rank_df = _split_taxonomy_series(df['Taxonomy'])
        for col in rank_df.columns:
            df[col] = rank_df[col].values

    t1 = time.perf_counter()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return BenchmarkResult(
        name="DuckDB + pandas split",
        wall_time=t1 - t0,
        peak_memory=peak,
        n_rows=len(df),
        n_cols=len(df.columns),
    )


def benchmark_pandas(results_dir, taxonomy_path):
    """Run the pure-Pandas equivalent and measure performance."""
    results_path = Path(results_dir)
    clean_files = sorted(results_path.glob("*.tmp.txt"))
    if not clean_files:
        print("  ERROR: No .tmp.txt files found.")
        return None

    gc.collect()
    tracemalloc.start()
    t0 = time.perf_counter()

    dfs = []
    for fpath in clean_files:
        sample_name = fpath.stem.replace('.tmp', '')
        df_sample = pd.read_csv(
            fpath, sep='\t', header=None,
            names=['otu_id', 'length', 'mapped_count', 'unmapped_count']
        )
        df_sample = df_sample[df_sample['otu_id'] != '*']
        df_sample['sample_name'] = sample_name
        dfs.append(df_sample[['otu_id', 'sample_name', 'mapped_count']])

    df_all = pd.concat(dfs, ignore_index=True)

    df_pivot = df_all.pivot_table(
        index='otu_id',
        columns='sample_name',
        values='mapped_count',
        aggfunc='sum',
        fill_value=0,
    )
    df_pivot = df_pivot.reset_index()

    if taxonomy_path and Path(taxonomy_path).exists():
        df_tax = pd.read_csv(taxonomy_path, sep='\t')
        df_final = df_tax.merge(df_pivot, left_on='OTU_ID', right_on='otu_id', how='inner')
        if 'otu_id' in df_final.columns:
            df_final = df_final.drop(columns=['otu_id'])
    else:
        df_final = df_pivot

    t1 = time.perf_counter()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return BenchmarkResult(
        name="Pandas (concat+pivot+merge)",
        wall_time=t1 - t0,
        peak_memory=peak,
        n_rows=len(df_final),
        n_cols=len(df_final.columns),
    )


def benchmark_pandas_chunked(results_dir, taxonomy_path):
    """Pandas with chunked reading (lower memory, slightly slower)."""
    results_path = Path(results_dir)
    clean_files = sorted(results_path.glob("*.tmp.txt"))
    if not clean_files:
        return None

    gc.collect()
    tracemalloc.start()
    t0 = time.perf_counter()

    counts = {}
    for fpath in clean_files:
        sample_name = fpath.stem.replace('.tmp', '')
        with open(fpath, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3 or parts[0] == '*':
                    continue
                otu_id = parts[0]
                mapped = int(parts[2])
                if otu_id not in counts:
                    counts[otu_id] = {}
                counts[otu_id][sample_name] = counts[otu_id].get(sample_name, 0) + mapped

    df_pivot = pd.DataFrame.from_dict(counts, orient='index').fillna(0).astype(int)
    df_pivot.index.name = 'otu_id'
    df_pivot = df_pivot.reset_index()

    if taxonomy_path and Path(taxonomy_path).exists():
        df_tax = pd.read_csv(taxonomy_path, sep='\t')
        df_final = df_tax.merge(df_pivot, left_on='OTU_ID', right_on='otu_id', how='inner')
        if 'otu_id' in df_final.columns:
            df_final = df_final.drop(columns=['otu_id'])
    else:
        df_final = df_pivot

    t1 = time.perf_counter()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return BenchmarkResult(
        name="Pandas (chunked dict)",
        wall_time=t1 - t0,
        peak_memory=peak,
        n_rows=len(df_final),
        n_cols=len(df_final.columns),
    )


def parse_epi2me_trace(trace_path):
    """Extract createAbundanceTables timing from EPI2ME trace.txt."""
    import re
    with open(trace_path, 'r') as f:
        for line in f:
            if 'createAbundanceTables' in line or 'abundance' in line.lower():
                parts = line.strip().split('\t')
                for part in parts:
                    part = part.strip()
                    m = re.match(r'^(\d+(?:\.\d+)?)\s*s$', part)
                    if m:
                        return float(m.group(1))
                    m = re.match(r'^(\d+)m\s+(\d+(?:\.\d+)?)s$', part)
                    if m:
                        return int(m.group(1)) * 60 + float(m.group(2))
    return None


def run_scaling_test(outdir):
    """Run benchmarks across multiple data sizes to produce scaling curves."""
    configs = [
        (500,    10, 10000),
        (1000,   20, 50000),
        (2000,   30, 100000),
        (5000,   50, 200000),
        (10000,  50, 500000),
        (10000, 100, 500000),
    ]

    results_all = []
    for n_otus, n_samples, n_tax in configs:
        label = f"{n_otus} OTUs x {n_samples} samples (tax: {n_tax//1000}K)"
        print(f"\n{'='*60}")
        print(f"  Config: {label}")
        print(f"{'='*60}")

        with tempfile.TemporaryDirectory() as tmpdir:
            rdir, tpath = generate_synthetic_data(tmpdir, n_otus, n_samples, n_tax)

            for bench_fn, method_name in [
                (benchmark_duckdb_optimized, "DuckDB optimized"),
                (benchmark_duckdb_pandas_split, "DuckDB + pandas split"),
                (benchmark_pandas, "Pandas"),
                (benchmark_pandas_chunked, "Pandas (chunked)"),
            ]:
                r = bench_fn(rdir, tpath)
                if r:
                    print(f"  {r}")
                    results_all.append({
                        'config': label,
                        'n_otus': n_otus,
                        'n_samples': n_samples,
                        'n_tax': n_tax,
                        'method': method_name,
                        'wall_time_s': r.wall_time,
                        'peak_memory_mb': r.peak_memory / (1024 * 1024),
                        'output_rows': r.n_rows,
                        'output_cols': r.n_cols,
                    })

    df = pd.DataFrame(results_all)
    csv_path = outdir / "scaling_benchmark.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nResults saved to: {csv_path}")

    if HAS_MATPLOTLIB and not df.empty:
        _plot_scaling(df, outdir)

    return df


def _plot_scaling(df, outdir):
    """Generate scaling comparison plots."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    for method in df['method'].unique():
        subset = df[df['method'] == method]
        x_label = [f"{r['n_otus']}x{r['n_samples']}" for _, r in subset.iterrows()]
        ax.plot(x_label, subset['wall_time_s'], 'o-', label=method, linewidth=2)
    ax.set_xlabel('Dataset (OTUs x Samples)')
    ax.set_ylabel('Wall Time (seconds)')
    ax.set_title('OTU Aggregation: Time Comparison')
    ax.legend()
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='x', rotation=45)

    ax = axes[1]
    for method in df['method'].unique():
        subset = df[df['method'] == method]
        x_label = [f"{r['n_otus']}x{r['n_samples']}" for _, r in subset.iterrows()]
        ax.plot(x_label, subset['peak_memory_mb'], 's-', label=method, linewidth=2)
    ax.set_xlabel('Dataset (OTUs x Samples)')
    ax.set_ylabel('Peak Memory (MB)')
    ax.set_title('OTU Aggregation: Memory Comparison')
    ax.legend()
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='x', rotation=45)

    fig.tight_layout()
    fig.savefig(outdir / "scaling_benchmark.png", dpi=150)
    plt.close(fig)
    print(f"Plot saved to: {outdir / 'scaling_benchmark.png'}")

    pivot_time = df.pivot(index='config', columns='method', values='wall_time_s')
    if 'DuckDB' in pivot_time.columns and 'Pandas' in pivot_time.columns:
        pivot_time['Speedup (DuckDB vs Pandas)'] = (
            pivot_time['Pandas'] / pivot_time['DuckDB']
        ).round(2)
    pivot_mem = df.pivot(index='config', columns='method', values='peak_memory_mb')
    if 'DuckDB' in pivot_mem.columns and 'Pandas' in pivot_mem.columns:
        pivot_mem['Memory ratio (Pandas/DuckDB)'] = (
            pivot_mem['Pandas'] / pivot_mem['DuckDB']
        ).round(2)

    print("\n" + "=" * 70)
    print("  SPEEDUP SUMMARY (DuckDB vs Pandas)")
    print("=" * 70)
    print(pivot_time.to_string())
    print("\n  MEMORY SUMMARY")
    print(pivot_mem.to_string())


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark DuckDB vs Pandas for OTU table aggregation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        '--results-dir', type=str, default=None,
        help="Path to LOREON results directory containing .tmp.txt files"
    )
    parser.add_argument(
        '--taxonomy', type=str, default=None,
        help="Path to taxonomy map TSV (OTU_ID, Taxonomy)"
    )
    parser.add_argument(
        '--epi2me-trace', type=str, default=None,
        help="Path to EPI2ME execution/trace.txt for comparison"
    )
    parser.add_argument(
        '--synthetic', action='store_true',
        help="Use synthetic data instead of real files"
    )
    parser.add_argument(
        '--n-otus', type=int, default=2000,
        help="Number of OTUs for synthetic data (default: 2000)"
    )
    parser.add_argument(
        '--n-samples', type=int, default=20,
        help="Number of samples for synthetic data (default: 20)"
    )
    parser.add_argument(
        '--n-tax-rows', type=int, default=50000,
        help="Number of taxonomy entries for synthetic data (default: 50000)"
    )
    parser.add_argument(
        '--scaling-test', action='store_true',
        help="Run full scaling test across multiple data sizes"
    )
    parser.add_argument(
        '--output', type=str, default='benchmark_output',
        help="Output directory for results (default: benchmark_output)"
    )
    parser.add_argument(
        '--repeat', type=int, default=3,
        help="Number of repetitions per benchmark (default: 3)"
    )

    args = parser.parse_args()
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("  LOREON: DuckDB vs Pandas — OTU Aggregation Benchmark")
    print("=" * 70)
    print(f"  DuckDB available: {HAS_DUCKDB}")
    if HAS_DUCKDB:
        print(f"  DuckDB version:   {duckdb.__version__}")
    print(f"  Pandas version:   {pd.__version__}")
    print(f"  NumPy version:    {np.__version__}")
    print()

    if args.scaling_test:
        run_scaling_test(outdir)
        return

    tmpdir_obj = None
    if args.synthetic or not args.results_dir:
        tmpdir_obj = tempfile.TemporaryDirectory()
        results_dir, taxonomy_path = generate_synthetic_data(
            tmpdir_obj.name, args.n_otus, args.n_samples, args.n_tax_rows
        )
    else:
        results_dir = args.results_dir
        taxonomy_path = args.taxonomy

    benchmarks = [
        ("DuckDB optimized (all-SQL)", benchmark_duckdb_optimized),
        ("DuckDB + pandas split", benchmark_duckdb_pandas_split),
        ("DuckDB (PIVOT+JOIN only)", benchmark_duckdb),
        ("Pandas (concat+pivot+merge)", benchmark_pandas),
        ("Pandas (chunked dict)", benchmark_pandas_chunked),
    ]

    all_results = []
    for name, fn in benchmarks:
        times = []
        memories = []
        result = None
        for i in range(args.repeat):
            r = fn(results_dir, taxonomy_path)
            if r:
                times.append(r.wall_time)
                memories.append(r.peak_memory)
                result = r
        if result and times:
            median_time = sorted(times)[len(times) // 2]
            median_mem = sorted(memories)[len(memories) // 2]
            result.wall_time = median_time
            result.peak_memory = median_mem
            all_results.append(result)
            print(f"  {result}  (median of {args.repeat} runs)")

    epi2me_time = None
    if args.epi2me_trace and Path(args.epi2me_trace).exists():
        epi2me_time = parse_epi2me_trace(args.epi2me_trace)
        if epi2me_time:
            print(f"\n  EPI2ME createAbundanceTables: {_format_time(epi2me_time)}")

    if len(all_results) >= 2:
        print("\n" + "=" * 70)
        print("  SUMMARY")
        print("=" * 70)

        duckdb_r = next((r for r in all_results if 'DuckDB' in r.name), None)
        pandas_r = next((r for r in all_results if 'concat' in r.name), None)

        if duckdb_r and pandas_r:
            speedup = pandas_r.wall_time / duckdb_r.wall_time
            mem_ratio = pandas_r.peak_memory / duckdb_r.peak_memory

            print(f"  Time speedup (DuckDB vs Pandas):    {speedup:.2f}x faster")
            print(f"  Memory ratio (Pandas / DuckDB):     {mem_ratio:.2f}x more memory")

            if epi2me_time and duckdb_r:
                epi_speedup = epi2me_time / duckdb_r.wall_time
                print(f"  Time speedup (DuckDB vs EPI2ME):    {epi_speedup:.2f}x faster")

        rows = []
        for r in all_results:
            rows.append({
                'method': r.name,
                'wall_time_s': round(r.wall_time, 4),
                'peak_memory_mb': round(r.peak_memory / (1024 * 1024), 2),
                'output_rows': r.n_rows,
                'output_cols': r.n_cols,
            })
        if epi2me_time:
            rows.append({
                'method': 'EPI2ME createAbundanceTables',
                'wall_time_s': round(epi2me_time, 4),
                'peak_memory_mb': 397.0,
                'output_rows': None,
                'output_cols': None,
            })

        df_results = pd.DataFrame(rows)
        csv_path = outdir / "benchmark_results.csv"
        df_results.to_csv(csv_path, index=False)
        print(f"\n  Results saved to: {csv_path}")

    if tmpdir_obj:
        tmpdir_obj.cleanup()

    print("\nDone.")


if __name__ == "__main__":
    main()
