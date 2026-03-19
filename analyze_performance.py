#!/usr/bin/env python3
"""
LOREON Performance Analyzer
============================
Scans all analysis directories under a root folder, collects
performance_report.json, mapping_stats.json, filter reports, and
pipeline.log to build a comprehensive benchmark comparison.

Output: an interactive HTML report with Plotly box-plots and tables.

Usage:
    python3 analyze_performance.py [--analyses-dir /path/to/ANALISI] [--output report.html]
"""

import argparse
import json
import os
import platform
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio


# ── System info ─────────────────────────────────────────────────────────────

def get_system_info():
    """Collect host machine specs (works on macOS, Linux, Windows)."""
    info = {
        'os': f"{platform.system()} {platform.release()}",
        'arch': platform.machine(),
        'python': platform.python_version(),
    }

    # CPU
    try:
        if platform.system() == 'Darwin':
            info['cpu'] = subprocess.check_output(
                ['sysctl', '-n', 'machdep.cpu.brand_string'],
                text=True
            ).strip()
            info['cores_physical'] = int(subprocess.check_output(
                ['sysctl', '-n', 'hw.physicalcpu'], text=True
            ).strip())
            info['cores_logical'] = int(subprocess.check_output(
                ['sysctl', '-n', 'hw.logicalcpu'], text=True
            ).strip())
            ram_bytes = int(subprocess.check_output(
                ['sysctl', '-n', 'hw.memsize'], text=True
            ).strip())
            info['ram_gb'] = round(ram_bytes / (1024 ** 3), 1)
        elif platform.system() == 'Linux':
            with open('/proc/cpuinfo') as f:
                for line in f:
                    if line.startswith('model name'):
                        info['cpu'] = line.split(':')[1].strip()
                        break
            info['cores_logical'] = os.cpu_count() or 0
            info['cores_physical'] = info['cores_logical']
            with open('/proc/meminfo') as f:
                for line in f:
                    if line.startswith('MemTotal'):
                        kb = int(line.split()[1])
                        info['ram_gb'] = round(kb / (1024 ** 2), 1)
                        break
        else:
            info['cpu'] = platform.processor() or 'Unknown'
            info['cores_logical'] = os.cpu_count() or 0
            info['cores_physical'] = info['cores_logical']
            info['ram_gb'] = 0
    except Exception:
        info.setdefault('cpu', 'Unknown')
        info.setdefault('cores_logical', os.cpu_count() or 0)
        info.setdefault('cores_physical', info.get('cores_logical', 0))
        info.setdefault('ram_gb', 0)

    return info


# ── Data collection ─────────────────────────────────────────────────────────

def parse_command_from_log(log_path):
    """Extract pipeline CLI arguments from pipeline.log."""
    params = {}
    try:
        with open(log_path, 'r', encoding='utf-8') as f:
            for line in f:
                if 'Command:' in line:
                    cmd = line.split('Command:')[1].strip()
                    # Database path
                    m = re.search(r'-d\s+(\S+)', cmd)
                    if m:
                        params['db_path'] = m.group(1)
                        params['db_name'] = Path(m.group(1)).stem
                    # Format
                    m = re.search(r'-f\s+(\S+)', cmd)
                    if m:
                        params['format'] = m.group(1)
                    # Threads
                    m = re.search(r'-T\s+(\d+)', cmd)
                    if m:
                        params['total_threads'] = int(m.group(1))
                    m = re.search(r'-t\s+(\d+)', cmd)
                    if m:
                        params['threads_per_job'] = int(m.group(1))
                    # Kmer / window
                    m = re.search(r'-k\s+(\d+)', cmd)
                    if m:
                        params['kmer_size'] = int(m.group(1))
                    m = re.search(r'-w\s+(\d+)', cmd)
                    if m:
                        params['window_size'] = int(m.group(1))
                    break
    except Exception:
        pass
    return params


def collect_analysis(analysis_dir):
    """Collect all data for a single analysis directory."""
    analysis_dir = Path(analysis_dir)
    name = analysis_dir.name
    record = {'analysis_name': name}

    # Performance report
    perf_path = analysis_dir / 'performance_report.json'
    if not perf_path.exists():
        return None
    with open(perf_path, 'r', encoding='utf-8') as f:
        perf = json.load(f)

    record['run_date'] = perf.get('run_date', '')
    record['total_wall_time_s'] = perf.get('total_wall_time_s', 0)
    record['peak_memory_mb'] = perf.get('peak_memory_mb', 0)

    # Per-step times
    for step_name, step_data in perf.get('steps', {}).items():
        key = step_name.replace('Step ', 'step').replace(': ', '_').replace(' ', '_').lower()
        record[f'{key}_time_s'] = step_data.get('wall_time_s', 0)
        record[f'{key}_mem_mb'] = step_data.get('mem_mb', 0)

        # Extra metrics
        if 'taxonomy_entries' in step_data:
            record['db_entries'] = step_data['taxonomy_entries']
        if 'sequences_total' in step_data:
            record['sequences_total'] = step_data['sequences_total']
        if 'sequences_good' in step_data:
            record['sequences_good'] = step_data['sequences_good']
        if 'throughput_seq_s' in step_data:
            record['filter_throughput_seq_s'] = step_data['throughput_seq_s']
        if 'files_to_map' in step_data:
            record['n_samples'] = step_data['files_to_map']

    # Mapping stats
    ms_path = analysis_dir / 'mapping_stats.json'
    if ms_path.exists():
        with open(ms_path, 'r', encoding='utf-8') as f:
            ms = json.load(f)
        summary = ms.get('summary', {})
        record['mapped_reads'] = summary.get('total_mapped_all', 0)
        record['mapping_rate_pct'] = summary.get('overall_mapping_rate_pct', 0)

    # Filter report (for sequence lengths)
    filter_reports = list(analysis_dir.glob('report_filtering_*.tsv'))
    if filter_reports:
        fr_path = filter_reports[0]
        # Extract min/max from filename
        m = re.search(r'report_filtering_(\d+)_(\d+)', fr_path.name)
        if m:
            record['filter_min_len'] = int(m.group(1))
            record['filter_max_len'] = int(m.group(2))
            is_disabled = record['filter_min_len'] == 0 and record['filter_max_len'] >= 999999
            record['len_filter_active'] = not is_disabled

        try:
            df_filter = pd.read_csv(fr_path, sep='\t')
            record['mean_seq_length_bp'] = round(df_filter['mean length'].mean(), 1)
            record['median_seq_length_bp'] = round(df_filter['median length'].mean(), 1)
        except Exception:
            pass

    # Pipeline log (for DB name, format, threads)
    log_path = analysis_dir / 'pipeline.log'
    if log_path.exists():
        log_params = parse_command_from_log(log_path)
        record.update(log_params)

    # Per-sample mapping stats for box plots
    if ms_path.exists():
        samples = ms.get('samples', [])
        record['_sample_details'] = samples

    return record


def collect_all(analyses_root):
    """Scan all subdirectories and collect analysis records."""
    analyses_root = Path(analyses_root)
    records = []
    for child in sorted(analyses_root.iterdir()):
        if child.is_dir():
            rec = collect_analysis(child)
            if rec is not None:
                records.append(rec)
    return records


# ── HTML Report Generation ──────────────────────────────────────────────────

def build_report(records, sys_info, output_path):
    """Generate the full HTML benchmark report with Plotly charts."""
    df = pd.DataFrame(records)

    # Drop internal fields for display
    df_display = df.drop(columns=['_sample_details'], errors='ignore')

    # ── Prepare per-sample dataframe for box plots ──
    sample_rows = []
    for rec in records:
        for s in rec.get('_sample_details', []):
            sample_rows.append({
                'analysis': rec['analysis_name'],
                'sample': s.get('sample', ''),
                'total_reads': s.get('total_reads', 0),
                'mapped_unfiltered': s.get('mapped_unfiltered', 0),
                'mapping_rate_pct': s.get('mapping_rate_pct', 0),
            })
    df_samples = pd.DataFrame(sample_rows) if sample_rows else pd.DataFrame()

    # ── Prepare per-barcode filter stats for box plots ──
    barcode_rows = []
    for rec in records:
        analysis_dir = Path(analyses_root) / rec['analysis_name']
        filter_reports = list(analysis_dir.glob('report_filtering_*.tsv'))
        if filter_reports:
            try:
                df_fr = pd.read_csv(filter_reports[0], sep='\t')
                for _, row in df_fr.iterrows():
                    barcode_rows.append({
                        'analysis': rec['analysis_name'],
                        'barcode': row.get('barcode', ''),
                        'mean_length': row.get('mean length', 0),
                        'median_length': row.get('median length', 0),
                        'good_count': row.get('good count', 0),
                        'bad_count': row.get('bad count', 0),
                        'total_count': row.get('total count', 0),
                    })
            except Exception:
                pass
    df_barcodes = pd.DataFrame(barcode_rows) if barcode_rows else pd.DataFrame()

    # ═══════════════════════════════════════════════════════════════════════
    # CHARTS
    # ═══════════════════════════════════════════════════════════════════════

    charts_html = []

    # ── 1. Step times stacked bar ──
    step_cols = [c for c in df.columns if c.endswith('_time_s') and c.startswith('step')]
    if step_cols:
        step_labels = {
            'step1_taxonomy_time_s': 'Taxonomy Parsing',
            'step2_filter_time_s': 'Sequence Filter',
            'step3_filter_report_time_s': 'Filter Report',
            'step4_mapping_time_s': 'Mapping & Tabeling',
            'step5_reformat_time_s': 'Index Reformat',
            'step6_otu_aggregation_time_s': 'OTU Aggregation',
        }
        fig_steps = go.Figure()
        for col in step_cols:
            label = step_labels.get(col, col)
            fig_steps.add_trace(go.Bar(
                name=label,
                x=df['analysis_name'],
                y=df[col],
                text=df[col].apply(lambda v: f'{v:.1f}s' if v < 100 else f'{v:.0f}s'),
                textposition='inside',
            ))
        fig_steps.update_layout(
            barmode='stack',
            title='Pipeline Step Execution Time (Stacked)',
            xaxis_title='Analysis',
            yaxis_title='Time (seconds)',
            yaxis_type='log',
            height=500,
            legend=dict(orientation='h', yanchor='bottom', y=1.02),
        )
        charts_html.append(pio.to_html(fig_steps, full_html=False, include_plotlyjs='cdn'))

    # ── 2. Box-plot: sequence lengths per barcode ──
    if not df_barcodes.empty:
        fig_len = px.box(
            df_barcodes, x='analysis', y='mean_length',
            color='analysis', points='all',
            title='Distribution of Mean Sequence Length per Barcode',
            labels={'mean_length': 'Mean Length (bp)', 'analysis': 'Analysis'},
        )
        fig_len.update_layout(height=450, showlegend=False)
        charts_html.append(pio.to_html(fig_len, full_html=False, include_plotlyjs=False))

    # ── 3. Box-plot: reads per sample ──
    if not df_samples.empty:
        fig_reads = px.box(
            df_samples, x='analysis', y='total_reads',
            color='analysis', points='all',
            title='Distribution of Reads per Sample',
            labels={'total_reads': 'Reads per Sample', 'analysis': 'Analysis'},
        )
        fig_reads.update_layout(height=450, showlegend=False)
        charts_html.append(pio.to_html(fig_reads, full_html=False, include_plotlyjs=False))

    # ── 4. Box-plot: mapping rate per sample ──
    if not df_samples.empty:
        fig_maprate = px.box(
            df_samples, x='analysis', y='mapping_rate_pct',
            color='analysis', points='all',
            title='Mapping Rate Distribution per Sample',
            labels={'mapping_rate_pct': 'Mapping Rate (%)', 'analysis': 'Analysis'},
        )
        fig_maprate.update_layout(height=450, showlegend=False, yaxis_range=[90, 101])
        charts_html.append(pio.to_html(fig_maprate, full_html=False, include_plotlyjs=False))

    # ── 5. Box-plot: good/bad sequences per barcode ──
    if not df_barcodes.empty:
        df_gb = df_barcodes.melt(
            id_vars=['analysis', 'barcode'],
            value_vars=['good_count', 'bad_count'],
            var_name='type', value_name='count'
        )
        df_gb['type'] = df_gb['type'].replace({'good_count': 'Passed', 'bad_count': 'Rejected'})
        fig_gb = px.box(
            df_gb, x='analysis', y='count', color='type',
            points='all',
            title='Passed vs Rejected Sequences per Barcode',
            labels={'count': 'Sequences', 'analysis': 'Analysis'},
            color_discrete_map={'Passed': '#2ca02c', 'Rejected': '#d62728'},
        )
        fig_gb.update_layout(height=450)
        charts_html.append(pio.to_html(fig_gb, full_html=False, include_plotlyjs=False))

    # ── 6. Scatter: mapping time vs input reads (colored by DB size) ──
    if 'sequences_good' in df.columns and 'step4_mapping_time_s' in df.columns:
        fig_scatter = px.scatter(
            df, x='sequences_good', y='step4_mapping_time_s',
            size='db_entries', color='analysis_name',
            hover_data=['db_name', 'n_samples', 'mapped_reads'],
            title='Mapping Time vs Input Reads (bubble size = DB entries)',
            labels={
                'sequences_good': 'Input Reads (after filter)',
                'step4_mapping_time_s': 'Mapping Time (s)',
            },
        )
        fig_scatter.update_layout(height=450)
        charts_html.append(pio.to_html(fig_scatter, full_html=False, include_plotlyjs=False))

    # ── 7. Memory usage per step ──
    mem_cols = [c for c in df.columns if c.endswith('_mem_mb') and c.startswith('step')]
    if mem_cols:
        mem_labels = {
            'step1_taxonomy_mem_mb': 'Taxonomy',
            'step2_filter_mem_mb': 'Filter',
            'step3_filter_report_mem_mb': 'Filter Report',
            'step4_mapping_mem_mb': 'Mapping',
            'step5_reformat_mem_mb': 'Reformat',
            'step6_otu_aggregation_mem_mb': 'OTU Aggregation',
        }
        mem_rows = []
        for _, row in df.iterrows():
            for col in mem_cols:
                mem_rows.append({
                    'analysis': row['analysis_name'],
                    'step': mem_labels.get(col, col),
                    'mem_mb': row.get(col, 0),
                })
        df_mem = pd.DataFrame(mem_rows)
        fig_mem = px.box(
            df_mem, x='step', y='mem_mb', color='analysis',
            points='all',
            title='Memory Usage per Pipeline Step',
            labels={'mem_mb': 'RSS Memory (MB)', 'step': 'Step'},
        )
        fig_mem.update_layout(height=450)
        charts_html.append(pio.to_html(fig_mem, full_html=False, include_plotlyjs=False))

    # ── 8. Throughput: seq/s per analysis ──
    if 'filter_throughput_seq_s' in df.columns:
        fig_tp = px.bar(
            df, x='analysis_name', y='filter_throughput_seq_s',
            color='analysis_name',
            title='Filter Throughput (sequences/second)',
            labels={'filter_throughput_seq_s': 'Throughput (seq/s)', 'analysis_name': 'Analysis'},
            text='filter_throughput_seq_s',
        )
        fig_tp.update_traces(texttemplate='%{text:,.0f}', textposition='outside')
        fig_tp.update_layout(height=400, showlegend=False)
        charts_html.append(pio.to_html(fig_tp, full_html=False, include_plotlyjs=False))

    # ═══════════════════════════════════════════════════════════════════════
    # BUILD HTML
    # ═══════════════════════════════════════════════════════════════════════

    # Summary table
    summary_cols = [
        'analysis_name', 'run_date', 'db_name', 'format', 'db_entries',
        'sequences_total', 'sequences_good', 'mean_seq_length_bp',
        'n_samples', 'mapped_reads', 'mapping_rate_pct',
        'total_wall_time_s', 'peak_memory_mb',
        'filter_min_len', 'filter_max_len',
        'total_threads', 'threads_per_job',
    ]
    available_cols = [c for c in summary_cols if c in df.columns]
    df_summary = df[available_cols].copy()

    col_rename = {
        'analysis_name': 'Analysis',
        'run_date': 'Date',
        'db_name': 'Database',
        'format': 'Format',
        'db_entries': 'DB Entries',
        'sequences_total': 'Total Seq',
        'sequences_good': 'Good Seq',
        'mean_seq_length_bp': 'Mean Len (bp)',
        'n_samples': 'Samples',
        'mapped_reads': 'Mapped',
        'mapping_rate_pct': 'Map Rate (%)',
        'total_wall_time_s': 'Total Time (s)',
        'peak_memory_mb': 'Peak RAM (MB)',
        'filter_min_len': 'Min Len',
        'filter_max_len': 'Max Len',
        'total_threads': 'Threads',
        'threads_per_job': 'Threads/Job',
    }
    df_summary = df_summary.rename(columns=col_rename)

    summary_html = df_summary.to_html(
        classes='table table-striped table-hover table-sm',
        index=False, border=0
    )

    # System info card
    sys_card = f"""
    <div class="row mb-4">
        <div class="col-md-4">
            <div class="card shadow-sm h-100">
                <div class="card-header"><strong>Host Machine</strong></div>
                <div class="card-body">
                    <table class="table table-sm mb-0">
                        <tr><td>CPU</td><td><strong>{sys_info.get('cpu', 'N/A')}</strong></td></tr>
                        <tr><td>Cores</td><td>{sys_info.get('cores_physical', '?')} physical / {sys_info.get('cores_logical', '?')} logical</td></tr>
                        <tr><td>RAM</td><td>{sys_info.get('ram_gb', '?')} GB</td></tr>
                        <tr><td>OS</td><td>{sys_info.get('os', 'N/A')}</td></tr>
                        <tr><td>Architecture</td><td>{sys_info.get('arch', 'N/A')}</td></tr>
                        <tr><td>Python</td><td>{sys_info.get('python', 'N/A')}</td></tr>
                    </table>
                </div>
            </div>
        </div>
        <div class="col-md-8">
            <div class="card shadow-sm h-100">
                <div class="card-header"><strong>Analyses Found: {len(records)}</strong></div>
                <div class="card-body">
                    <table class="table table-sm mb-0">
                        <thead><tr><th>Analysis</th><th>Database</th><th>Format</th><th>Sequences</th><th>Time</th></tr></thead>
                        <tbody>
                        {"".join(f'<tr><td><strong>{r["analysis_name"]}</strong></td><td>{r.get("db_name","?")}</td><td>{r.get("format","?")}</td><td>{r.get("sequences_total",0):,} seq</td><td>{r.get("total_wall_time_s",0):.1f}s</td></tr>' for r in records)}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
    </div>
    """

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>LOREON Performance Benchmark</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {{ background-color: #f8f9fa; }}
        .container {{ max-width: 1400px; }}
        h1 {{ border-bottom: 2px solid #dee2e6; padding-bottom: 10px; margin-bottom: 20px; }}
        h2 {{ margin-top: 30px; color: #495057; }}
        .card {{ margin-bottom: 1rem; }}
        .table-wrapper {{ overflow-x: auto; }}
        .table th, .table td {{ text-align: left; }}
    </style>
</head>
<body>
    <div class="container mt-4">
        <h1>LOREON Performance Benchmark</h1>
        <p class="text-muted">Comparative analysis across {len(records)} pipeline runs</p>

        {sys_card}

        <h2>Analysis Summary</h2>
        <div class="card shadow-sm">
            <div class="card-body table-wrapper">
                {summary_html}
            </div>
        </div>

        <h2>Performance Charts</h2>
        {"".join(f'<div class="card shadow-sm"><div class="card-body">{c}</div></div>' for c in charts_html)}

    </div>
</body>
</html>"""

    output_path = Path(output_path)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)

    print(f"\nBenchmark report generated: {output_path.resolve()}")
    print(f"  Analyses: {len(records)}")
    print(f"  Charts: {len(charts_html)}")


# ── Main ────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='LOREON Performance Analyzer')
    parser.add_argument(
        '--analyses-dir', '-d',
        default='/Users/robm4/Documents/AIDocuments/dottorato/ANALISI',
        help='Root directory containing analysis folders'
    )
    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Output HTML report path (default: <analyses-dir>/performance_benchmark.html)'
    )
    args = parser.parse_args()

    analyses_root = Path(args.analyses_dir)
    if not analyses_root.exists():
        print(f"ERROR: Directory not found: {analyses_root}", file=sys.stderr)
        sys.exit(1)

    output_path = args.output or str(analyses_root / 'performance_benchmark.html')

    print(f"Scanning: {analyses_root}")
    records = collect_all(analyses_root)

    if not records:
        print("No analyses with performance_report.json found.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(records)} analyses:")
    for r in records:
        print(f"  - {r['analysis_name']}: {r.get('sequences_total', 0):,} sequences, "
              f"{r.get('total_wall_time_s', 0):.1f}s")

    sys_info = get_system_info()
    print(f"\nHost: {sys_info.get('cpu', '?')} | {sys_info.get('ram_gb', '?')} GB RAM | "
          f"{sys_info.get('cores_physical', '?')} cores")

    build_report(records, sys_info, output_path)
