#!/usr/bin/env python3
"""
LOREON vs Epi2me Performance Comparison
=========================================
Parses Epi2me trace.txt + params.json and LOREON performance_report.json
to produce a side-by-side benchmark HTML report.

Usage:
    python3 compare_epi2me_loreon.py \
        --epi2me /path/to/epi2me_sample \
        --loreon /path/to/ANALISI/Synthetic16S_Final \
        --output comparison_report.html
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
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio


# ── System info ─────────────────────────────────────────────────────────────

def get_system_info():
    """Collect host machine specs (works on macOS, Linux, Windows)."""
    info = {
        'os': f"{platform.system()} {platform.release()}",
        'arch': platform.machine(),
        'python': platform.python_version(),
    }
    try:
        if platform.system() == 'Darwin':
            info['cpu'] = subprocess.check_output(
                ['sysctl', '-n', 'machdep.cpu.brand_string'], text=True
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


# ── Duration parsing ────────────────────────────────────────────────────────

def parse_duration_to_seconds(duration_str):
    """Convert Nextflow duration strings like '58m 38s', '27.5s', '403ms' to seconds."""
    duration_str = duration_str.strip()
    total = 0.0

    # Hours
    m = re.search(r'(\d+(?:\.\d+)?)\s*h', duration_str)
    if m:
        total += float(m.group(1)) * 3600

    # Minutes
    m = re.search(r'(\d+(?:\.\d+)?)\s*m(?!s)', duration_str)
    if m:
        total += float(m.group(1)) * 60

    # Seconds
    m = re.search(r'(\d+(?:\.\d+)?)\s*s(?!$|[a-z])', duration_str)
    if m:
        total += float(m.group(1))
    # Also match standalone "Xs" at end
    m = re.search(r'(\d+(?:\.\d+?)?)\s*s$', duration_str)
    if m and total == 0.0:
        total = float(m.group(1))

    # Milliseconds only
    m = re.search(r'(\d+(?:\.\d+)?)\s*ms', duration_str)
    if m and total == 0.0:
        total = float(m.group(1)) / 1000.0

    return total


def parse_memory_to_mb(mem_str):
    """Convert Nextflow memory strings like '5.2 GB', '93.1 MB' to MB."""
    mem_str = mem_str.strip()
    m = re.match(r'([\d.]+)\s*(KB|MB|GB|TB)', mem_str, re.IGNORECASE)
    if not m:
        return 0.0
    val = float(m.group(1))
    unit = m.group(2).upper()
    if unit == 'KB':
        return val / 1024
    elif unit == 'MB':
        return val
    elif unit == 'GB':
        return val * 1024
    elif unit == 'TB':
        return val * 1024 * 1024
    return val


# ── Epi2me data collection ──────────────────────────────────────────────────

def collect_epi2me(epi2me_dir):
    """Parse Epi2me trace.txt and params.json."""
    epi2me_dir = Path(epi2me_dir)
    data = {'tool': 'Epi2me (wf-16s)'}

    # params.json
    params_path = epi2me_dir / 'params.json'
    if params_path.exists():
        with open(params_path, 'r', encoding='utf-8') as f:
            params = json.load(f)
        data['classifier'] = params.get('classifier', '?')
        data['database_set'] = params.get('database_set', '?')
        data['min_len'] = params.get('min_len', '?')
        data['max_len'] = params.get('max_len', '?')
        data['threads'] = params.get('threads', '?')
        data['min_percent_identity'] = params.get('min_percent_identity', '?')
        data['min_ref_coverage'] = params.get('min_ref_coverage', '?')
        data['agent'] = params.get('wf', {}).get('agent', '?')

    # versions_all.txt
    versions_path = epi2me_dir / 'versions_all.txt'
    if versions_path.exists():
        versions = {}
        with open(versions_path, 'r', encoding='utf-8') as f:
            for line in f:
                if ',' in line:
                    k, v = line.strip().split(',', 1)
                    versions[k.strip()] = v.strip()
        data['minimap2_version'] = versions.get('minimap2', '?')
        data['samtools_version'] = versions.get('samtools', '?')
        data['versions'] = versions

    # trace.txt (Nextflow execution trace)
    trace_path = epi2me_dir / 'execution' / 'trace.txt'
    tasks = []
    if trace_path.exists():
        df_trace = pd.read_csv(trace_path, sep='\t')
        for _, row in df_trace.iterrows():
            task = {
                'name': row.get('name', ''),
                'status': row.get('status', ''),
                'duration_s': parse_duration_to_seconds(str(row.get('duration', '0s'))),
                'realtime_s': parse_duration_to_seconds(str(row.get('realtime', '0s'))),
                'cpu_pct': float(str(row.get('%cpu', '0')).replace('%', '')),
                'peak_rss_mb': parse_memory_to_mb(str(row.get('peak_rss', '0 MB'))),
                'peak_vmem_mb': parse_memory_to_mb(str(row.get('peak_vmem', '0 MB'))),
            }
            tasks.append(task)
        data['tasks'] = tasks
        data['total_realtime_s'] = sum(t['realtime_s'] for t in tasks)
        data['peak_rss_mb'] = max((t['peak_rss_mb'] for t in tasks), default=0)

        # Wall time: from first submit to last completion
        if 'submit' in df_trace.columns:
            try:
                submits = pd.to_datetime(df_trace['submit'])
                durations = df_trace['duration'].apply(lambda d: parse_duration_to_seconds(str(d)))
                end_times = submits + pd.to_timedelta(durations, unit='s')
                data['wall_time_s'] = (end_times.max() - submits.min()).total_seconds()
            except Exception:
                data['wall_time_s'] = data['total_realtime_s']

    # abundance table (for number of taxa found)
    abund_path = epi2me_dir / 'abundance_table_genus.tsv'
    if abund_path.exists():
        df_abund = pd.read_csv(abund_path, sep='\t')
        data['n_taxa_genus'] = len(df_abund)
        if 'total' in df_abund.columns:
            data['total_assigned_reads'] = int(df_abund['total'].sum())

    return data


# ── LOREON data collection ──────────────────────────────────────────────────

def collect_loreon(loreon_dir):
    """Parse LOREON performance data."""
    loreon_dir = Path(loreon_dir)
    data = {'tool': 'LOREON'}

    # performance_report.json
    perf_path = loreon_dir / 'performance_report.json'
    if perf_path.exists():
        with open(perf_path, 'r', encoding='utf-8') as f:
            perf = json.load(f)
        data['run_date'] = perf.get('run_date', '')
        data['total_wall_time_s'] = perf.get('total_wall_time_s', 0)
        data['peak_memory_mb'] = perf.get('peak_memory_mb', 0)

        steps = []
        for step_name, step_data in perf.get('steps', {}).items():
            steps.append({
                'name': step_name,
                'wall_time_s': step_data.get('wall_time_s', 0),
                'mem_mb': step_data.get('mem_mb', 0),
            })
            # Extra fields
            if 'taxonomy_entries' in step_data:
                data['db_entries'] = step_data['taxonomy_entries']
            if 'sequences_total' in step_data:
                data['sequences_total'] = step_data['sequences_total']
            if 'sequences_good' in step_data:
                data['sequences_good'] = step_data['sequences_good']
            if 'files_to_map' in step_data:
                data['n_samples'] = step_data['files_to_map']
        data['steps'] = steps

    # mapping_stats.json
    ms_path = loreon_dir / 'mapping_stats.json'
    if ms_path.exists():
        with open(ms_path, 'r', encoding='utf-8') as f:
            ms = json.load(f)
        summary = ms.get('summary', {})
        data['mapped_reads'] = summary.get('total_mapped_all', 0)
        data['mapping_rate_pct'] = summary.get('overall_mapping_rate_pct', 0)

    # pipeline.log for extra info
    log_path = loreon_dir / 'pipeline.log'
    if log_path.exists():
        with open(log_path, 'r', encoding='utf-8') as f:
            for line in f:
                if 'Command:' in line and 'metaGenomics' in line:
                    cmd = line.split('Command:')[1].strip()
                    m = re.search(r'-d\s+(\S+)', cmd)
                    if m:
                        data['db_name'] = Path(m.group(1)).stem
                    m = re.search(r'-f\s+(\S+)', cmd)
                    if m:
                        data['format'] = m.group(1)
                    m = re.search(r'-T\s+(\d+)', cmd)
                    if m:
                        data['total_threads'] = int(m.group(1))
                    m = re.search(r'-t\s+(\d+)', cmd)
                    if m:
                        data['threads_per_job'] = int(m.group(1))
                    m = re.search(r'-k\s+(\d+)', cmd)
                    if m:
                        data['kmer_size'] = int(m.group(1))
                    m = re.search(r'-w\s+(\d+)', cmd)
                    if m:
                        data['window_size'] = int(m.group(1))
                    # minimap2 peak RSS from log
                for line2 in f:
                    pass
        # Re-read for minimap2 RSS
        with open(log_path, 'r', encoding='utf-8') as f:
            for line in f:
                m = re.search(r'Peak RSS:\s+([\d.]+)\s*(GB|MB)', line)
                if m:
                    val = float(m.group(1))
                    unit = m.group(2)
                    data['minimap2_peak_rss_mb'] = val * 1024 if unit == 'GB' else val

    # Filter report
    filter_reports = list(loreon_dir.glob('report_filtering_*.tsv'))
    if filter_reports:
        m = re.search(r'report_filtering_(\d+)_(\d+)', filter_reports[0].name)
        if m:
            data['min_len'] = int(m.group(1))
            data['max_len'] = int(m.group(2))
        try:
            df_fr = pd.read_csv(filter_reports[0], sep='\t')
            data['mean_seq_length_bp'] = round(df_fr['mean length'].mean(), 1)
        except Exception:
            pass

    return data


# ── Map Epi2me tasks to pipeline phases ─────────────────────────────────────

def map_epi2me_to_phases(epi_data):
    """Group Epi2me Nextflow tasks into comparable pipeline phases."""
    phases = {}
    tasks = epi_data.get('tasks', [])

    phase_mapping = {
        'Setup & DB Download': ['getParams', 'getVersions', 'getVersionsCommon',
                                'prepare_databases:download_reference_ref2taxid',
                                'prepare_databases:download_unpack_taxonomy'],
        'Read Ingestion (fastcat)': ['fastcat'],
        'Mapping (minimap2)': ['minimap_pipeline:minimap'],
        'Taxonomy Assignment': ['minimap_pipeline:minimapTaxonomy',
                                'minimap_pipeline:filter_references'],
        'Abundance Tables': ['minimap_pipeline:createAbundanceTables'],
        'Report Generation': ['makeReport'],
        'IGV Config': ['minimap_pipeline:configure_igv'],
    }

    for phase_name, task_prefixes in phase_mapping.items():
        matched = []
        for t in tasks:
            for prefix in task_prefixes:
                if t['name'].startswith(prefix) or prefix in t['name']:
                    matched.append(t)
                    break
        if matched:
            phases[phase_name] = {
                'realtime_s': sum(t['realtime_s'] for t in matched),
                'peak_rss_mb': max(t['peak_rss_mb'] for t in matched),
                'tasks': matched,
            }

    return phases


def map_loreon_to_phases(lor_data):
    """Map LOREON steps to comparable pipeline phases."""
    phases = {}
    steps = lor_data.get('steps', [])

    phase_mapping = {
        'Taxonomy Parsing': ['Step 1'],
        'Sequence Filter': ['Step 2', 'Step 3'],
        'Mapping (minimap2)': ['Step 4'],
        'Index Reformat': ['Step 5'],
        'OTU Aggregation': ['Step 6'],
    }

    for phase_name, step_prefixes in phase_mapping.items():
        matched = [s for s in steps if any(s['name'].startswith(p) for p in step_prefixes)]
        if matched:
            phases[phase_name] = {
                'realtime_s': sum(s['wall_time_s'] for s in matched),
                'peak_rss_mb': max(s['mem_mb'] for s in matched),
            }

    return phases


# ── HTML Report ─────────────────────────────────────────────────────────────

def build_comparison_report(epi_data, lor_data, sys_info, output_path):
    """Generate the HTML comparison report."""

    epi_phases = map_epi2me_to_phases(epi_data)
    lor_phases = map_loreon_to_phases(lor_data)

    charts_html = []

    # ── Color palette ──
    COLOR_LOREON = '#2196F3'
    COLOR_EPI2ME = '#FF9800'

    # ═══════════════════════════════════════════════════════════════════════
    # CHART 1: Total wall time comparison (horizontal bar)
    # ═══════════════════════════════════════════════════════════════════════
    epi_wall = epi_data.get('wall_time_s', sum(t['realtime_s'] for t in epi_data.get('tasks', [])))
    lor_wall = lor_data.get('total_wall_time_s', 0)

    fig_total = go.Figure()
    fig_total.add_trace(go.Bar(
        y=['LOREON', 'Epi2me'],
        x=[lor_wall, epi_wall],
        orientation='h',
        marker_color=[COLOR_LOREON, COLOR_EPI2ME],
        text=[f'{lor_wall:.1f}s ({lor_wall/60:.1f} min)', f'{epi_wall:.1f}s ({epi_wall/60:.1f} min)'],
        textposition='outside',
    ))
    fig_total.update_layout(
        title='Total Wall Time',
        xaxis_title='Time (seconds)',
        height=250,
        margin=dict(l=100, r=150),
        showlegend=False,
    )
    charts_html.append(pio.to_html(fig_total, full_html=False, include_plotlyjs='cdn'))

    # ═══════════════════════════════════════════════════════════════════════
    # CHART 2: Peak RSS Memory comparison
    # ═══════════════════════════════════════════════════════════════════════
    lor_peak = lor_data.get('minimap2_peak_rss_mb', lor_data.get('peak_memory_mb', 0))
    epi_peak = epi_data.get('peak_rss_mb', 0)

    fig_mem = go.Figure()
    fig_mem.add_trace(go.Bar(
        y=['LOREON', 'Epi2me'],
        x=[lor_peak, epi_peak],
        orientation='h',
        marker_color=[COLOR_LOREON, COLOR_EPI2ME],
        text=[f'{lor_peak:.0f} MB ({lor_peak/1024:.1f} GB)', f'{epi_peak:.0f} MB ({epi_peak/1024:.1f} GB)'],
        textposition='outside',
    ))
    fig_mem.update_layout(
        title='Peak Memory Usage (RSS)',
        xaxis_title='Memory (MB)',
        height=250,
        margin=dict(l=100, r=180),
        showlegend=False,
    )
    charts_html.append(pio.to_html(fig_mem, full_html=False, include_plotlyjs=False))

    # ═══════════════════════════════════════════════════════════════════════
    # CHART 3: Phase-by-phase time comparison (grouped bar)
    # ═══════════════════════════════════════════════════════════════════════

    # Comparable phases
    comparable = [
        ('Mapping (minimap2)', 'Mapping (minimap2)', 'Mapping (minimap2)'),
    ]

    # Build all phases for each tool
    all_phase_names = []
    lor_times = []
    epi_times = []

    # LOREON phases
    for phase_name, phase_data in lor_phases.items():
        if phase_name not in all_phase_names:
            all_phase_names.append(f'L: {phase_name}')
            lor_times.append(phase_data['realtime_s'])
            epi_times.append(0)

    # Epi2me phases
    for phase_name, phase_data in epi_phases.items():
        if phase_name not in all_phase_names:
            all_phase_names.append(f'E: {phase_name}')
            lor_times.append(0)
            epi_times.append(phase_data['realtime_s'])

    # Direct comparisons (mapping)
    mapping_comparison = ['Mapping (minimap2)']
    compare_names = []
    compare_lor = []
    compare_epi = []
    for phase in mapping_comparison:
        if phase in lor_phases and phase in epi_phases:
            compare_names.append(phase)
            compare_lor.append(lor_phases[phase]['realtime_s'])
            compare_epi.append(epi_phases[phase]['realtime_s'])

    fig_phases = go.Figure()
    fig_phases.add_trace(go.Bar(
        name='LOREON',
        x=compare_names,
        y=compare_lor,
        marker_color=COLOR_LOREON,
        text=[f'{v:.1f}s' for v in compare_lor],
        textposition='outside',
    ))
    fig_phases.add_trace(go.Bar(
        name='Epi2me',
        x=compare_names,
        y=compare_epi,
        marker_color=COLOR_EPI2ME,
        text=[f'{v:.1f}s' for v in compare_epi],
        textposition='outside',
    ))
    fig_phases.update_layout(
        barmode='group',
        title='Mapping Step: Direct Comparison (same minimap2, same data)',
        yaxis_title='Time (seconds)',
        height=400,
        legend=dict(orientation='h', yanchor='bottom', y=1.02),
    )
    charts_html.append(pio.to_html(fig_phases, full_html=False, include_plotlyjs=False))

    # ═══════════════════════════════════════════════════════════════════════
    # CHART 4: LOREON step breakdown (stacked bar)
    # ═══════════════════════════════════════════════════════════════════════
    lor_step_names = [s['name'] for s in lor_data.get('steps', [])]
    lor_step_times = [s['wall_time_s'] for s in lor_data.get('steps', [])]

    fig_lor_steps = go.Figure()
    fig_lor_steps.add_trace(go.Bar(
        x=lor_step_names,
        y=lor_step_times,
        marker_color=COLOR_LOREON,
        text=[f'{v:.1f}s' if v < 100 else f'{v:.0f}s' for v in lor_step_times],
        textposition='outside',
    ))
    fig_lor_steps.update_layout(
        title='LOREON: Step Breakdown',
        yaxis_title='Time (seconds)',
        yaxis_type='log',
        height=400,
        showlegend=False,
    )
    charts_html.append(pio.to_html(fig_lor_steps, full_html=False, include_plotlyjs=False))

    # ═══════════════════════════════════════════════════════════════════════
    # CHART 5: Epi2me task breakdown (stacked bar)
    # ═══════════════════════════════════════════════════════════════════════
    epi_tasks = epi_data.get('tasks', [])
    epi_task_names = [t['name'] for t in epi_tasks]
    epi_task_times = [t['realtime_s'] for t in epi_tasks]

    fig_epi_steps = go.Figure()
    fig_epi_steps.add_trace(go.Bar(
        x=epi_task_names,
        y=epi_task_times,
        marker_color=COLOR_EPI2ME,
        text=[f'{v:.1f}s' if v < 100 else f'{v:.0f}s' for v in epi_task_times],
        textposition='outside',
    ))
    fig_epi_steps.update_layout(
        title='Epi2me: Task Breakdown',
        yaxis_title='Time (seconds)',
        yaxis_type='log',
        height=400,
        showlegend=False,
        xaxis_tickangle=-30,
    )
    charts_html.append(pio.to_html(fig_epi_steps, full_html=False, include_plotlyjs=False))

    # ═══════════════════════════════════════════════════════════════════════
    # CHART 6: Memory per phase (grouped bar)
    # ═══════════════════════════════════════════════════════════════════════
    mem_phases = []
    mem_lor = []
    mem_epi = []

    for phase_name, phase_data in lor_phases.items():
        mem_phases.append(f'L: {phase_name}')
        mem_lor.append(phase_data['peak_rss_mb'])
        mem_epi.append(0)

    for phase_name, phase_data in epi_phases.items():
        mem_phases.append(f'E: {phase_name}')
        mem_lor.append(0)
        mem_epi.append(phase_data['peak_rss_mb'])

    fig_mem_phases = go.Figure()
    fig_mem_phases.add_trace(go.Bar(
        name='LOREON', x=mem_phases, y=mem_lor,
        marker_color=COLOR_LOREON,
    ))
    fig_mem_phases.add_trace(go.Bar(
        name='Epi2me', x=mem_phases, y=mem_epi,
        marker_color=COLOR_EPI2ME,
    ))
    fig_mem_phases.update_layout(
        barmode='group',
        title='Memory Usage per Pipeline Phase',
        yaxis_title='Peak RSS (MB)',
        height=450,
        xaxis_tickangle=-30,
        legend=dict(orientation='h', yanchor='bottom', y=1.02),
    )
    charts_html.append(pio.to_html(fig_mem_phases, full_html=False, include_plotlyjs=False))

    # ═══════════════════════════════════════════════════════════════════════
    # CHART 7: Resource efficiency radar
    # ═══════════════════════════════════════════════════════════════════════
    categories = ['Speed (1/time)', 'Low Memory', 'Low Disk Overhead',
                  'Ease of Setup', 'HPC Support', 'Resume Support']

    # Normalize: higher = better (0-10 scale)
    max_time = max(lor_wall, epi_wall)
    lor_scores = [
        10 * (1 - lor_wall / max_time) + 1,        # speed
        10 * (1 - lor_peak / max(lor_peak, epi_peak)),  # low memory
        9,                                           # no work/ directory
        9,                                           # GUI install
        2,                                           # no HPC
        2,                                           # no resume
    ]
    epi_scores = [
        10 * (1 - epi_wall / max_time) + 1,
        10 * (1 - epi_peak / max(lor_peak, epi_peak)),
        5,                                           # work/ directory
        5,                                           # Nextflow + Java + Docker
        9,                                           # Nextflow executors
        8,                                           # -resume
    ]

    fig_radar = go.Figure()
    fig_radar.add_trace(go.Scatterpolar(
        r=lor_scores + [lor_scores[0]],
        theta=categories + [categories[0]],
        fill='toself', name='LOREON',
        line_color=COLOR_LOREON, fillcolor='rgba(33,150,243,0.2)',
    ))
    fig_radar.add_trace(go.Scatterpolar(
        r=epi_scores + [epi_scores[0]],
        theta=categories + [categories[0]],
        fill='toself', name='Epi2me',
        line_color=COLOR_EPI2ME, fillcolor='rgba(255,152,0,0.2)',
    ))
    fig_radar.update_layout(
        polar=dict(radialaxis=dict(visible=True, range=[0, 10])),
        title='Resource Efficiency Comparison',
        height=500,
        legend=dict(orientation='h', yanchor='bottom', y=1.05),
    )
    charts_html.append(pio.to_html(fig_radar, full_html=False, include_plotlyjs=False))

    # ═══════════════════════════════════════════════════════════════════════
    # BUILD HTML
    # ═══════════════════════════════════════════════════════════════════════

    # Parameter comparison table
    param_rows = [
        ('Input Data', 'mbarc (Synthetic 16S)', 'mbarc (Synthetic 16S)'),
        ('Total Sequences', f'{lor_data.get("sequences_total", "?"):,}', f'{lor_data.get("sequences_total", "?"):,}'),
        ('Classifier', 'minimap2', epi_data.get('classifier', '?')),
        ('Database', lor_data.get('db_name', '?'), epi_data.get('database_set', '?')),
        ('DB Entries', f'{lor_data.get("db_entries", "?"):,}', 'N/A (NCBI pre-built)'),
        ('Min Length Filter', f'{lor_data.get("min_len", "?")} bp', f'{epi_data.get("min_len", "?")} bp'),
        ('Max Length Filter', f'{lor_data.get("max_len", "?")} bp', f'{epi_data.get("max_len", "?")} bp'),
        ('Min % Identity', f'{lor_data.get("min_percent_identity", 0)}', f'{epi_data.get("min_percent_identity", "?")}'),
        ('Min Ref Coverage', f'{lor_data.get("min_ref_coverage", 0)}', f'{epi_data.get("min_ref_coverage", "?")}'),
        ('Threads', f'{lor_data.get("total_threads", "?")} total, {lor_data.get("threads_per_job", "?")} per job', f'{epi_data.get("threads", "?")}'),
        ('minimap2 Version', '2.28-r1221', epi_data.get('minimap2_version', '?')),
        ('samtools Version', '1.21', epi_data.get('samtools_version', '?')),
        ('Mapped Reads', f'{lor_data.get("mapped_reads", "?"):,}', f'{epi_data.get("total_assigned_reads", "?"):,}'),
        ('Taxa Found (Genus)', '—', f'{epi_data.get("n_taxa_genus", "?")}'),
        ('Total Wall Time', f'{lor_wall:.1f}s ({lor_wall/60:.1f} min)', f'{epi_wall:.1f}s ({epi_wall/60:.1f} min)'),
        ('Peak Memory (RSS)', f'{lor_peak:.0f} MB', f'{epi_peak:.0f} MB'),
    ]

    param_table = '\n'.join(
        f'<tr><td><strong>{label}</strong></td><td>{loreon_val}</td><td>{epi_val}</td></tr>'
        for label, loreon_val, epi_val in param_rows
    )

    # Epi2me task detail table
    epi_task_rows = '\n'.join(
        f'<tr><td>{t["name"]}</td><td>{t["realtime_s"]:.1f}s</td>'
        f'<td>{t["cpu_pct"]:.1f}%</td><td>{t["peak_rss_mb"]:.0f} MB</td></tr>'
        for t in epi_tasks
    )

    # LOREON step detail table
    lor_step_rows = '\n'.join(
        f'<tr><td>{s["name"]}</td><td>{s["wall_time_s"]:.1f}s</td>'
        f'<td>—</td><td>{s["mem_mb"]:.0f} MB</td></tr>'
        for s in lor_data.get('steps', [])
    )

    # Speed difference
    if epi_wall > lor_wall:
        speed_diff = f'LOREON is <strong>{epi_wall/lor_wall:.1f}x faster</strong> in total wall time'
    else:
        speed_diff = f'Epi2me is <strong>{lor_wall/epi_wall:.1f}x faster</strong> in total wall time'

    if epi_peak > lor_peak:
        mem_diff = f'LOREON uses <strong>{epi_peak/lor_peak:.1f}x less peak memory</strong>'
    else:
        mem_diff = f'Epi2me uses <strong>{lor_peak/epi_peak:.1f}x less peak memory</strong>'

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>LOREON vs Epi2me — Performance Comparison</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {{ background-color: #f8f9fa; }}
        .container {{ max-width: 1400px; }}
        h1 {{ border-bottom: 2px solid #dee2e6; padding-bottom: 10px; margin-bottom: 20px; }}
        h2 {{ margin-top: 30px; color: #495057; }}
        .card {{ margin-bottom: 1rem; }}
        .table-wrapper {{ overflow-x: auto; }}
        .table th, .table td {{ text-align: left; }}
        .badge-loreon {{ background-color: {COLOR_LOREON}; color: #fff; font-size: 1rem; padding: 6px 14px; }}
        .badge-epi2me {{ background-color: {COLOR_EPI2ME}; color: #fff; font-size: 1rem; padding: 6px 14px; }}
        .kpi-card {{ text-align: center; padding: 20px; }}
        .kpi-card h3 {{ font-size: 2rem; margin-bottom: 5px; }}
        .kpi-card p {{ color: #6c757d; margin: 0; }}
    </style>
</head>
<body>
    <div class="container mt-4">
        <h1>LOREON vs Epi2me — Performance Comparison</h1>
        <p class="text-muted">Same input data (mbarc, Synthetic 16S mock community) &mdash; different pipelines</p>

        <!-- KPI Cards -->
        <div class="row mb-4">
            <div class="col-md-3">
                <div class="card shadow-sm kpi-card">
                    <span class="badge badge-loreon mb-2" style="display:inline-block">LOREON</span>
                    <h3>{lor_wall/60:.1f} min</h3>
                    <p>Total Wall Time</p>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card shadow-sm kpi-card">
                    <span class="badge badge-epi2me mb-2" style="display:inline-block">Epi2me</span>
                    <h3>{epi_wall/60:.1f} min</h3>
                    <p>Total Wall Time</p>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card shadow-sm kpi-card">
                    <span class="badge badge-loreon mb-2" style="display:inline-block">LOREON</span>
                    <h3>{lor_peak/1024:.1f} GB</h3>
                    <p>Peak Memory</p>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card shadow-sm kpi-card">
                    <span class="badge badge-epi2me mb-2" style="display:inline-block">Epi2me</span>
                    <h3>{epi_peak/1024:.1f} GB</h3>
                    <p>Peak Memory</p>
                </div>
            </div>
        </div>

        <!-- Key Findings -->
        <div class="alert alert-info">
            <strong>Key Findings:</strong> {speed_diff}. {mem_diff}.
            <br><em>Note: different databases were used (LOREON: SILVA 138.2 with {lor_data.get('db_entries', '?'):,} entries;
            Epi2me: NCBI 16S/18S pre-built). Both used minimap2 as classifier on the same input reads.</em>
        </div>

        <!-- Host Machine -->
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
                        </table>
                    </div>
                </div>
            </div>
            <div class="col-md-8">
                <div class="card shadow-sm h-100">
                    <div class="card-header"><strong>Parameter Comparison</strong></div>
                    <div class="card-body table-wrapper">
                        <table class="table table-sm table-striped mb-0">
                            <thead><tr><th>Parameter</th><th>LOREON</th><th>Epi2me</th></tr></thead>
                            <tbody>{param_table}</tbody>
                        </table>
                    </div>
                </div>
            </div>
        </div>

        <!-- Charts -->
        <h2>Performance Charts</h2>
        {"".join(f'<div class="card shadow-sm"><div class="card-body">{c}</div></div>' for c in charts_html)}

        <!-- Detailed task tables side by side -->
        <h2>Detailed Step/Task Breakdown</h2>
        <div class="row">
            <div class="col-md-6">
                <div class="card shadow-sm">
                    <div class="card-header"><strong>LOREON Steps</strong></div>
                    <div class="card-body table-wrapper">
                        <table class="table table-sm table-striped mb-0">
                            <thead><tr><th>Step</th><th>Time</th><th>CPU</th><th>Memory</th></tr></thead>
                            <tbody>
                                {lor_step_rows}
                                <tr class="table-active"><td><strong>TOTAL</strong></td>
                                <td><strong>{lor_wall:.1f}s</strong></td><td>—</td>
                                <td><strong>{lor_data.get('peak_memory_mb', 0):.0f} MB</strong></td></tr>
                            </tbody>
                        </table>
                    </div>
                </div>
            </div>
            <div class="col-md-6">
                <div class="card shadow-sm">
                    <div class="card-header"><strong>Epi2me Tasks (Nextflow trace)</strong></div>
                    <div class="card-body table-wrapper">
                        <table class="table table-sm table-striped mb-0">
                            <thead><tr><th>Task</th><th>Time</th><th>CPU</th><th>Memory</th></tr></thead>
                            <tbody>
                                {epi_task_rows}
                                <tr class="table-active"><td><strong>TOTAL</strong></td>
                                <td><strong>{epi_wall:.1f}s</strong></td><td>—</td>
                                <td><strong>{epi_peak:.0f} MB</strong></td></tr>
                            </tbody>
                        </table>
                    </div>
                </div>
            </div>
        </div>

    </div>
</body>
</html>"""

    output_path = Path(output_path)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)

    print(f"\nComparison report generated: {output_path.resolve()}")
    print(f"  LOREON wall time: {lor_wall:.1f}s ({lor_wall/60:.1f} min)")
    print(f"  Epi2me wall time: {epi_wall:.1f}s ({epi_wall/60:.1f} min)")
    print(f"  LOREON peak RSS: {lor_peak:.0f} MB")
    print(f"  Epi2me peak RSS: {epi_peak:.0f} MB")


# ── Main ────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='LOREON vs Epi2me Performance Comparison')
    parser.add_argument(
        '--epi2me', '-e',
        default='/Users/robm4/Documents/AIDocuments/dottorato/epi2me_sample',
        help='Epi2me output directory (contains params.json, execution/trace.txt)'
    )
    parser.add_argument(
        '--loreon', '-l',
        default='/Users/robm4/Documents/AIDocuments/dottorato/ANALISI/Synthetic16S_Final',
        help='LOREON analysis directory (contains performance_report.json)'
    )
    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Output HTML report path'
    )
    args = parser.parse_args()

    epi2me_dir = Path(args.epi2me)
    loreon_dir = Path(args.loreon)

    if not epi2me_dir.exists():
        print(f"ERROR: Epi2me directory not found: {epi2me_dir}", file=sys.stderr)
        sys.exit(1)
    if not loreon_dir.exists():
        print(f"ERROR: LOREON directory not found: {loreon_dir}", file=sys.stderr)
        sys.exit(1)

    output_path = args.output or str(
        Path(loreon_dir).parent / 'loreon_vs_epi2me_comparison.html'
    )

    print(f"Epi2me data: {epi2me_dir}")
    print(f"LOREON data: {loreon_dir}")

    epi_data = collect_epi2me(epi2me_dir)
    lor_data = collect_loreon(loreon_dir)
    sys_info = get_system_info()

    print(f"\nHost: {sys_info.get('cpu', '?')} | {sys_info.get('ram_gb', '?')} GB RAM")

    build_comparison_report(epi_data, lor_data, sys_info, output_path)
