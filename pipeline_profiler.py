#!/usr/bin/env python3

import time
import json
from datetime import datetime
from pathlib import Path


class PipelineProfiler:

    def __init__(self):
        self._steps = {}
        self._order = []
        self._start_wall = None
        self._start_cpu = None
        self._peak_mem_mb = 0.0
        self._started_at = datetime.now()
        # OPT-12: detect memory backend once at init
        self._mem_fn = self._detect_mem_backend()

    @staticmethod
    def _detect_mem_backend():
        """Detect the best available memory measurement method once."""
        try:
            import resource
            import sys as _sys
            platform = _sys.platform

            def _fn():
                ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                return ru / (1024 * 1024) if platform == 'darwin' else ru / 1024

            _fn()  # test call
            return _fn
        except ImportError:
            pass
        try:
            import psutil
            import os
            proc = psutil.Process(os.getpid())

            def _fn():
                return proc.memory_info().rss / (1024 * 1024)

            _fn()  # test call
            return _fn
        except Exception:
            pass
        return lambda: 0.0

    def start_step(self, name: str):
        self._current = name
        self._start_wall = time.perf_counter()
        self._start_cpu = time.process_time()
        if name not in self._steps:
            self._order.append(name)
            self._steps[name] = {}

    def end_step(self, name: str, **metrics):
        wall = time.perf_counter() - self._start_wall
        cpu = time.process_time() - self._start_cpu
        mem_mb = self._mem_fn()
        if mem_mb > self._peak_mem_mb:
            self._peak_mem_mb = mem_mb
        self._steps[name].update({
            'wall_time_s': round(wall, 3),
            'cpu_time_s': round(cpu, 3),
            'mem_mb': round(mem_mb, 1),
            **metrics
        })

    def save(self, output_dir: Path) -> Path:
        output_dir = Path(output_dir)
        total_wall = round(sum(
            self._steps[n].get('wall_time_s', 0) for n in self._order
        ), 3)
        report = {
            'run_date': self._started_at.strftime('%Y-%m-%d %H:%M:%S'),
            'steps': {name: self._steps[name] for name in self._order},
            'total_wall_time_s': total_wall,
            'peak_memory_mb': round(self._peak_mem_mb, 1),
        }

        json_path = output_dir / 'profile_report.json'
        with open(json_path, 'w') as f:
            json.dump(report, f, indent=2)

        self._print_summary(report)
        print(f"[PROFILER] JSON report saved to: {json_path}")

        html_path = self._save_html(output_dir, report)
        if html_path:
            print(f"[PROFILER] HTML report saved to: {html_path}")

        return json_path

    # ------------------------------------------------------------------
    # HTML report (Bootstrap 5 + Plotly, consistent with main report)
    # ------------------------------------------------------------------

    def _save_html(self, output_dir: Path, report: dict) -> Path:
        run_date = report['run_date']
        total_wall = report['total_wall_time_s']
        peak_mem = report['peak_memory_mb']
        num_steps = len(self._order)

        # Format total wall time as human-readable
        def _fmt_time(seconds: float) -> str:
            if seconds < 60:
                return f"{seconds:.1f} s"
            m, s = divmod(seconds, 60)
            if m < 60:
                return f"{int(m)}m {s:.0f}s"
            h, m = divmod(int(m), 60)
            return f"{int(h)}h {int(m)}m {s:.0f}s"

        total_wall_fmt = _fmt_time(total_wall)

        # Collect all extra metric keys across all steps
        base_keys = {'wall_time_s', 'cpu_time_s', 'mem_mb'}
        extra_keys = []
        for name in self._order:
            for k in self._steps[name]:
                if k not in base_keys and k not in extra_keys:
                    extra_keys.append(k)

        # Build table rows HTML
        table_rows = []
        for idx, name in enumerate(self._order):
            s = self._steps[name]
            row_class = 'table-light' if idx % 2 == 0 else ''
            extra_cells = ''.join(
                f'<td class="text-center">{s.get(k, "—")}</td>'
                for k in extra_keys
            )
            table_rows.append(
                f'<tr class="{row_class}">'
                f'<td class="fw-bold">{name}</td>'
                f'<td class="text-center">{s.get("wall_time_s", 0):.3f}</td>'
                f'<td class="text-center">{s.get("cpu_time_s", 0):.3f}</td>'
                f'<td class="text-center">{s.get("mem_mb", 0):.1f}</td>'
                f'{extra_cells}'
                f'</tr>'
            )

        extra_headers = ''.join(
            f'<th class="text-center">{k.replace("_", " ").title()}</th>'
            for k in extra_keys
        )
        table_rows_html = '\n'.join(table_rows)

        extra_total_cells = ''.join(
            '<td class="text-center">—</td>' for _ in extra_keys
        )

        # Build progress bar segments (EPI2ME style)
        progress_segments = []
        colors = ['#198754', '#0d6efd', '#6f42c1', '#fd7e14',
                  '#20c997', '#0dcaf0', '#d63384', '#ffc107']
        for i, name in enumerate(self._order):
            s = self._steps[name]
            pct = (s.get('wall_time_s', 0) / total_wall * 100) if total_wall > 0 else 0
            color = colors[i % len(colors)]
            label = name.replace('Step ', '').split(':')[0] if ':' in name else name
            progress_segments.append(
                f'<div class="progress-bar" role="progressbar" '
                f'style="width:{pct:.1f}%; background-color:{color};" '
                f'data-bs-toggle="tooltip" title="{name}: {_fmt_time(s.get("wall_time_s",0))}">'
                f'<span class="text-truncate px-1" style="font-size:0.75rem;">'
                f'{label}</span></div>'
            )
        progress_bar_html = '\n'.join(progress_segments)

        # JSON data for Plotly charts
        step_names_json = json.dumps([n for n in self._order])
        wall_times_json = json.dumps([
            self._steps[n].get('wall_time_s', 0) for n in self._order
        ])
        cpu_times_json = json.dumps([
            self._steps[n].get('cpu_time_s', 0) for n in self._order
        ])
        mem_json = json.dumps([
            self._steps[n].get('mem_mb', 0) for n in self._order
        ])
        colors_json = json.dumps(
            [colors[i % len(colors)] for i in range(num_steps)]
        )

        # Waterfall/Gantt data: cumulative start offsets
        cumulative = 0.0
        gantt_starts = []
        gantt_durations = []
        for name in self._order:
            gantt_starts.append(round(cumulative, 3))
            dur = self._steps[name].get('wall_time_s', 0)
            gantt_durations.append(round(dur, 3))
            cumulative += dur
        gantt_starts_json = json.dumps(gantt_starts)
        gantt_durations_json = json.dumps(gantt_durations)

        # Throughput data (if 'throughput_seq_s' exists in any step)
        has_throughput = any(
            'throughput_seq_s' in self._steps[n] for n in self._order
        )
        throughput_names_json = json.dumps([
            n for n in self._order if 'throughput_seq_s' in self._steps[n]
        ])
        throughput_vals_json = json.dumps([
            self._steps[n]['throughput_seq_s']
            for n in self._order if 'throughput_seq_s' in self._steps[n]
        ])

        # Build the complete HTML document
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>LOREON Pipeline — Performance Report</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{ background-color: #f8f9fa; }}
        .container {{ max-width: 1400px; }}
        .card {{ margin-bottom: 1.5rem; }}
        .kpi-card {{ text-align: center; }}
        .kpi-value {{ font-size: 2.5rem; font-weight: bold; }}
        .kpi-label {{ font-size: 1rem; color: #6c757d; }}
        h1, h2 {{ border-bottom: 2px solid #dee2e6; padding-bottom: 10px; margin-top: 1.5rem; }}
        .table-wrapper {{ overflow-x: auto; white-space: nowrap; }}
    </style>
</head>
<body>
    <div class="container mt-4">

        <!-- Header -->
        <div class="card shadow-sm">
            <div class="card-body">
                <h1 class="display-5">LOREON Pipeline &mdash; Performance Report</h1>
                <p class="lead mb-2">Run date: {run_date}</p>
                <!-- EPI2ME-style progress bar -->
                <div class="progress" style="height: 1.6rem; border-radius: 0.25rem;">
                    {progress_bar_html}
                </div>
                <small class="text-muted">Total duration: <strong>{total_wall_fmt}</strong></small>
            </div>
        </div>

        <!-- KPI Cards -->
        <div class="row">
            <div class="col-md col-sm-6">
                <div class="card kpi-card shadow-sm"><div class="card-body">
                    <div class="kpi-value">{total_wall_fmt}</div>
                    <div class="kpi-label">Total Wall Time</div>
                </div></div>
            </div>
            <div class="col-md col-sm-6">
                <div class="card kpi-card shadow-sm"><div class="card-body">
                    <div class="kpi-value">{peak_mem:.1f} MB</div>
                    <div class="kpi-label">Peak Memory (RSS)</div>
                </div></div>
            </div>
            <div class="col-md col-sm-6">
                <div class="card kpi-card shadow-sm"><div class="card-body">
                    <div class="kpi-value">{num_steps}</div>
                    <div class="kpi-label">Pipeline Steps</div>
                </div></div>
            </div>
        </div>

        <!-- Section 1: Resource Usage -->
        <h2>Resource Usage</h2>

        <!-- Timeline / Gantt chart -->
        <div class="card shadow-sm">
            <div class="card-header"><h5 class="mb-0">Pipeline Timeline</h5></div>
            <div class="card-body">
                <div id="gantt-chart"></div>
            </div>
        </div>

        <!-- CPU & Wall Time (grouped bar) -->
        <div class="row">
            <div class="col-lg-8">
                <div class="card shadow-sm">
                    <div class="card-header"><h5 class="mb-0">Job Duration</h5></div>
                    <div class="card-body"><div id="time-chart"></div></div>
                </div>
            </div>
            <div class="col-lg-4">
                <div class="card shadow-sm">
                    <div class="card-header"><h5 class="mb-0">Time Distribution</h5></div>
                    <div class="card-body"><div id="pie-chart"></div></div>
                </div>
            </div>
        </div>

        <!-- Memory chart -->
        <div class="card shadow-sm">
            <div class="card-header"><h5 class="mb-0">Memory Usage per Step</h5></div>
            <div class="card-body"><div id="mem-chart"></div></div>
        </div>

        <!-- Throughput chart (only if data exists) -->
        <div id="throughput-section" style="display:none;">
            <div class="card shadow-sm">
                <div class="card-header"><h5 class="mb-0">Throughput</h5></div>
                <div class="card-body"><div id="throughput-chart"></div></div>
            </div>
        </div>

        <!-- Section 2: Step Details Table -->
        <h2>Step Details</h2>
        <div class="card shadow-sm">
            <div class="card-body table-wrapper">
                <table class="table table-striped table-hover table-sm">
                    <thead class="table-dark">
                        <tr>
                            <th>Step</th>
                            <th class="text-center">Wall Time (s)</th>
                            <th class="text-center">CPU Time (s)</th>
                            <th class="text-center">Mem (MB)</th>
                            {extra_headers}
                        </tr>
                    </thead>
                    <tbody>
                        {table_rows_html}
                    </tbody>
                    <tfoot class="table-dark">
                        <tr>
                            <td class="fw-bold">TOTAL</td>
                            <td class="text-center fw-bold">{total_wall:.3f}</td>
                            <td class="text-center">&mdash;</td>
                            <td class="text-center fw-bold">{peak_mem:.1f} (peak)</td>
                            {extra_total_cells}
                        </tr>
                    </tfoot>
                </table>
            </div>
        </div>

        <div class="alert alert-secondary mt-2" role="alert">
            <small>
                <strong>Wall Time:</strong> elapsed clock time per step. &nbsp;
                <strong>CPU Time:</strong> actual CPU seconds consumed by the main process. &nbsp;
                <strong>Mem:</strong> resident set size at step end. &nbsp;
                <strong>Peak Mem:</strong> highest RSS across all steps.
            </small>
        </div>

    </div>

    <script>
        const stepNames = {step_names_json};
        const wallTimes = {wall_times_json};
        const cpuTimes  = {cpu_times_json};
        const memMB     = {mem_json};
        const stepColors = {colors_json};
        const ganttStarts    = {gantt_starts_json};
        const ganttDurations = {gantt_durations_json};

        // --- Chart 1: Pipeline Timeline (horizontal Gantt) ---
        const ganttTrace = {{
            y: stepNames.slice().reverse(),
            x: ganttDurations.slice().reverse(),
            base: ganttStarts.slice().reverse(),
            type: 'bar',
            orientation: 'h',
            marker: {{ color: stepColors.slice().reverse() }},
            text: ganttDurations.slice().reverse().map(d => d.toFixed(1) + 's'),
            textposition: 'inside',
            insidetextanchor: 'middle',
            hovertemplate: '%{{y}}<br>Start: %{{base:.1f}}s<br>Duration: %{{x:.1f}}s<extra></extra>'
        }};
        Plotly.newPlot('gantt-chart', [ganttTrace], {{
            title: 'Step Execution Timeline',
            xaxis: {{ title: 'Elapsed Time (s)', rangemode: 'tozero' }},
            height: 60 + stepNames.length * 45,
            margin: {{ l: 200 }},
            showlegend: false
        }}, {{ responsive: true }});

        // --- Chart 2: Job Duration (grouped bar, EPI2ME style) ---
        Plotly.newPlot('time-chart', [
            {{ x: stepNames, y: wallTimes, name: 'Wall Time (s)',
               type: 'bar', marker: {{ color: '#1f77b4' }} }},
            {{ x: stepNames, y: cpuTimes, name: 'CPU Time (s)',
               type: 'bar', marker: {{ color: '#aec7e8' }} }}
        ], {{
            title: 'Wall Time vs CPU Time per Step',
            barmode: 'group', height: 420,
            yaxis: {{ title: 'Seconds', rangemode: 'tozero' }},
            xaxis: {{ tickangle: -25 }}
        }}, {{ responsive: true }});

        // --- Chart 3: Pie/Donut time distribution ---
        Plotly.newPlot('pie-chart', [{{
            labels: stepNames, values: wallTimes,
            type: 'pie', hole: 0.45,
            marker: {{ colors: stepColors }},
            textinfo: 'percent',
            hovertemplate: '%{{label}}<br>%{{value:.1f}}s (%{{percent}})<extra></extra>'
        }}], {{
            title: 'Time Distribution',
            height: 420, showlegend: true,
            legend: {{ orientation: 'h', y: -0.15 }}
        }}, {{ responsive: true }});

        // --- Chart 4: Memory usage (bar + peak line) ---
        const peakMem = {peak_mem};
        Plotly.newPlot('mem-chart', [
            {{ x: stepNames, y: memMB, type: 'bar', name: 'RSS at step end',
               marker: {{ color: '#2ca02c' }} }},
            {{ x: stepNames,
               y: Array(stepNames.length).fill(peakMem),
               type: 'scatter', mode: 'lines',
               name: 'Peak RSS (' + peakMem.toFixed(1) + ' MB)',
               line: {{ color: '#d62728', dash: 'dash', width: 2 }} }}
        ], {{
            title: 'Memory Usage',
            yaxis: {{ title: 'MB', rangemode: 'tozero' }},
            xaxis: {{ tickangle: -25 }},
            height: 380
        }}, {{ responsive: true }});

        // --- Chart 5: Throughput (if data exists) ---
        const throughputNames = {throughput_names_json};
        const throughputVals  = {throughput_vals_json};
        if (throughputNames.length > 0) {{
            document.getElementById('throughput-section').style.display = 'block';
            Plotly.newPlot('throughput-chart', [{{
                x: throughputNames, y: throughputVals,
                type: 'bar',
                marker: {{ color: '#ff7f0e' }},
                text: throughputVals.map(v => v.toLocaleString() + ' seq/s'),
                textposition: 'outside'
            }}], {{
                title: 'Processing Throughput',
                yaxis: {{ title: 'Sequences / second', rangemode: 'tozero' }},
                xaxis: {{ tickangle: -25 }},
                height: 380
            }}, {{ responsive: true }});
        }}
    </script>
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js"></script>
    <script>
        // Enable tooltips (Bootstrap 5) — must run AFTER bootstrap.bundle.min.js is loaded
        document.querySelectorAll('[data-bs-toggle="tooltip"]')
            .forEach(el => new bootstrap.Tooltip(el));
    </script>

</body>
</html>"""

        path = output_dir / 'profile_report.html'
        try:
            with open(path, 'w', encoding='utf-8') as f:
                f.write(html)
            return path
        except Exception as e:
            print(f"[PROFILER] ERROR saving HTML report: {e}")
            return None

    def _print_summary(self, report: dict):
        col = 22
        print("\n=== PERFORMANCE PROFILE ===")
        print(f"{'Step':<{col}} {'Wall (s)':>10} {'CPU (s)':>10} {'Mem (MB)':>10}  Metrics")
        print("-" * 80)
        for name in self._order:
            s = self._steps[name]
            extras = {
                k: v for k, v in s.items()
                if k not in ('wall_time_s', 'cpu_time_s', 'mem_mb')
            }
            extras_str = '  ' + ', '.join(f"{k}={v}" for k, v in extras.items()) if extras else ''
            print(
                f"{name:<{col}} "
                f"{s.get('wall_time_s', 0):>10.3f} "
                f"{s.get('cpu_time_s', 0):>10.3f} "
                f"{s.get('mem_mb', 0):>10.1f}"
                f"{extras_str}"
            )
        print("-" * 80)
        print(
            f"{'TOTAL':<{col}} "
            f"{report['total_wall_time_s']:>10.3f} "
            f"{'':>10} "
            f"{report['peak_memory_mb']:>10.1f}"
        )
        print("===========================\n")
