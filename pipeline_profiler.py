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
        mem_mb = self._get_mem_mb()
        if mem_mb > self._peak_mem_mb:
            self._peak_mem_mb = mem_mb
        self._steps[name].update({
            'wall_time_s': round(wall, 3),
            'cpu_time_s': round(cpu, 3),
            'mem_mb': round(mem_mb, 1),
            **metrics
        })

    def _get_mem_mb(self) -> float:
        try:
            import resource
            return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        except ImportError:
            pass
        try:
            import psutil
            import os
            return psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)
        except Exception:
            pass
        try:
            with open('/proc/self/status') as f:
                for line in f:
                    if line.startswith('VmRSS:'):
                        return int(line.split()[1]) / 1024
        except Exception:
            pass
        return 0.0

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

        docx_path = self._save_word(output_dir, report)
        if docx_path:
            print(f"[PROFILER] Word report saved to: {docx_path}")

        return json_path

    def _save_word(self, output_dir: Path, report: dict) -> Path:
        try:
            from docx import Document
            from docx.shared import Pt, RGBColor, Cm
            from docx.enum.text import WD_ALIGN_PARAGRAPH
            from docx.enum.table import WD_ALIGN_VERTICAL
            from docx.oxml.ns import qn
            from docx.oxml import OxmlElement
        except ImportError:
            print("[PROFILER] WARNING: python-docx not installed. Word report skipped.")
            return None

        def set_cell_bg(cell, hex_color: str):
            tc = cell._tc
            tcPr = tc.get_or_add_tcPr()
            shd = OxmlElement('w:shd')
            shd.set(qn('w:fill'), hex_color)
            shd.set(qn('w:val'), 'clear')
            tcPr.append(shd)

        def set_cell_font(cell, bold=False, color=None, size=10):
            for para in cell.paragraphs:
                for run in para.runs:
                    run.bold = bold
                    run.font.size = Pt(size)
                    if color:
                        run.font.color.rgb = RGBColor(*color)

        doc = Document()

        # Page margins
        for section in doc.sections:
            section.top_margin = Cm(2)
            section.bottom_margin = Cm(2)
            section.left_margin = Cm(2.5)
            section.right_margin = Cm(2.5)

        # Title
        title = doc.add_heading('LOREON Pipeline — Performance Report', level=1)
        title.alignment = WD_ALIGN_PARAGRAPH.CENTER

        # Subtitle
        sub = doc.add_paragraph(f"Run date: {report['run_date']}")
        sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
        sub.runs[0].font.size = Pt(10)
        sub.runs[0].font.color.rgb = RGBColor(0x55, 0x55, 0x55)

        doc.add_paragraph()

        # Collect all extra metric keys across all steps
        extra_keys = []
        for name in self._order:
            for k in self._steps[name]:
                if k not in ('wall_time_s', 'cpu_time_s', 'mem_mb') and k not in extra_keys:
                    extra_keys.append(k)

        base_headers = ['Step', 'Wall Time (s)', 'CPU Time (s)', 'Mem (MB)']
        extra_headers = [k.replace('_', ' ').title() for k in extra_keys]
        all_headers = base_headers + extra_headers

        table = doc.add_table(rows=1, cols=len(all_headers))
        table.style = 'Table Grid'

        # Header row
        hdr_row = table.rows[0]
        for i, h in enumerate(all_headers):
            cell = hdr_row.cells[i]
            cell.text = h
            set_cell_bg(cell, '1F3864')
            cell.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER
            set_cell_font(cell, bold=True, color=(0xFF, 0xFF, 0xFF), size=10)

        # Data rows
        for idx, name in enumerate(self._order):
            s = self._steps[name]
            row = table.add_row()
            bg = 'DCE6F1' if idx % 2 == 0 else 'FFFFFF'

            values = [
                name,
                f"{s.get('wall_time_s', 0):.3f}",
                f"{s.get('cpu_time_s', 0):.3f}",
                f"{s.get('mem_mb', 0):.1f}",
            ] + [str(s.get(k, '—')) for k in extra_keys]

            for i, val in enumerate(values):
                cell = row.cells[i]
                cell.text = val
                set_cell_bg(cell, bg)
                align = WD_ALIGN_PARAGRAPH.LEFT if i == 0 else WD_ALIGN_PARAGRAPH.CENTER
                cell.paragraphs[0].alignment = align
                set_cell_font(cell, bold=(i == 0), size=9)

        # Total row
        total_row = table.add_row()
        total_values = [
            'TOTAL',
            f"{report['total_wall_time_s']:.3f}",
            '—',
            f"{report['peak_memory_mb']:.1f} (peak)",
        ] + ['—'] * len(extra_keys)

        for i, val in enumerate(total_values):
            cell = total_row.cells[i]
            cell.text = val
            set_cell_bg(cell, '2E4057')
            cell.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER if i > 0 else WD_ALIGN_PARAGRAPH.LEFT
            set_cell_font(cell, bold=True, color=(0xFF, 0xFF, 0xFF), size=10)

        doc.add_paragraph()
        note = doc.add_paragraph(
            'Wall Time: elapsed clock time per step.  '
            'CPU Time: actual CPU seconds consumed.  '
            'Mem: resident memory at step end.  '
            'Peak Mem: highest RSS across all steps.'
        )
        note.runs[0].font.size = Pt(8)
        note.runs[0].font.color.rgb = RGBColor(0x77, 0x77, 0x77)

        path = output_dir / 'profile_report.docx'
        doc.save(str(path))
        return path

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
