import subprocess
import sys
import os
import signal
from pathlib import Path

from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot, QTimer

from venv_manager import get_venv_python

DOCKER_IMAGE = "loreon:latest"


def _app_root() -> Path:
    """Return the application root directory.

    When running from a PyInstaller bundle, __file__ points inside
    _internal/ but the .py data files live one level up next to the
    executable.  When running from source, __file__'s parent is correct.
    """
    if getattr(sys, 'frozen', False):
        # PyInstaller: executable is /opt/loreon/LOREON
        return Path(sys.executable).parent
    return Path(__file__).parent


class PipelineWorker(QObject):

    log_signal = pyqtSignal(str)
    progress_signal = pyqtSignal(int, str)
    finished_signal = pyqtSignal(bool, str)

    def __init__(self):
        super().__init__()
        self.is_running = False
        self.process = None
        self.settings = {}

    @pyqtSlot(dict)
    def run(self, settings):
        if self.is_running:
            return

        self.is_running = True
        self.settings = settings

        try:
            use_docker = settings.get("use_docker", False)

            cmd_pipeline = self._build_docker_pipeline_cmd() if use_docker else self._build_native_pipeline_cmd()
            success_pipeline = self.run_process(cmd_pipeline)

            if not self.is_running or not success_pipeline:
                self.finished_signal.emit(False, None)
                self.is_running = False
                return

            self.progress_signal.emit(95, "Step 7: HTML Report generation...")

            cmd_report = self._build_docker_report_cmd() if use_docker else self._build_native_report_cmd()
            success_report = self.run_process(cmd_report)

            if not success_report:
                self.log_signal.emit("[ERROR] report_generator.py failed.")
                self.finished_signal.emit(True, None)
                self.is_running = False
                return

            self.progress_signal.emit(100, "Completed!")
            self.finished_signal.emit(True, str(self._html_report_path()))

        except Exception as e:
            self.log_signal.emit("--- WORKER CRITICAL ERROR ---")
            self.log_signal.emit(str(e))
            self.finished_signal.emit(False, None)

        finally:
            self.is_running = False

    # --- Command builders ---

    def _pipeline_args(self, input_dir, output_dir, db_path):
        cmd = [
            "-i", str(input_dir),
            "-o", str(output_dir),
            "-d", str(db_path),
            "-f", self.settings["format"],
            "-T", str(self.settings["total_threads"]),
            "-t", str(self.settings["threads_per_job"]),
            "-k", str(self.settings["kmer_size"]),
            "-w", str(self.settings["window_size"]),
        ]
        if self.settings["debug"]:
            cmd.append("--debug")
        if self.settings["enable_filter"]:
            cmd.extend(["-min", str(self.settings["min_len"]),
                        "-max", str(self.settings["max_len"])])
        else:
            self.log_signal.emit("[INFO] Length filter disabled.")
            cmd.extend(["-min", "0", "-max", "999999"])
        cmd.extend(["--min-percent-identity", str(self.settings.get("min_percent_identity", 95.0))])
        cmd.extend(["--min-ref-coverage", str(self.settings.get("min_ref_coverage", 90.0))])
        if self.settings.get("force_tax_map", False):
            cmd.append("--force-tax-map")
        if self.settings.get("profile", False):
            cmd.append("--profile")
        if self.settings.get("delete_bam", False):
            cmd.append("--delete-bam")
        if self.settings.get("flat_sh", False):
            cmd.append("--flat-sh")
        return cmd

    def _report_names(self):
        output_dir = Path(self.settings["output_dir"])
        if self.settings.get("use_docker", False):
            analysis_name = "output"
        else:
            analysis_name = output_dir.name
        db_name = Path(self.settings["db_path"]).stem
        min_len = 0 if not self.settings["enable_filter"] else self.settings["min_len"]
        max_len = 999999 if not self.settings["enable_filter"] else self.settings["max_len"]
        return {
            # OPT-7: pipeline now writes TSV; report_generator auto-detects both
            "filter_report": f"report_filtering_{min_len}_{max_len}.tsv",
            "otu_table": f"OTU_Table_{analysis_name}_{db_name}.xlsx",
            "html_report": f"Report_{analysis_name}_{db_name}.html",
            "mapping_stats": "mapping_stats.json",
            "project_title": f"Report: {analysis_name} (DB: {db_name})",
        }

    def _html_report_path(self):
        names = self._report_names()
        return Path(self.settings["output_dir"]) / names["html_report"]

    def _get_python(self):
        """Return the Python interpreter to use for running pipeline scripts.

        When frozen (PyInstaller), uses the venv interpreter managed by
        venv_manager.  When running from source, uses sys.executable.
        """
        return get_venv_python()

    def _build_native_pipeline_cmd(self):
        script_path = str(_app_root() / "metaGenomics_new.py")
        # -u: force unbuffered stdout/stderr so log lines appear in real time
        return [self._get_python(), "-u", script_path] + self._pipeline_args(
            self.settings["input_dir"],
            self.settings["output_dir"],
            self.settings["db_path"],
        )

    def _build_native_report_cmd(self):
        script_path = str(_app_root() / "report_generator.py")
        output_dir = Path(self.settings["output_dir"])
        names = self._report_names()
        mapping_stats_path = output_dir / names["mapping_stats"]
        min_len = 0 if not self.settings["enable_filter"] else self.settings["min_len"]
        max_len = 999999 if not self.settings["enable_filter"] else self.settings["max_len"]
        cmd = [
            self._get_python(), "-u", script_path,
            "-f", str(output_dir / names["filter_report"]),
            "-otu", str(output_dir / names["otu_table"]),
            "-o", str(output_dir / names["html_report"]),
            "-pn", names["project_title"],
            "--min-len", str(min_len),
            "--max-len", str(max_len),
        ]
        if mapping_stats_path.exists():
            cmd.extend(["-ms", str(mapping_stats_path)])
        return cmd

    @staticmethod
    def _docker_path(p: Path) -> str:
        """Convert a host path to Docker-compatible format (forward slashes)."""
        return str(p).replace('\\', '/')

    def _docker_base_cmd(self):
        input_dir = Path(self.settings["input_dir"])
        output_dir = Path(self.settings["output_dir"])
        db_path = Path(self.settings["db_path"])
        return [
            "docker", "run", "--rm",
            "-v", f"{self._docker_path(input_dir)}:/data/input:ro",
            "-v", f"{self._docker_path(output_dir)}:/data/output",
            "-v", f"{self._docker_path(db_path.parent)}:/data/db",
            DOCKER_IMAGE,
        ]

    def _build_docker_pipeline_cmd(self):
        db_filename = Path(self.settings["db_path"]).name
        args = self._pipeline_args(
            input_dir="/data/input",
            output_dir="/data/output",
            db_path=f"/data/db/{db_filename}",
        )
        return self._docker_base_cmd() + ["python3", "/app/metaGenomics_new.py"] + args

    def _build_docker_report_cmd(self):
        names = self._report_names()
        output_dir = Path(self.settings["output_dir"])
        mapping_stats_path = output_dir / names["mapping_stats"]
        min_len = 0 if not self.settings["enable_filter"] else self.settings["min_len"]
        max_len = 999999 if not self.settings["enable_filter"] else self.settings["max_len"]
        cmd = self._docker_base_cmd() + [
            "python3", "/app/report_generator.py",
            "-f", f"/data/output/{names['filter_report']}",
            "-otu", f"/data/output/{names['otu_table']}",
            "-o", f"/data/output/{names['html_report']}",
            "-pn", names["project_title"],
            "--min-len", str(min_len),
            "--max-len", str(max_len),
        ]
        if mapping_stats_path.exists():
            cmd.extend(["-ms", f"/data/output/{names['mapping_stats']}"])
        return cmd

    # --- Process execution ---

    def run_process(self, command_list):
        if not self.is_running:
            return False
        if self.settings.get('debug', False):
            cmd_str = ' '.join(f'"{arg}"' if ' ' in arg else arg for arg in command_list)
            self.log_signal.emit(f"[DEBUG] Command: {cmd_str}")
        local_return_code = -1
        try:
            popen_kwargs = {}
            if sys.platform != 'win32':
                popen_kwargs['start_new_session'] = True

            self.process = subprocess.Popen(
                command_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, encoding='utf-8', errors='replace', bufsize=1,
                **popen_kwargs
            )
            with self.process.stdout:
                for line in iter(self.process.stdout.readline, ''):
                    line = line.strip()
                    if line:
                        self.log_signal.emit(line)
                        self.update_progress_from_log(line)
                    if not self.is_running:
                        break
            # OPT-15: no timeout — stdout loop already blocks until process
            # closes its pipe; a timeout here would kill long-running jobs
            self.process.wait()
            local_return_code = self.process.returncode
        except Exception as e:
            self.log_signal.emit(f"[SUBPROCESS Error] {e}")
            self.is_running = False
            return False
        finally:
            self.process = None
        if not self.is_running:
            self.log_signal.emit(f"Process ended (code: {local_return_code}). Request interrupt.")
            return False
        return local_return_code == 0

    # Step names updated to integer numbering
    def update_progress_from_log(self, line):
        if "STEP 1:" in line:
            self.progress_signal.emit(5, "Step 1: Parsing Taxonomy...")
        elif "STEP 2:" in line:
            self.progress_signal.emit(10, "Step 2: Filter...")
        elif "STEP 3:" in line:
            self.progress_signal.emit(40, "Step 3: Filter Report...")
        elif "STEP 4:" in line:
            self.progress_signal.emit(50, "Step 4: Mapping/Tabeling...")
        elif "STEP 5:" in line:
            self.progress_signal.emit(85, "Step 5: Reformatting...")
        elif "STEP 6:" in line:
            self.progress_signal.emit(90, "Step 6: OTU Aggregation...")

    @pyqtSlot()
    def stop(self):
        """Stop the pipeline. Invoked via Qt.DirectConnection so this runs on
        the GUI thread even while the worker thread is blocked in a readline
        loop inside run_process().

        Strategy: send SIGTERM to the whole process group (so minimap2,
        samtools and mp.Pool workers all receive it), then escalate to
        SIGKILL after 3 s if the tree is still alive — multiprocessing
        workers blocked on subprocess.communicate() do not always react to
        SIGTERM in time.
        """
        if not self.is_running:
            self.log_signal.emit("--- Request Interrupt (not running) ---")
            return

        self.is_running = False
        proc = self.process  # snapshot (worker thread may null it)
        if proc is None:
            self.log_signal.emit("--- Request Interrupt (no process active) ---")
            return

        self.log_signal.emit("--- STOPPING PIPELINE (killing process tree) ---")
        pid = proc.pid

        if sys.platform == 'win32':
            try:
                subprocess.run(
                    ['taskkill', '/F', '/T', '/PID', str(pid)],
                    capture_output=True
                )
            except Exception as e:
                self.log_signal.emit(f"Error during stop (taskkill): {e}")
            return

        # POSIX: send SIGTERM to the process group, then SIGKILL after a grace period
        try:
            pgid = os.getpgid(pid)
        except (ProcessLookupError, OSError):
            return

        try:
            os.killpg(pgid, signal.SIGTERM)
        except (ProcessLookupError, OSError):
            return
        except Exception as e:
            self.log_signal.emit(f"Error during stop (SIGTERM): {e}")
            try:
                proc.terminate()
            except Exception:
                pass
            return

        def _force_kill():
            # Escalation: if the tree is still alive after 3 s, send SIGKILL.
            try:
                os.killpg(pgid, 0)  # probe (0 = no-op, raises if dead)
            except (ProcessLookupError, OSError):
                return              # already dead, nothing to do
            try:
                os.killpg(pgid, signal.SIGKILL)
                self.log_signal.emit(
                    "--- SIGKILL sent (process tree did not exit on SIGTERM) ---"
                )
            except (ProcessLookupError, OSError):
                pass

        # Runs on the GUI thread (same thread that called stop via
        # DirectConnection), so QTimer has a live event loop.
        QTimer.singleShot(3000, _force_kill)
