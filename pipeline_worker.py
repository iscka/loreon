import re
import subprocess
import sys
import os
import signal
import threading
import time
import uuid
from pathlib import Path

from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot, QTimer

from venv_manager import get_venv_python

DOCKER_IMAGE = "loreon:latest"
PARTIAL_FAILURE_EXIT_CODE = 3
_STEP4_PROGRESS = re.compile(r'Step 4 progress: (\d+)/(\d+) samples completed')
_SAMPLES_FAILED = re.compile(r'(\d+ of \d+ samples failed)')


def _app_root() -> Path:
    """Return the application root directory.

    When running from a PyInstaller bundle, __file__ points inside
    _internal/ but the .py data files live one level up next to the
    executable.  When running from source, __file__'s parent is correct.
    """
    if getattr(sys, 'frozen', False):
        return Path(sys.executable).parent
    return Path(__file__).parent


class PipelineWorker(QObject):

    log_signal = pyqtSignal(str)
    progress_signal = pyqtSignal(int, str)
    finished_signal = pyqtSignal(bool, str)
    warning_signal = pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.is_running = False
        self.process = None
        self.settings = {}
        self.last_return_code = None
        self._container_name = None
        self._kill_thread = None
        self._failed_samples_text = None

    @pyqtSlot(dict)
    def run(self, settings):
        if self.is_running:
            return

        self.is_running = True
        self.settings = settings
        self._failed_samples_text = None

        try:
            use_docker = settings.get("use_docker", False)

            cmd_pipeline = self._build_docker_pipeline_cmd() if use_docker else self._build_native_pipeline_cmd()
            success_pipeline = self.run_process(cmd_pipeline)
            partial_failure = (not success_pipeline
                               and self.last_return_code == PARTIAL_FAILURE_EXIT_CODE)

            if not self.is_running or not (success_pipeline or partial_failure):
                self.finished_signal.emit(False, None)
                self.is_running = False
                return

            if partial_failure:
                self.warning_signal.emit(
                    self._failed_samples_text or "Some samples failed"
                )

            self.progress_signal.emit(95, "Step 7: HTML Report generation...")

            cmd_report = self._build_docker_report_cmd() if use_docker else self._build_native_report_cmd()
            success_report = self.run_process(cmd_report)

            if not self.is_running:
                self.finished_signal.emit(False, None)
                return

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
        cmd.extend(["--min-percent-identity", str(self.settings.get("min_percent_identity", 0.0))])
        cmd.extend(["--min-ref-coverage", str(self.settings.get("min_ref_coverage", 0.0))])
        sample_timeout = self.settings.get("sample_timeout_min", 0)
        if sample_timeout and sample_timeout > 0:
            cmd.extend(["--sample-timeout", str(sample_timeout)])
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
        self._container_name = f"loreon-{uuid.uuid4().hex[:12]}"
        return [
            "docker", "run", "--rm",
            "--name", self._container_name,
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


    def run_process(self, command_list):
        if not self.is_running:
            return False
        if self.settings.get('debug', False):
            cmd_str = ' '.join(f'"{arg}"' if ' ' in arg else arg for arg in command_list)
            self.log_signal.emit(f"[DEBUG] Command: {cmd_str}")
        local_return_code = -1
        self.last_return_code = None
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
            if not self.is_running and self.process.poll() is None:
                self._request_container_stop()
                self._terminate_local_process(self.process)
            self.process.wait()
            local_return_code = self.process.returncode
            self.last_return_code = local_return_code
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

    def update_progress_from_log(self, line):
        if "STEP 1:" in line:
            self.progress_signal.emit(5, "Step 1: Parsing Taxonomy...")
        elif "STEP 2:" in line:
            self.progress_signal.emit(10, "Step 2: Filter...")
        elif "STEP 3:" in line:
            self.progress_signal.emit(40, "Step 3: Filter Report...")
        elif "STEP 4:" in line:
            self.progress_signal.emit(50, "Step 4: Mapping/Tabeling...")
        elif "Step 4 progress:" in line:
            match = _STEP4_PROGRESS.search(line)
            if match and int(match.group(2)) > 0:
                done, total = int(match.group(1)), int(match.group(2))
                self.progress_signal.emit(
                    50 + int(35 * done / total),
                    f"Step 4: Mapping/Tabeling... {done}/{total} samples"
                )
        elif "samples failed" in line:
            match = _SAMPLES_FAILED.search(line)
            if match:
                self._failed_samples_text = match.group(1)
        elif "STEP 5:" in line:
            self.progress_signal.emit(85, "Step 5: Reformatting...")
        elif "STEP 6:" in line:
            self.progress_signal.emit(90, "Step 6: OTU Aggregation...")

    def _kill_container(self, name):
        """Stop the Docker container that runs the pipeline.

        Ending the local 'docker run' client does not stop the container, so
        it is killed by name. The attempt is repeated for a few seconds to
        cover a container that was still being created when Stop was pressed.
        """
        extra = {}
        if sys.platform == 'win32':
            extra['creationflags'] = getattr(subprocess, 'CREATE_NO_WINDOW', 0)
        for _ in range(5):
            try:
                result = subprocess.run(
                    ['docker', 'kill', name],
                    capture_output=True, text=True, timeout=30, **extra
                )
            except (OSError, subprocess.TimeoutExpired) as e:
                self.log_signal.emit(f"Error during stop (docker kill): {e}")
                return
            if result.returncode == 0:
                self.log_signal.emit(f"--- Docker container {name} stopped ---")
                return
            time.sleep(1.0)

    def _request_container_stop(self):
        """Start stopping the pipeline container in the background, unless already requested."""
        if not (self.settings.get("use_docker", False) and self._container_name):
            return
        if self._kill_thread is not None and self._kill_thread.is_alive():
            return
        self._kill_thread = threading.Thread(
            target=self._kill_container, args=(self._container_name,), daemon=True
        )
        self._kill_thread.start()

    def wait_for_container_stop(self, timeout=None):
        """Wait until a requested container stop has completed. Returns False on timeout."""
        thread = self._kill_thread
        if thread is None:
            return True
        thread.join(timeout)
        return not thread.is_alive()

    @staticmethod
    def _terminate_local_process(proc):
        """End a pipeline process that is still running after a stop request.

        Used by the worker thread when the process was started while the stop
        request was being handled, so stop() could not reach it.
        """
        if sys.platform == 'win32':
            try:
                subprocess.run(
                    ['taskkill', '/F', '/T', '/PID', str(proc.pid)],
                    capture_output=True,
                    creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0)
                )
            except OSError:
                pass
            return
        try:
            pgid = os.getpgid(proc.pid)
            os.killpg(pgid, signal.SIGTERM)
        except (ProcessLookupError, OSError):
            return
        try:
            proc.wait(timeout=3)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(pgid, signal.SIGKILL)
            except (ProcessLookupError, OSError):
                pass

    @pyqtSlot()
    def stop(self):
        """Stop the pipeline. Invoked via Qt.DirectConnection so this runs on
        the GUI thread even while the worker thread is blocked in a readline
        loop inside run_process().

        Strategy: in Docker mode the container is killed by name from a
        background thread, so the interface is not blocked. Then SIGTERM is
        sent to the whole local process group (so minimap2, samtools and
        mp.Pool workers all receive it), escalating to SIGKILL after 3 s if
        the tree is still alive.
        """
        if not self.is_running:
            self.log_signal.emit("--- Request Interrupt (not running) ---")
            return

        self.is_running = False
        proc = self.process
        if proc is None:
            self.log_signal.emit("--- Request Interrupt (no process active) ---")
            return

        self.log_signal.emit("--- STOPPING PIPELINE (killing process tree) ---")
        pid = proc.pid

        self._request_container_stop()

        if sys.platform == 'win32':
            try:
                subprocess.run(
                    ['taskkill', '/F', '/T', '/PID', str(pid)],
                    capture_output=True,
                    creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0)
                )
            except Exception as e:
                self.log_signal.emit(f"Error during stop (taskkill): {e}")
            return

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
            try:
                os.killpg(pgid, 0)
            except (ProcessLookupError, OSError):
                return
            try:
                os.killpg(pgid, signal.SIGKILL)
                self.log_signal.emit(
                    "--- SIGKILL sent (process tree did not exit on SIGTERM) ---"
                )
            except (ProcessLookupError, OSError):
                pass

        QTimer.singleShot(3000, _force_kill)
