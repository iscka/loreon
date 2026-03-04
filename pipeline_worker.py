import subprocess
import sys
from pathlib import Path

from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot

DOCKER_IMAGE = "loreon:latest"


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

            self.progress_signal.emit(95, "Step 4: HTML Report generation...")

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
        if self.settings.get("force_tax_map", False):
            cmd.append("--force-tax-map")
        return cmd

    def _report_names(self):
        output_dir = Path(self.settings["output_dir"])
        analysis_name = output_dir.name
        db_name = Path(self.settings["db_path"]).stem
        min_len = 0 if not self.settings["enable_filter"] else self.settings["min_len"]
        max_len = 999999 if not self.settings["enable_filter"] else self.settings["max_len"]
        return {
            "filter_report": f"report_filtering_{min_len}_{max_len}.xlsx",
            "otu_table": f"OTU_Table_{analysis_name}_{db_name}.xlsx",
            "html_report": f"Report_{analysis_name}_{db_name}.html",
            "project_title": f"Report: {analysis_name} (DB: {db_name})",
        }

    def _html_report_path(self):
        names = self._report_names()
        return Path(self.settings["output_dir"]) / names["html_report"]

    def _build_native_pipeline_cmd(self):
        script_path = str(Path(__file__).parent / "metaGenomics_new.py")
        return [sys.executable, script_path] + self._pipeline_args(
            self.settings["input_dir"],
            self.settings["output_dir"],
            self.settings["db_path"],
        )

    def _build_native_report_cmd(self):
        script_path = str(Path(__file__).parent / "report_generator.py")
        output_dir = Path(self.settings["output_dir"])
        names = self._report_names()
        return [
            sys.executable, script_path,
            "-f", str(output_dir / names["filter_report"]),
            "-otu", str(output_dir / names["otu_table"]),
            "-o", str(output_dir / names["html_report"]),
            "-pn", names["project_title"],
        ]

    def _docker_base_cmd(self):
        input_dir = Path(self.settings["input_dir"])
        output_dir = Path(self.settings["output_dir"])
        db_path = Path(self.settings["db_path"])
        return [
            "docker", "run", "--rm",
            "-v", f"{input_dir}:/data/input:ro",
            "-v", f"{output_dir}:/data/output",
            "-v", f"{db_path.parent}:/data/db:ro",
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
        return self._docker_base_cmd() + [
            "python3", "/app/report_generator.py",
            "-f", f"/data/output/{names['filter_report']}",
            "-otu", f"/data/output/{names['otu_table']}",
            "-o", f"/data/output/{names['html_report']}",
            "-pn", names["project_title"],
        ]

    # --- Process execution (unchanged) ---

    def run_process(self, command_list):
        if not self.is_running:
            return False
        if self.settings.get('debug', False):
            cmd_str = ' '.join(f'"{arg}"' if ' ' in arg else arg for arg in command_list)
            self.log_signal.emit(f"[DEBUG] Command: {cmd_str}")
        local_return_code = -1
        try:
            self.process = subprocess.Popen(
                command_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, encoding='utf-8', errors='replace', bufsize=1
            )
            with self.process.stdout:
                for line in iter(self.process.stdout.readline, ''):
                    line = line.strip()
                    if line:
                        self.log_signal.emit(line)
                        self.update_progress_from_log(line)
                    if not self.is_running:
                        break
            try:
                self.process.wait(timeout=600)
            except subprocess.TimeoutExpired:
                self.log_signal.emit("[WARNING] Process timed out after 600 seconds. Terminating.")
                self.process.terminate()
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

    def update_progress_from_log(self, line):
        if "--- Step 0:" in line:
            self.progress_signal.emit(5, "Step 0: Parsing Taxonomy...")
        elif "--- Step 1:" in line:
            self.progress_signal.emit(10, "Step 1: Filter...")
        elif "--- Step 1.5:" in line:
            self.progress_signal.emit(40, "Step 1.5: Filter Report...")
        elif "--- Step 2:" in line:
            self.progress_signal.emit(50, "Step 2: Mapping/Tabeling...")
        elif "--- Step 2.5:" in line:
            self.progress_signal.emit(85, "Step 2.5: Reformatting...")
        elif "--- Step 3:" in line:
            self.progress_signal.emit(90, "Step 3: OTU Aggregation...")

    @pyqtSlot()
    def stop(self):
        if self.is_running:
            self.is_running = False
            if self.process:
                self.log_signal.emit("--- SENDING SIGTERM SIGNAL ---")
                try:
                    self.process.terminate()
                except Exception as e:
                    self.log_signal.emit(f"Error during sigterm: {e}")
            else:
                self.log_signal.emit("--- Request Interrupt (no process active) ---")
        else:
            self.log_signal.emit("--- Request Interrupt (not running) ---")
