import subprocess
import sys
from pathlib import Path

from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot


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
            python_exe = sys.executable
            script_path = str(Path(__file__).parent / "metaGenomics_new.py")

            cmd_list = [
                python_exe,
                script_path,
                "-i", self.settings["input_dir"],
                "-o", self.settings["output_dir"],
                "-d", self.settings["db_path"],
                "-f", self.settings["format"],
                "-T", str(self.settings["total_threads"]),
                "-t", str(self.settings["threads_per_job"]),
                "-k", str(self.settings["kmer_size"]),
                "-w", str(self.settings["window_size"]),
            ]

            if self.settings["debug"]:
                cmd_list.append("--debug")

            if self.settings["enable_filter"]:
                cmd_list.extend(["-min", str(self.settings["min_len"])])
                cmd_list.extend(["-max", str(self.settings["max_len"])])
            else:
                self.log_signal.emit("[INFO] Length filter disabled.")
                cmd_list.extend(["-min", "0"])
                cmd_list.extend(["-max", "999999"])

            if self.settings.get("force_tax_map", False):
                cmd_list.append("--force-tax-map")

            success_pipeline = self.run_process(cmd_list)

            if not self.is_running or not success_pipeline:
                self.finished_signal.emit(False, None)
                self.is_running = False
                return

            self.progress_signal.emit(95, "Step 4: HTML Report generation...")

            report_script_path = str(Path(__file__).parent / "report_generator.py")

            analysis_name = Path(self.settings["output_dir"]).name
            db_name = Path(self.settings["db_path"]).stem

            min_len_report = 0 if not self.settings["enable_filter"] else self.settings["min_len"]
            max_len_report = 999999 if not self.settings["enable_filter"] else self.settings["max_len"]
            filter_report = Path(
                self.settings["output_dir"]) / f"report_filtering_{min_len_report}_{max_len_report}.xlsx"
            otu_filename = f"OTU_Table_{analysis_name}_{db_name}.xlsx"
            otu_table = Path(self.settings["output_dir"]) / otu_filename

            html_report_name = f"Report_{analysis_name}_{db_name}.html"
            html_report_path = Path(self.settings["output_dir"]) / html_report_name
            project_title = f"Report: {analysis_name} (DB: {db_name})"

            cmd_report = [
                python_exe,
                report_script_path,
                "-f", str(filter_report),
                "-otu", str(otu_table),
                "-o", str(html_report_path),
                "-pn", project_title
            ]

            success_report = self.run_process(cmd_report)

            if not success_report:
                self.log_signal.emit("[ERROR] report_generator.py failed.")
                self.finished_signal.emit(True, None)
                self.is_running = False
                return

            self.progress_signal.emit(100, "Completed!")
            self.finished_signal.emit(True, str(html_report_path))

        except Exception as e:
            self.log_signal.emit("--- WORKER CRITICAL ERROR ---")
            self.log_signal.emit(str(e))
            self.finished_signal.emit(False, None)

        finally:
            self.is_running = False

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
