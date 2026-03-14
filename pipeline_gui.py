#!/usr/bin/env python3

import shutil
import subprocess
import sys
import multiprocessing
from datetime import datetime
from pathlib import Path
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QGroupBox, QLabel, QLineEdit, QPushButton,
    QComboBox, QCheckBox, QPlainTextEdit, QProgressBar, QMessageBox,
    QFileDialog, QSpinBox, QDoubleSpinBox, QDialog, QDialogButtonBox,
    QSystemTrayIcon, QMenu, QSplashScreen, QSizePolicy, QTextBrowser
)
from PyQt5.QtCore import QThread, QUrl, pyqtSignal, QProcess, Qt, QTimer, QEventLoop, QSize
from PyQt5.QtGui import QDesktopServices, QIcon, QPixmap


def _get_base_dir():
    if getattr(sys, 'frozen', False):
        return Path(sys.executable).parent
    return Path(__file__).parent


def _load_icon():
    base = _get_base_dir()
    for name in ('loreon_app_icon_2.ico', 'loreon_app_icon_2.png'):
        p = base / name
        if p.exists():
            return QIcon(str(p))
    return QIcon()

try:
    from pipeline_worker import PipelineWorker, DOCKER_IMAGE
except ImportError:
    print("ERROR: Cannot find 'pipeline_worker.py'.")
    sys.exit(1)


class DockerBuildDialog(QDialog):
    """Dialog that runs 'docker build' and shows live output."""

    def __init__(self, build_dir, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Building Docker Image")
        self.setMinimumSize(700, 400)
        self.build_dir = build_dir

        layout = QVBoxLayout()
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setStyleSheet(
            "font-family: 'Courier New', monospace; background-color: #2b2b2b; color: #f0f0f0;"
        )
        layout.addWidget(self.log)

        self.buttons = QDialogButtonBox(QDialogButtonBox.Close)
        self.buttons.button(QDialogButtonBox.Close).setEnabled(False)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)
        self.setLayout(layout)

        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.MergedChannels)
        self.process.readyRead.connect(self._on_output)
        self.process.finished.connect(self._on_finished)

    def start_build(self):
        self.log.appendPlainText(f"[INFO] Running: docker build -t {DOCKER_IMAGE} {self.build_dir}\n")
        self.process.start("docker", ["build", "-t", DOCKER_IMAGE, self.build_dir])

    def _on_output(self):
        data = self.process.readAll().data().decode("utf-8", errors="replace")
        self.log.appendPlainText(data)

    def _on_finished(self, exit_code, _exit_status):
        if exit_code == 0:
            self.log.appendPlainText(f"\n[SUCCESS] Image '{DOCKER_IMAGE}' built successfully.")
        else:
            self.log.appendPlainText(f"\n[ERROR] Build failed (exit code {exit_code}).")
        self.buttons.button(QDialogButtonBox.Close).setEnabled(True)


class _ScaledPixmapLabel(QLabel):
    """QLabel that fills the full width up to an optional max_height, and is clickable."""
    clicked = pyqtSignal()

    def __init__(self, pixmap, max_height=120, parent=None):
        super().__init__(parent)
        self._src = pixmap
        self._max_height = max_height
        self.setAlignment(Qt.AlignCenter)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.setCursor(Qt.PointingHandCursor)

    def hasHeightForWidth(self):
        return not self._src.isNull()

    def heightForWidth(self, width):
        if self._src.isNull() or self._src.width() == 0:
            return 0
        natural = int(self._src.height() * width / self._src.width())
        return min(natural, self._max_height)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self._src.isNull():
            return
        available_h = self.heightForWidth(self.width())
        # Scale to fit within (width × max_height), preserving aspect ratio
        scaled = self._src.scaled(
            self.width(), available_h,
            Qt.KeepAspectRatio, Qt.SmoothTransformation
        )
        super().setPixmap(scaled)

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.clicked.emit()
        super().mousePressEvent(event)


class CreditsDialog(QDialog):
    """Dialog showing license and credits read from credits.md."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("License & Credits — LOREON")
        self.setMinimumSize(620, 520)

        layout = QVBoxLayout()
        browser = QTextBrowser()
        browser.setOpenExternalLinks(True)

        credits_path = _get_base_dir() / 'credits.md'
        if credits_path.exists():
            browser.setMarkdown(credits_path.read_text(encoding='utf-8'))
        else:
            browser.setPlainText("Credits file not found.")

        layout.addWidget(browser)

        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self.setLayout(layout)


class MainWindow(QMainWindow):
    start_pipeline_signal = pyqtSignal(dict)
    stop_pipeline_signal = pyqtSignal()

    def __init__(self):
        super().__init__()
        self.setWindowTitle("LOREON Metagenomic Pipeline (OTU Gen v3.3 GUI)")
        self.setGeometry(100, 100, 900, 750)
        self.is_running = False
        self._log_file = None
        icon = _load_icon()
        self.setWindowIcon(icon)
        self.init_ui()
        self.setup_worker_thread()
        self._setup_tray(icon)

    def init_ui(self):
        outer_layout = QVBoxLayout()
        outer_layout.setContentsMargins(0, 0, 0, 0)
        outer_layout.setSpacing(0)

        # --- Logo banner (edge-to-edge, clickable) ---
        logo_path = _get_base_dir() / 'loreon.jpeg'
        if logo_path.exists():
            logo_label = _ScaledPixmapLabel(QPixmap(str(logo_path)), max_height=120)
            logo_label.clicked.connect(self._show_credits)
            outer_layout.addWidget(logo_label)
        # --- End logo banner ---

        inner_widget = QWidget()
        main_layout = QVBoxLayout()
        inner_widget.setLayout(main_layout)
        outer_layout.addWidget(inner_widget)

        path_group = QGroupBox("1. Paths (Mandatory)")
        path_layout = QGridLayout()
        self.input_dir_btn = QPushButton("Select Input folder...")
        self.input_dir_label = QLineEdit()
        self.input_dir_label.setPlaceholderText("e.g., /path/to/raw/data (barcodeXX, unclassified)")
        self.db_btn = QPushButton("Ref DB FASTA...")
        self.db_label = QLineEdit()
        self.db_label.setPlaceholderText("e.g., /path/to/unite.fasta")
        self.output_dir_btn = QPushButton("Select Output Folder...")
        self.output_dir_label = QLineEdit()
        self.output_dir_label.setPlaceholderText("e.g., /path/to/pipeline_results")
        path_layout.addWidget(self.input_dir_btn, 0, 0)
        path_layout.addWidget(self.input_dir_label, 0, 1)
        path_layout.addWidget(self.db_btn, 1, 0)
        path_layout.addWidget(self.db_label, 1, 1)
        path_layout.addWidget(self.output_dir_btn, 2, 0)
        path_layout.addWidget(self.output_dir_label, 2, 1)
        path_group.setLayout(path_layout)
        main_layout.addWidget(path_group)

        options_group = QGroupBox("2. Pipeline Parameters")
        options_layout = QHBoxLayout()
        options_left_layout = QVBoxLayout()
        self.format_combo = QComboBox()
        self.format_combo.addItems(['unite', 'silva', 'eukariome', 'cbs', 'none'])
        self.format_combo.setCurrentText('unite')
        options_left_layout.addWidget(QLabel("Database Format:"))
        options_left_layout.addWidget(self.format_combo)
        self.debug_check = QCheckBox("Debug Mode (verbose)")
        options_left_layout.addWidget(self.debug_check)
        self.open_report_check = QCheckBox("Open HTML report when done")
        self.open_report_check.setChecked(True)
        options_left_layout.addWidget(self.open_report_check)
        self.force_tax_map_check = QCheckBox("Force Taxonomy Map Recreation")
        self.force_tax_map_check.setToolTip(
            "If checked, ignores cache and recreates the taxonomy map.\n"
            "Useful if the database has been updated."
        )
        options_left_layout.addWidget(self.force_tax_map_check)
        self.profile_check = QCheckBox("Enable Performance Profiling")
        self.profile_check.setToolTip(
            "Saves profile_report.json in the output folder with wall time,\n"
            "CPU time, memory and domain metrics for each pipeline step."
        )
        options_left_layout.addWidget(self.profile_check)
        self.delete_bam_check = QCheckBox("Delete BAM/SAM after analysis")
        self.delete_bam_check.setToolTip(
            "Delete intermediate BAM and SAM files after tabeling\n"
            "to free disk space. The OTU results are NOT affected."
        )
        options_left_layout.addWidget(self.delete_bam_check)
        options_layout.addLayout(options_left_layout)

        filter_box = QGroupBox("Read & Alignment Filter")
        filter_layout = QGridLayout()
        self.filter_check = QCheckBox("Enable Length Filter")
        self.filter_check.setChecked(True)
        self.min_len_input = QSpinBox()
        self.min_len_input.setRange(0, 99999)
        self.min_len_input.setValue(200)
        self.max_len_input = QSpinBox()
        self.max_len_input.setRange(0, 99999)
        self.max_len_input.setValue(300)
        self.min_identity_input = QDoubleSpinBox()
        self.min_identity_input.setRange(0.0, 100.0)
        self.min_identity_input.setSingleStep(0.5)
        self.min_identity_input.setDecimals(1)
        self.min_identity_input.setValue(95.0)
        self.min_identity_input.setToolTip(
            "Minimum percent identity for mapped reads (0 = disabled).\n"
            "E.g. 97.0 keeps only alignments with ≥97% identity."
        )
        self.min_ref_cov_input = QDoubleSpinBox()
        self.min_ref_cov_input.setRange(0.0, 100.0)
        self.min_ref_cov_input.setSingleStep(1.0)
        self.min_ref_cov_input.setDecimals(1)
        self.min_ref_cov_input.setValue(90.0)
        self.min_ref_cov_input.setToolTip(
            "Minimum reference coverage for mapped reads (0 = disabled).\n"
            "E.g. 80.0 keeps only alignments covering ≥80% of the reference."
        )
        filter_layout.addWidget(self.filter_check, 0, 0, 1, 2)
        filter_layout.addWidget(QLabel("Min Len:"), 1, 0)
        filter_layout.addWidget(self.min_len_input, 1, 1)
        filter_layout.addWidget(QLabel("Max Len:"), 2, 0)
        filter_layout.addWidget(self.max_len_input, 2, 1)
        filter_layout.addWidget(QLabel("Min Identity (%):"), 3, 0)
        filter_layout.addWidget(self.min_identity_input, 3, 1)
        filter_layout.addWidget(QLabel("Min Ref Cov (%):"), 4, 0)
        filter_layout.addWidget(self.min_ref_cov_input, 4, 1)
        filter_box.setLayout(filter_layout)
        options_layout.addWidget(filter_box)

        perf_box = QGroupBox("Performance (Threads)")
        perf_layout = QGridLayout()
        self.total_threads_input = QSpinBox()
        self.total_threads_input.setRange(1, multiprocessing.cpu_count())
        self.total_threads_input.setValue(8)
        self.job_threads_input = QSpinBox()
        self.job_threads_input.setRange(1, multiprocessing.cpu_count())
        self.job_threads_input.setValue(2)
        self.kmer_input = QSpinBox()
        self.kmer_input.setValue(15)
        self.window_input = QSpinBox()
        self.window_input.setValue(10)
        perf_layout.addWidget(QLabel("Total Threads (-T):"), 0, 0)
        perf_layout.addWidget(self.total_threads_input, 0, 1)
        perf_layout.addWidget(QLabel("Threads/Job (-t):"), 1, 0)
        perf_layout.addWidget(self.job_threads_input, 1, 1)
        perf_layout.addWidget(QLabel("K-mer (-k):"), 2, 0)
        perf_layout.addWidget(self.kmer_input, 2, 1)
        perf_layout.addWidget(QLabel("Window (-w):"), 3, 0)
        perf_layout.addWidget(self.window_input, 3, 1)
        perf_box.setLayout(perf_layout)
        options_layout.addWidget(perf_box)
        options_group.setLayout(options_layout)
        main_layout.addWidget(options_group)

        # --- Docker section ---
        docker_group = QGroupBox("3. Docker (required on Windows)")
        docker_layout = QHBoxLayout()

        self.docker_check = QCheckBox("Use Docker  (minimap2 + samtools run inside container)")
        if sys.platform == 'win32':
            self.docker_check.setChecked(True)
            self.docker_check.setToolTip("Docker is required on Windows and cannot be disabled.")
        docker_layout.addWidget(self.docker_check)

        self.docker_status_label = QLabel("Status: unknown")
        docker_layout.addWidget(self.docker_status_label)

        docker_layout.addStretch()

        self.docker_check_btn = QPushButton("Check Image")
        self.docker_check_btn.setToolTip(f"Verify that the Docker image '{DOCKER_IMAGE}' exists")
        self.docker_check_btn.clicked.connect(self.check_docker_image)
        docker_layout.addWidget(self.docker_check_btn)

        self.docker_build_btn = QPushButton("Build Image")
        self.docker_build_btn.setToolTip(
            f"Build the Docker image '{DOCKER_IMAGE}' from the Dockerfile in the app folder"
        )
        self.docker_build_btn.clicked.connect(self.build_docker_image)
        docker_layout.addWidget(self.docker_build_btn)

        docker_group.setLayout(docker_layout)
        main_layout.addWidget(docker_group)
        # --- End Docker section ---

        button_layout = QHBoxLayout()
        self.run_btn = QPushButton("START PIPELINE")
        self.run_btn.setStyleSheet(
            "font-size: 18px; background-color: #2ca02c; color: white; font-weight: bold; padding: 10px;")
        self.stop_btn = QPushButton("STOP")
        self.stop_btn.setStyleSheet(
            "font-size: 18px; background-color: #d62728; color: white; font-weight: bold; padding: 10px;")
        self.stop_btn.setEnabled(False)
        button_layout.addWidget(self.run_btn)
        button_layout.addWidget(self.stop_btn)
        main_layout.addLayout(button_layout)

        log_group = QGroupBox("Execution Log")
        log_layout = QVBoxLayout()
        self.log_output = QPlainTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setStyleSheet(
            "font-family: 'Courier New', monospace; background-color: #2b2b2b; color: #f0f0f0;")
        self.status_bar = QProgressBar()
        self.status_bar.setRange(0, 100)
        self.status_bar.setValue(0)
        self.status_bar.setTextVisible(True)
        self.status_bar.setFormat("Waiting...")
        log_layout.addWidget(self.log_output)
        log_layout.addWidget(self.status_bar)
        log_group.setLayout(log_layout)
        main_layout.addWidget(log_group)

        self.input_dir_btn.clicked.connect(lambda: self.select_directory(self.input_dir_label))
        self.output_dir_btn.clicked.connect(lambda: self.select_directory(self.output_dir_label))
        self.db_btn.clicked.connect(lambda: self.select_file(self.db_label, "Database FASTA (*.fasta *.fa)"))
        self.filter_check.stateChanged.connect(self.toggle_filter_inputs)
        self.run_btn.clicked.connect(self.run_pipeline)
        self.stop_btn.clicked.connect(self.stop_pipeline)

        central_widget = QWidget()
        central_widget.setLayout(outer_layout)
        self.setCentralWidget(central_widget)

    def _show_credits(self):
        CreditsDialog(self).exec_()

    def setup_worker_thread(self):
        self.worker_thread = QThread()
        self.worker = PipelineWorker()
        self.worker.moveToThread(self.worker_thread)
        self.start_pipeline_signal.connect(self.worker.run)
        self.stop_pipeline_signal.connect(self.worker.stop)
        self.worker.log_signal.connect(self.append_log)
        self.worker.progress_signal.connect(self.update_progress)
        self.worker.finished_signal.connect(self.pipeline_finished)
        self.worker_thread.start()

    def select_directory(self, label_widget):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Folder")
        if dir_path:
            label_widget.setText(dir_path)

    def select_file(self, label_widget, file_filter):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", "", file_filter)
        if file_path:
            label_widget.setText(file_path)

    def toggle_filter_inputs(self, state):
        is_enabled = (state == 2)
        self.min_len_input.setEnabled(is_enabled)
        self.max_len_input.setEnabled(is_enabled)

    def append_log(self, text):
        stamped = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {text}"
        # OPT-13: appendPlainText auto-scrolls when cursor is at bottom;
        # explicit setValue(maximum()) was redundant and caused GUI lag
        self.log_output.appendPlainText(stamped)
        if self._log_file and not self._log_file.closed:
            try:
                self._log_file.write(stamped + "\n")
                self._log_file.flush()
            except Exception:
                pass

    def update_progress(self, percentage, text):
        self.status_bar.setValue(percentage)
        self.status_bar.setFormat(text)

    # --- Docker helpers ---

    def _docker_available(self):
        return shutil.which("docker") is not None

    def check_docker_image(self):
        if not self._docker_available():
            self.docker_status_label.setText("Status: Docker NOT found")
            QMessageBox.warning(self, "Docker not found",
                                "Docker is not installed or not in PATH.\n"
                                "Install Docker Desktop: https://www.docker.com/products/docker-desktop")
            return
        try:
            result = subprocess.run(
                ["docker", "images", "-q", DOCKER_IMAGE],
                capture_output=True, text=True, timeout=10
            )
            if result.stdout.strip():
                self.docker_status_label.setText(f"Status: image '{DOCKER_IMAGE}' ready")
                QMessageBox.information(self, "Docker OK", f"Image '{DOCKER_IMAGE}' found and ready.")
            else:
                self.docker_status_label.setText("Status: image NOT found")
                QMessageBox.warning(self, "Image not found",
                                    f"Image '{DOCKER_IMAGE}' not found.\n"
                                    "Click 'Build Image' to create it.")
        except Exception as e:
            self.docker_status_label.setText("Status: error")
            QMessageBox.critical(self, "Error", f"Cannot check Docker image:\n{e}")

    def build_docker_image(self):
        if not self._docker_available():
            QMessageBox.warning(self, "Docker not found",
                                "Docker is not installed or not in PATH.")
            return
        if getattr(sys, 'frozen', False):
            build_dir = str(Path(sys.executable).parent)
        else:
            build_dir = str(Path(__file__).parent)
        dlg = DockerBuildDialog(build_dir, parent=self)
        dlg.show()
        dlg.start_build()
        dlg.exec_()
        self.check_docker_image()

    # --- Pipeline execution ---

    def run_pipeline(self):
        if self.is_running:
            QMessageBox.warning(self, "Warning", "Pipeline is already running.")
            return

        use_docker = self.docker_check.isChecked()
        if sys.platform == 'win32' and not use_docker:
            QMessageBox.critical(self, "Docker required",
                                 "Docker must be enabled on Windows, please enable it.")
            return
        if use_docker and not self._docker_available():
            QMessageBox.critical(self, "Docker not found",
                                 "Docker mode is enabled but Docker is not in PATH.\n"
                                 "Install Docker Desktop or disable Docker mode.")
            return

        self.log_output.clear()

        # Open log file in the output directory
        try:
            log_path = Path(self.output_dir_label.text()) / "pipeline.log"
            self._log_file = open(log_path, "w", encoding="utf-8")
        except Exception as e:
            self._log_file = None
            self.log_output.appendPlainText(f"[WARNING] Cannot create pipeline.log: {e}")

        self.append_log("--- STARTING PIPELINE ---")
        if use_docker:
            self.append_log(f"[INFO] Docker mode: pipeline will run inside '{DOCKER_IMAGE}'")

        settings = {
            "input_dir": self.input_dir_label.text(),
            "output_dir": self.output_dir_label.text(),
            "db_path": self.db_label.text(),
            "format": self.format_combo.currentText(),
            "debug": self.debug_check.isChecked(),
            "enable_filter": self.filter_check.isChecked(),
            "min_len": self.min_len_input.value(),
            "max_len": self.max_len_input.value(),
            "min_percent_identity": self.min_identity_input.value(),
            "min_ref_coverage": self.min_ref_cov_input.value(),
            "total_threads": self.total_threads_input.value(),
            "threads_per_job": self.job_threads_input.value(),
            "kmer_size": self.kmer_input.value(),
            "window_size": self.window_input.value(),
            "open_report": self.open_report_check.isChecked(),
            "force_tax_map": self.force_tax_map_check.isChecked(),
            "profile": self.profile_check.isChecked(),
            "delete_bam": self.delete_bam_check.isChecked(),
            "use_docker": use_docker,
        }
        if not all([settings["input_dir"], settings["output_dir"], settings["db_path"]]):
            QMessageBox.critical(self, "Error", "Input, Output and Database paths are mandatory.")
            self.append_log("--- ERROR: Missing paths. Pipeline aborted. ---")
            return
        self.set_running_state(True)
        self.start_pipeline_signal.emit(settings)

    def stop_pipeline(self):
        if not self.is_running:
            return
        response = QMessageBox.question(self, "Confirm Interruption",
                                        "Are you sure you want to stop the pipeline?",
                                        QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if response == QMessageBox.Yes:
            self.append_log("--- USER REQUESTED INTERRUPTION ---")
            self.status_bar.setFormat("Stopping...")
            self.stop_pipeline_signal.emit()

    def pipeline_finished(self, success, html_report_path):
        self.append_log("--- PIPELINE FINISHED ---")

        # Close log file
        if self._log_file and not self._log_file.closed:
            try:
                self._log_file.close()
            except Exception:
                pass
            self._log_file = None

        self.set_running_state(False)
        if success:
            self.status_bar.setFormat("Completed successfully!")
            self.status_bar.setValue(100)
            if self.open_report_check.isChecked() and html_report_path:
                self.append_log(f"Opening report: {html_report_path}")
                try:
                    QDesktopServices.openUrl(QUrl.fromLocalFile(html_report_path))
                except Exception as e:
                    self.append_log(f"Cannot open report automatically: {e}")
            elif self.open_report_check.isChecked():
                self.append_log("HTML report not generated (check log).")
        else:
            if "Stopping" in self.status_bar.format():
                self.status_bar.setFormat("Process stopped by user.")
            else:
                self.status_bar.setFormat("Failed! Check log.")
                QMessageBox.critical(self, "Pipeline Error",
                                     "Execution failed. Check the log for details.")
            self.status_bar.setValue(0)

    def set_running_state(self, is_running):
        self.is_running = is_running
        self.run_btn.setEnabled(not is_running)
        self.stop_btn.setEnabled(is_running)
        self.run_btn.setText("RUNNING..." if is_running else "START PIPELINE")
        for widget in [self.input_dir_btn, self.input_dir_label, self.db_btn,
                       self.db_label, self.output_dir_btn, self.output_dir_label,
                       self.format_combo, self.debug_check,
                       self.open_report_check, self.filter_check, self.min_len_input,
                       self.max_len_input, self.min_identity_input, self.min_ref_cov_input,
                       self.total_threads_input,
                       self.job_threads_input, self.kmer_input, self.window_input,
                       self.force_tax_map_check, self.profile_check,
                       self.delete_bam_check, self.docker_check,
                       self.docker_check_btn, self.docker_build_btn]:
            widget.setEnabled(not is_running)

    def _setup_tray(self, icon):
        self.tray_icon = QSystemTrayIcon(icon, self)
        tray_menu = QMenu()
        show_action = tray_menu.addAction("Mostra")
        show_action.triggered.connect(self.show)
        quit_action = tray_menu.addAction("Esci")
        quit_action.triggered.connect(QApplication.quit)
        self.tray_icon.setContextMenu(tray_menu)
        self.tray_icon.setToolTip("LOREON Pipeline")
        self.tray_icon.activated.connect(self._on_tray_activated)
        self.tray_icon.show()

    def _on_tray_activated(self, reason):
        if reason == QSystemTrayIcon.DoubleClick:
            self.show()
            self.raise_()
            self.activateWindow()

    def closeEvent(self, event):
        if self.is_running:
            response = QMessageBox.question(self, "Confirm Exit",
                                            "Pipeline is still running. Exit anyway?",
                                            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if response == QMessageBox.Yes:
                self.stop_pipeline_signal.emit()
                self.worker_thread.quit()
                self.worker_thread.wait(5000)
                event.accept()
            else:
                event.ignore()
        else:
            self.worker_thread.quit()
            self.worker_thread.wait(5000)
            event.accept()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = QApplication(sys.argv)

    icon = _load_icon()
    app.setWindowIcon(icon)

    jpeg_path = _get_base_dir() / 'loreon.jpeg'
    if jpeg_path.exists():
        pixmap = QPixmap(str(jpeg_path)).scaledToWidth(
            500, Qt.SmoothTransformation
        )
        splash = QSplashScreen(pixmap, Qt.WindowStaysOnTopHint)
        splash.show()
        app.processEvents()
        loop = QEventLoop()
        QTimer.singleShot(1000, loop.quit)
        loop.exec_()
        splash.close()

    window = MainWindow()
    window.show()

    sys.exit(app.exec_())
