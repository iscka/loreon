"""Graphical interface defaults and the worker that launches the pipeline."""
import os
import sys
import threading
import time

import pytest

pytest.importorskip('PyQt5')

from PyQt5.QtWidgets import QApplication, QDialog, QMessageBox

import OtuUtils
import pipeline_gui
import pipeline_worker
import venv_manager
from helpers import SYSTEM_PATH, make_tool


@pytest.fixture(scope='session')
def app():
    instance = QApplication.instance() or QApplication([])
    yield instance


@pytest.fixture
def dialogs(monkeypatch):
    """Replace modal message boxes, which would wait for a user, with recorders."""
    shown = []
    for kind in ('warning', 'critical', 'information'):
        monkeypatch.setattr(QMessageBox, kind,
                            lambda *args, _kind=kind, **kwargs: shown.append((_kind, args[2])))
    monkeypatch.setattr(QMessageBox, 'question', lambda *args, **kwargs: QMessageBox.Yes)
    return shown


@pytest.fixture
def window(app, dialogs):
    main_window = pipeline_gui.MainWindow()
    yield main_window
    main_window.is_running = False
    main_window.close()
    app.processEvents()


def base_settings(**overrides):
    settings = {
        'input_dir': '/data/in', 'output_dir': '/data/out', 'db_path': '/data/db/ref.fasta',
        'format': 'unite', 'debug': False, 'enable_filter': False, 'min_len': 0, 'max_len': 0,
        'total_threads': 2, 'threads_per_job': 1, 'kmer_size': 15, 'window_size': 10,
        'use_docker': False,
    }
    settings.update(overrides)
    return settings


def option_value(arguments, flag):
    return arguments[arguments.index(flag) + 1]


def test_quality_filters_default_to_zero(window):
    assert window.min_identity_input.value() == 0.0
    assert window.min_ref_cov_input.value() == 0.0


def test_length_filter_is_off_by_default(window):
    settings = window.collect_settings()
    assert settings['enable_filter'] is False
    assert not window.min_len_input.isEnabled()
    assert not window.max_len_input.isEnabled()


def test_length_menu_lists_every_preset(window):
    combo = window.length_mode_combo
    entries = [combo.itemData(i) for i in range(combo.count())]
    assert entries[:2] == [None, 'custom']
    assert entries[2:] == list(OtuUtils.LENGTH_PRESETS)


def test_choosing_a_preset_fills_the_length_window(window):
    window.length_mode_combo.setCurrentIndex(window.length_mode_combo.findData('its-lsu'))
    settings = window.collect_settings()
    assert settings['enable_filter'] is True
    assert (settings['min_len'], settings['max_len']) == (900, 1800)


def test_editing_a_preset_value_switches_to_custom_range(window):
    window.length_mode_combo.setCurrentIndex(window.length_mode_combo.findData('its-lsu'))
    window.max_len_input.setValue(2000)
    assert window.length_mode_combo.currentData() == 'custom'
    assert window.collect_settings()['max_len'] == 2000


def test_sample_time_limit_defaults_to_no_limit(window):
    assert window.collect_settings()['sample_timeout_min'] == 0


def test_worker_disables_quality_filters_when_settings_omit_them(app):
    worker = pipeline_worker.PipelineWorker()
    worker.settings = base_settings()
    arguments = worker._pipeline_args('in', 'out', 'db.fasta')
    assert float(option_value(arguments, '--min-percent-identity')) == 0.0
    assert float(option_value(arguments, '--min-ref-coverage')) == 0.0


def test_worker_passes_the_sample_time_limit_only_when_set(app):
    worker = pipeline_worker.PipelineWorker()
    worker.settings = base_settings(sample_timeout_min=45)
    assert option_value(worker._pipeline_args('in', 'out', 'db.fasta'), '--sample-timeout') == '45'
    worker.settings = base_settings(sample_timeout_min=0)
    assert '--sample-timeout' not in worker._pipeline_args('in', 'out', 'db.fasta')


def test_retain_unmapped_reads_is_one_checkbox_ticked_by_default(window):
    assert window.retain_unmapped_check.text() == 'Retain unmapped reads'
    assert window.retain_unmapped_check.isChecked()
    assert window.collect_settings()['retain_unmapped'] is True
    window.retain_unmapped_check.setChecked(False)
    assert window.collect_settings()['retain_unmapped'] is False


def test_profiling_and_bam_removal_are_ticked_by_default(window):
    settings = window.collect_settings()
    assert window.profile_check.isChecked() and settings['profile'] is True
    assert window.delete_bam_check.isChecked() and settings['delete_bam'] is True


def test_retain_unmapped_reads_cannot_be_changed_during_a_run(window):
    window.set_running_state(True)
    assert not window.retain_unmapped_check.isEnabled()
    window.set_running_state(False)
    assert window.retain_unmapped_check.isEnabled()


def test_worker_passes_retain_unmapped_only_when_set(app):
    worker = pipeline_worker.PipelineWorker()
    worker.settings = base_settings(retain_unmapped=True)
    assert '--retain-unmapped' in worker._pipeline_args('in', 'out', 'db.fasta')
    worker.settings = base_settings()
    assert '--retain-unmapped' not in worker._pipeline_args('in', 'out', 'db.fasta')


def test_docker_mode_receives_retain_unmapped_too(app):
    worker = pipeline_worker.PipelineWorker()
    worker.settings = base_settings(retain_unmapped=True, use_docker=True)
    assert '--retain-unmapped' in worker._build_docker_pipeline_cmd()


def test_each_container_gets_its_own_name(app):
    worker = pipeline_worker.PipelineWorker()
    worker.settings = base_settings(use_docker=True)
    first = option_value(worker._docker_base_cmd(), '--name')
    second = option_value(worker._docker_base_cmd(), '--name')
    assert first.startswith('loreon-')
    assert first != second


@pytest.mark.fake_tools
def test_stop_kills_the_running_container(app, tmp_path, monkeypatch):
    bin_dir = tmp_path / 'bin'
    log = tmp_path / 'docker.log'
    make_tool(bin_dir, 'docker', f"""
        echo "$@" >> "{log}"
        if [ "$1" = "run" ]; then
            echo "--- STEP 4: Mapping and Tabeling ---"
            exec sleep 300
        fi
        exit 0
    """)
    monkeypatch.setenv('PATH', f"{bin_dir}:{SYSTEM_PATH}")
    worker = pipeline_worker.PipelineWorker()
    worker.settings = base_settings(use_docker=True)
    worker.is_running = True
    command = worker._build_docker_pipeline_cmd()
    name = option_value(command, '--name')
    outcome = {}
    runner = threading.Thread(target=lambda: outcome.update(result=worker.run_process(command)))
    runner.start()
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline and not (log.exists() and 'run' in log.read_text()):
        time.sleep(0.05)
    started = time.monotonic()
    worker.stop()
    assert time.monotonic() - started < 2
    runner.join(10)
    assert not runner.is_alive()
    assert outcome['result'] is False
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline and f"kill {name}" not in log.read_text():
        time.sleep(0.05)
    assert f"kill {name}" in log.read_text()


def run_worker_with_exit_codes(app, monkeypatch, pipeline_code, report_code=0):
    worker = pipeline_worker.PipelineWorker()
    monkeypatch.setattr(worker, '_build_native_pipeline_cmd',
                        lambda: ['bash', '-c', f'echo pipeline output; exit {pipeline_code}'])
    monkeypatch.setattr(worker, '_build_native_report_cmd',
                        lambda: ['bash', '-c', f'exit {report_code}'])
    events = {'finished': [], 'warnings': []}
    worker.finished_signal.connect(lambda ok, path: events['finished'].append(ok))
    worker.warning_signal.connect(lambda text: events['warnings'].append(text))
    worker.run(base_settings())
    return events


@pytest.mark.fake_tools
def test_worker_reports_partial_failure_and_still_builds_the_report(app, monkeypatch):
    events = run_worker_with_exit_codes(app, monkeypatch, pipeline_worker.PARTIAL_FAILURE_EXIT_CODE)
    assert events['finished'] == [True]
    assert len(events['warnings']) == 1


@pytest.mark.fake_tools
def test_worker_reports_failure_when_the_pipeline_fails(app, monkeypatch):
    events = run_worker_with_exit_codes(app, monkeypatch, 1)
    assert events['finished'] == [False]
    assert events['warnings'] == []


def test_window_tells_the_user_when_some_samples_failed(window, dialogs):
    window.set_running_state(True)
    window.on_worker_warning('2 of 12 samples failed')
    window.pipeline_finished(True, None)
    assert ('warning', '2 of 12 samples failed') in [(kind, text[:22]) for kind, text in dialogs]
    assert 'error' in window.status_bar.format().lower()


def write_reads(folder, length, count=40):
    folder.mkdir(parents=True)
    with open(folder / 'reads.fastq', 'w') as handle:
        for i in range(count):
            handle.write(f"@r{i}\n{'A' * length}\n+\n{'I' * length}\n")


def test_scan_shows_the_lengths_and_applies_the_suggested_window(window, tmp_path, monkeypatch):
    write_reads(tmp_path / 'barcode01', 1210)
    write_reads(tmp_path / 'barcode02', 1390)
    reports = []

    def accept(dialog):
        reports.append(dialog.report_text)
        return QDialog.Accepted

    monkeypatch.setattr(pipeline_gui.ReadLengthDialog, 'exec_', accept)
    window.input_dir_label.setText(str(tmp_path))
    window.scan_read_lengths()
    assert 'median' in reports[0]
    settings = window.collect_settings()
    assert window.length_mode_combo.currentData() == 'custom'
    assert (settings['min_len'], settings['max_len']) == (1200, 1400)


def test_scan_leaves_the_settings_alone_when_the_dialog_is_closed(window, tmp_path, monkeypatch):
    write_reads(tmp_path / 'barcode01', 1210)
    monkeypatch.setattr(pipeline_gui.ReadLengthDialog, 'exec_', lambda dialog: QDialog.Rejected)
    window.input_dir_label.setText(str(tmp_path))
    window.scan_read_lengths()
    assert window.collect_settings()['enable_filter'] is False


def test_scan_asks_for_the_input_folder_when_it_is_missing(window, dialogs):
    window.input_dir_label.setText('')
    window.scan_read_lengths()
    assert dialogs and dialogs[0][0] == 'warning'


def pump_until(app, condition, seconds=10):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline and not condition():
        app.processEvents()
        time.sleep(0.05)
    app.processEvents()
    return condition()


@pytest.mark.fake_tools
def test_stop_during_the_report_step_is_not_reported_as_success(app, monkeypatch):
    worker = pipeline_worker.PipelineWorker()
    monkeypatch.setattr(worker, '_build_native_pipeline_cmd', lambda: ['bash', '-c', 'exit 0'])
    monkeypatch.setattr(worker, '_build_native_report_cmd',
                        lambda: ['bash', '-c', 'echo report started; exec sleep 30'])
    logs = []
    finished = []
    worker.log_signal.connect(logs.append)
    worker.finished_signal.connect(lambda ok, path: finished.append(ok))
    runner = threading.Thread(target=worker.run, args=(base_settings(),))
    runner.start()
    assert pump_until(app, lambda: 'report started' in logs)
    worker.stop()
    assert pump_until(app, lambda: not runner.is_alive())
    assert pump_until(app, lambda: bool(finished))
    assert finished == [False]


@pytest.mark.fake_tools
def test_process_started_after_a_stop_request_is_terminated(app, monkeypatch):
    import subprocess
    worker = pipeline_worker.PipelineWorker()
    worker.settings = base_settings()
    worker.is_running = True
    real_popen = subprocess.Popen

    def popen_then_stop(*args, **kwargs):
        process = real_popen(*args, **kwargs)
        worker.is_running = False
        return process

    monkeypatch.setattr(pipeline_worker.subprocess, 'Popen', popen_then_stop)
    started = time.monotonic()
    result = worker.run_process(['bash', '-c', 'echo started; exec sleep 30'])
    assert result is False
    assert time.monotonic() - started < 10


@pytest.mark.fake_tools
def test_container_stop_can_be_awaited(app, tmp_path, monkeypatch):
    bin_dir = tmp_path / 'bin'
    log = tmp_path / 'docker.log'
    make_tool(bin_dir, 'docker', f"""
        if [ "$1" = "run" ]; then
            echo "run" >> "{log}"
            echo "--- STEP 4: Mapping and Tabeling ---"
            exec sleep 300
        fi
        if [ "$1" = "kill" ]; then
            sleep 1
            echo "killed $2" >> "{log}"
        fi
        exit 0
    """)
    monkeypatch.setenv('PATH', f"{bin_dir}:{SYSTEM_PATH}")
    worker = pipeline_worker.PipelineWorker()
    worker.settings = base_settings(use_docker=True)
    worker.is_running = True
    command = worker._build_docker_pipeline_cmd()
    name = option_value(command, '--name')
    runner = threading.Thread(target=worker.run_process, args=(command,))
    runner.start()
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline and not (log.exists() and 'run' in log.read_text()):
        time.sleep(0.05)
    worker.stop()
    assert worker.wait_for_container_stop(10) is True
    assert f"killed {name}" in log.read_text()
    runner.join(10)
    assert not runner.is_alive()


def test_waiting_for_a_container_stop_returns_at_once_when_nothing_is_pending(app):
    worker = pipeline_worker.PipelineWorker()
    started = time.monotonic()
    assert worker.wait_for_container_stop(10) is True
    assert time.monotonic() - started < 1


def test_closing_the_window_during_a_run_waits_for_the_container_to_stop(window, monkeypatch):
    waited = []
    monkeypatch.setattr(window.worker, 'wait_for_container_stop',
                        lambda timeout=None: waited.append(timeout) or True)
    window.is_running = True
    window.close()
    assert len(waited) == 1


def test_start_is_refused_with_instructions_when_components_are_missing(window, dialogs, tmp_path, monkeypatch):
    monkeypatch.setattr(pipeline_gui, 'environment_problems',
                        lambda use_docker: ['samtools was not found: install it with conda'])
    window.input_dir_label.setText(str(tmp_path))
    window.output_dir_label.setText(str(tmp_path))
    window.db_label.setText(str(tmp_path / 'db.fasta'))
    window.run_pipeline()
    assert window.is_running is False
    assert dialogs and dialogs[-1][0] == 'critical'
    assert 'samtools was not found' in dialogs[-1][1]
    assert not (tmp_path / 'pipeline.log').exists()


def test_window_reports_missing_components_when_it_opens(app, dialogs, monkeypatch):
    monkeypatch.setattr(pipeline_gui, 'environment_problems',
                        lambda use_docker: ['minimap2 was not found: install it with conda'])
    main_window = pipeline_gui.MainWindow()
    try:
        assert 'minimap2 was not found' in main_window.log_output.toPlainText()
    finally:
        main_window.close()
        app.processEvents()


def test_example_button_fills_the_form_and_starts(window, tmp_path, monkeypatch):
    started = []
    monkeypatch.setattr(window, 'run_pipeline', lambda: started.append(window.collect_settings()))
    monkeypatch.setattr(pipeline_gui, 'default_example_output_dir', lambda: tmp_path / 'example_out')
    window.length_mode_combo.setCurrentIndex(window.length_mode_combo.findData('its1'))
    window.min_identity_input.setValue(97.0)
    window.example_btn.click()
    assert len(started) == 1
    settings = started[0]
    assert settings['input_dir'] == str(OtuUtils.EXAMPLE_READS)
    assert settings['db_path'] == str(OtuUtils.EXAMPLE_DATABASE)
    assert settings['output_dir'] == str(tmp_path / 'example_out')
    assert (tmp_path / 'example_out').is_dir()
    assert settings['format'] == 'unite'
    assert settings['enable_filter'] is False
    assert settings['min_percent_identity'] == 0.0 and settings['min_ref_coverage'] == 0.0


def test_example_run_writes_its_taxonomy_map_to_the_output_folder(window, tmp_path, monkeypatch):
    sent = []
    window.start_pipeline_signal.disconnect()
    window.start_pipeline_signal.connect(sent.append)
    monkeypatch.setattr(pipeline_gui, 'environment_problems', lambda use_docker: [])
    monkeypatch.setattr(pipeline_gui, 'default_example_output_dir', lambda: tmp_path / 'example_out')
    monkeypatch.setattr(window, '_docker_available', lambda: True)
    window.docker_check.setChecked(True)
    window.example_btn.click()
    assert len(sent) == 1
    assert sent[0]['force_tax_map'] is True
    assert window.force_tax_map_check.isChecked() is False
    window.pipeline_finished(False, None)
    window.run_pipeline()
    assert len(sent) == 2
    assert sent[1]['force_tax_map'] is False


def finish_example_run(window, tmp_path, table_text):
    window.output_dir_label.setText(str(tmp_path))
    window.db_label.setText(str(OtuUtils.EXAMPLE_DATABASE))
    (tmp_path / f"OTU_Table_{tmp_path.name}_unite_subset.tsv").write_text(table_text)
    window._example_run = True
    window.set_running_state(True)
    window.pipeline_finished(True, None)


def test_matching_example_result_is_announced(window, dialogs, tmp_path):
    finish_example_run(window, tmp_path, OtuUtils.EXAMPLE_EXPECTED_TABLE.read_text())
    assert ('information', True) in [(kind, 'installation works' in text) for kind, text in dialogs]
    assert window._example_run is False


def test_differing_example_result_is_announced(window, dialogs, tmp_path):
    lines = OtuUtils.EXAMPLE_EXPECTED_TABLE.read_text().splitlines()
    finish_example_run(window, tmp_path, '\n'.join(lines[:-1]) + '\n')
    assert ('warning', True) in [(kind, 'differs' in text) for kind, text in dialogs]


def test_installed_windows_application_does_not_create_a_virtual_environment(app, monkeypatch):
    used = []
    monkeypatch.setattr(sys, 'frozen', True, raising=False)
    monkeypatch.setattr(sys, 'platform', 'win32')
    monkeypatch.setattr(venv_manager, 'is_venv_ready', lambda: used.append('check') or False)
    monkeypatch.setattr(venv_manager, 'create_venv', lambda progress_callback=None: used.append('create'))
    assert pipeline_gui._ensure_venv(app) is True
    assert used == []


def test_installed_application_still_creates_it_on_other_systems(app, monkeypatch):
    used = []
    monkeypatch.setattr(sys, 'frozen', True, raising=False)
    monkeypatch.setattr(sys, 'platform', 'linux')
    monkeypatch.setattr(venv_manager, 'is_venv_ready', lambda: used.append('check') or False)
    monkeypatch.setattr(venv_manager, 'create_venv', lambda progress_callback=None: used.append('create'))
    assert pipeline_gui._ensure_venv(app) is True
    assert used == ['check', 'create']


def test_packaged_macos_application_adds_the_tool_folders_to_path(monkeypatch, tmp_path):
    folders = (str(tmp_path / 'homebrew'), str(tmp_path / 'local'))
    for folder in folders:
        os.mkdir(folder)
    monkeypatch.setattr(venv_manager, 'MACOS_TOOL_DIRS', folders)
    monkeypatch.setattr(sys, 'frozen', True, raising=False)
    monkeypatch.setattr(sys, 'platform', 'darwin')
    monkeypatch.setenv('PATH', '/usr/bin:/bin')
    pipeline_gui._add_tool_dirs_to_path()
    assert os.environ['PATH'] == ':'.join(('/usr/bin', '/bin') + folders)


def test_application_run_from_the_source_code_keeps_its_path(monkeypatch, tmp_path):
    monkeypatch.setattr(venv_manager, 'MACOS_TOOL_DIRS', (str(tmp_path),))
    monkeypatch.setattr(sys, 'platform', 'darwin')
    monkeypatch.setenv('PATH', '/usr/bin:/bin')
    pipeline_gui._add_tool_dirs_to_path()
    assert os.environ['PATH'] == '/usr/bin:/bin'


@pytest.fixture
def opened(monkeypatch):
    """Record the files that the window asks the system to open, instead of opening them."""
    files = []
    monkeypatch.setattr(pipeline_gui.QDesktopServices, 'openUrl',
                        lambda url: files.append(url.toLocalFile()) or True)
    return files


def start_run(window, tmp_path, monkeypatch):
    window.start_pipeline_signal.disconnect()
    window.start_pipeline_signal.connect(lambda settings: None)
    monkeypatch.setattr(pipeline_gui, 'environment_problems', lambda use_docker: [])
    monkeypatch.setattr(window, '_docker_available', lambda: True)
    window.docker_check.setChecked(True)
    window.input_dir_label.setText(str(tmp_path))
    window.output_dir_label.setText(str(tmp_path))
    window.db_label.setText(str(tmp_path / 'db.fasta'))
    window.run_pipeline()
    assert window.is_running


def test_performance_report_opens_when_profiling_is_ticked(window, opened, tmp_path, monkeypatch):
    start_run(window, tmp_path, monkeypatch)
    (tmp_path / pipeline_gui.PERFORMANCE_REPORT_NAME).write_text('<html></html>')
    window.pipeline_finished(True, None)
    assert opened == [str(tmp_path / 'performance_report.html')]


def test_performance_report_is_not_opened_when_profiling_is_unticked(window, opened, tmp_path, monkeypatch):
    window.profile_check.setChecked(False)
    start_run(window, tmp_path, monkeypatch)
    (tmp_path / pipeline_gui.PERFORMANCE_REPORT_NAME).write_text('<html></html>')
    window.pipeline_finished(True, None)
    assert opened == []


def test_performance_report_of_an_earlier_run_is_not_opened(window, opened, tmp_path, monkeypatch):
    (tmp_path / pipeline_gui.PERFORMANCE_REPORT_NAME).write_text('<html>earlier run</html>')
    start_run(window, tmp_path, monkeypatch)
    window.pipeline_finished(True, None)
    assert opened == []
    assert 'Performance report not generated' in window.log_output.toPlainText()


def test_performance_report_is_not_opened_when_the_run_fails(window, opened, dialogs, tmp_path, monkeypatch):
    start_run(window, tmp_path, monkeypatch)
    (tmp_path / pipeline_gui.PERFORMANCE_REPORT_NAME).write_text('<html></html>')
    window.pipeline_finished(False, None)
    assert opened == []


def test_analysis_report_is_opened_after_the_performance_report(window, opened, tmp_path, monkeypatch):
    report = tmp_path / 'Report_run_db.html'
    report.write_text('<html></html>')
    start_run(window, tmp_path, monkeypatch)
    (tmp_path / pipeline_gui.PERFORMANCE_REPORT_NAME).write_text('<html></html>')
    window.pipeline_finished(True, str(report))
    assert opened == [str(tmp_path / 'performance_report.html'), str(report)]


def test_profiler_writes_the_report_that_the_interface_opens(tmp_path):
    from pipeline_profiler import PipelineProfiler
    profiler = PipelineProfiler()
    profiler.start_step('Step 1: Taxonomy')
    profiler.end_step('Step 1: Taxonomy')
    profiler.save(tmp_path)
    assert (tmp_path / pipeline_gui.PERFORMANCE_REPORT_NAME).is_file()
