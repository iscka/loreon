"""Checks made when the interface starts, and the messages shown when something is missing."""
import os
import subprocess
import sys

import pytest

import venv_manager
from helpers import SOURCES_DIR

LDD_OUTPUT = """\
\tlinux-vdso.so.1 (0x0000ffff9a2c4000)
\tlibxcb-xinerama.so.0 => not found
\tlibQt5XcbQpa.so.5 => /work/.venv/lib/python3.10/site-packages/PyQt5/Qt5/lib/libQt5XcbQpa.so.5 (0x0000ffff99e00000)
\tlibxcb-icccm.so.4 => not found
\tlibxcb-xinerama.so.0 => not found
\tlibc.so.6 => /lib/aarch64-linux-gnu/libc.so.6 (0x0000ffff99c50000)
"""


def run_gui_with_broken_pyqt(tmp_path, error_statement):
    package = tmp_path / 'broken' / 'PyQt5'
    package.mkdir(parents=True)
    (package / '__init__.py').write_text(error_statement + '\n')
    env = dict(os.environ, PYTHONPATH=f"{tmp_path / 'broken'}{os.pathsep}{SOURCES_DIR}")
    return subprocess.run([sys.executable, str(SOURCES_DIR / 'pipeline_gui.py')],
                          capture_output=True, text=True, env=env, timeout=60)


def test_missing_pyqt_is_reported_with_install_instructions(tmp_path):
    result = run_gui_with_broken_pyqt(
        tmp_path, "raise ModuleNotFoundError(\"No module named 'PyQt5'\", name='PyQt5')")
    assert result.returncode == 1
    assert 'Traceback' not in result.stderr
    assert 'environment.yml' in result.stderr
    assert 'pip install -r requirements.txt' in result.stderr


def test_missing_system_library_is_reported_with_its_package(tmp_path):
    result = run_gui_with_broken_pyqt(
        tmp_path,
        "raise ImportError('libGL.so.1: cannot open shared object file: No such file or directory')")
    assert result.returncode == 1
    assert 'Traceback' not in result.stderr
    assert 'libGL.so.1' in result.stderr
    assert 'sudo apt install libgl1' in result.stderr


def test_missing_libraries_are_read_from_the_linker_output():
    assert venv_manager.missing_libraries_from_ldd(LDD_OUTPUT) == \
        ['libxcb-icccm.so.4', 'libxcb-xinerama.so.0']
    assert venv_manager.missing_libraries_from_ldd('\tlibc.so.6 => /lib/libc.so.6 (0x1)\n') == []


@pytest.mark.parametrize('library,package', [
    ('libxcb-xinerama.so.0', 'libxcb-xinerama0'),
    ('libxcb-icccm.so.4', 'libxcb-icccm4'),
    ('libxkbcommon-x11.so.0', 'libxkbcommon-x11-0'),
    ('libGL.so.1', 'libgl1'),
    ('libdbus-1.so.3', 'libdbus-1-3'),
    ('libX11-xcb.so.1', 'libx11-xcb1'),
    ('libfontconfig.so.1', 'libfontconfig1'),
    ('libX11.so.6', 'libx11-6'),
    ('libglib-2.0.so.0', 'libglib2.0-0'),
    ('libgthread-2.0.so.0', 'libglib2.0-0'),
])
def test_debian_package_name_is_derived_from_the_library_name(library, package):
    assert venv_manager.debian_package_for(library) == package


def test_qt_library_problem_names_libraries_and_packages(tmp_path):
    plugin = tmp_path / 'platforms' / 'libqxcb.so'
    plugin.parent.mkdir()
    plugin.write_text('')
    message = venv_manager.qt_library_problem(tmp_path, linker=lambda path: LDD_OUTPUT)
    assert 'libxcb-icccm.so.4' in message and 'libxcb-xinerama.so.0' in message
    assert 'sudo apt install libxcb-icccm4 libxcb-xinerama0' in message
    assert 'environment.yml' in message


def test_qt_library_problem_is_silent_when_nothing_is_missing(tmp_path):
    plugin = tmp_path / 'platforms' / 'libqxcb.so'
    plugin.parent.mkdir()
    plugin.write_text('')
    assert venv_manager.qt_library_problem(tmp_path, linker=lambda path: '\tlibc.so.6 => /lib/libc.so.6\n') is None
    assert venv_manager.qt_library_problem(tmp_path / 'absent', linker=lambda path: LDD_OUTPUT) is None


def test_missing_display_is_explained_on_linux_only():
    message = venv_manager.display_problem('linux', {})
    assert 'no graphical session' in message
    assert 'metaGenomics_new.py' in message
    assert venv_manager.display_problem('linux', {'DISPLAY': ':0'}) is None
    assert venv_manager.display_problem('linux', {'WAYLAND_DISPLAY': 'wayland-0'}) is None
    assert venv_manager.display_problem('linux', {'QT_QPA_PLATFORM': 'offscreen'}) is None
    assert venv_manager.display_problem('darwin', {}) is None
    assert venv_manager.display_problem('win32', {}) is None


def test_environment_report_names_missing_tools_and_modules():
    problems = venv_manager.environment_problems(
        use_docker=False, platform='linux',
        which=lambda name: None if name == 'samtools' else f'/usr/bin/{name}',
        find_spec=lambda name: None if name == 'duckdb' else object())
    text = '\n'.join(problems)
    assert 'samtools' in text and 'minimap2 was not found' not in text
    assert 'conda' in text
    assert 'duckdb' in text and 'pip install -r requirements.txt' in text


def test_environment_report_is_empty_when_everything_is_present():
    assert venv_manager.environment_problems(
        use_docker=False, platform='linux',
        which=lambda name: f'/usr/bin/{name}', find_spec=lambda name: object()) == []


def test_docker_mode_needs_docker_only():
    problems = venv_manager.environment_problems(
        use_docker=True, platform='win32', which=lambda name: None, find_spec=lambda name: None)
    assert len(problems) == 1 and 'Docker' in problems[0]
    assert venv_manager.environment_problems(
        use_docker=True, platform='win32',
        which=lambda name: 'C:/docker.exe' if name == 'docker' else None,
        find_spec=lambda name: None) == []


def test_packages_are_listed_once_each(tmp_path):
    plugin = tmp_path / 'platforms' / 'libqxcb.so'
    plugin.parent.mkdir()
    plugin.write_text('')
    output = "\tlibglib-2.0.so.0 => not found\n\tlibgthread-2.0.so.0 => not found\n"
    message = venv_manager.qt_library_problem(tmp_path, linker=lambda path: output)
    assert 'sudo apt install libglib2.0-0\n' in message


def test_import_problem_lists_every_missing_library_at_once(tmp_path):
    binary = tmp_path / 'libQt5Gui.so.5'
    binary.write_text('')
    error = ImportError('libGL.so.1: cannot open shared object file: No such file or directory')
    message = venv_manager.pyqt_import_problem(error, binaries=[binary], linker=lambda path: LDD_OUTPUT)
    for library in ('libGL.so.1', 'libxcb-icccm.so.4', 'libxcb-xinerama.so.0'):
        assert library in message
    assert 'sudo apt install libgl1 libxcb-icccm4 libxcb-xinerama0' in message


def test_optional_spreadsheet_writer_is_not_required():
    assert venv_manager.environment_problems(
        use_docker=False, platform='linux', which=lambda name: f'/usr/bin/{name}',
        find_spec=lambda name: None if name == 'xlsxwriter' else object()) == []


def test_packaged_macos_application_gets_the_tool_folders_of_a_terminal():
    path = venv_manager.path_with_tool_dirs('/usr/bin:/bin:/usr/sbin:/sbin', platform='darwin',
                                            isdir=lambda folder: True)
    assert path == '/usr/bin:/bin:/usr/sbin:/sbin:/opt/homebrew/bin:/usr/local/bin'


def test_tool_folders_are_added_once_and_only_when_they_exist():
    assert venv_manager.path_with_tool_dirs(
        '/usr/local/bin:/usr/bin', platform='darwin',
        isdir=lambda folder: True) == '/usr/local/bin:/usr/bin:/opt/homebrew/bin'
    assert venv_manager.path_with_tool_dirs(
        '/usr/bin:/bin', platform='darwin',
        isdir=lambda folder: folder == '/usr/local/bin') == '/usr/bin:/bin:/usr/local/bin'
    assert venv_manager.path_with_tool_dirs(
        '', platform='darwin', isdir=lambda folder: True) == '/opt/homebrew/bin:/usr/local/bin'


def test_path_is_not_changed_on_linux_and_windows():
    for platform in ('linux', 'win32'):
        assert venv_manager.path_with_tool_dirs(
            '/usr/bin', platform=platform, isdir=lambda folder: True) == '/usr/bin'
