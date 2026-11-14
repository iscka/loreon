#!/usr/bin/env python3
"""LOREON Virtual Environment Manager.

When the application is distributed via PyInstaller (frozen), the bundled
Python cannot run the pipeline scripts (metaGenomics_new.py, report_generator.py)
because those scripts need full packages (pandas, duckdb, plotly, etc.) that
are NOT bundled by PyInstaller (only the GUI dependencies are bundled).

This module:
  1. Detects whether a suitable virtualenv already exists.
  2. Creates one (using the system ``python3``) on first launch if needed.
  3. Installs ``requirements.txt`` into the venv via pip.
  4. Exposes ``get_venv_python()`` so that ``pipeline_worker.py`` can use
     the venv interpreter to run pipeline scripts.

The venv is stored in a platform-appropriate user data directory so it
persists across application updates.
"""

import importlib.util
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

PIPELINE_MODULES = ('pandas', 'numpy', 'duckdb', 'openpyxl', 'plotly', 'jinja2')
CONDA_ENVIRONMENT_HINT = "conda env create -f environment.yml && conda activate loreon"
_IRREGULAR_PACKAGES = {
    'libglib-2.0.so.0': 'libglib2.0-0',
    'libgthread-2.0.so.0': 'libglib2.0-0',
    'libgobject-2.0.so.0': 'libglib2.0-0',
    'libgio-2.0.so.0': 'libglib2.0-0',
    'libz.so.1': 'zlib1g',
    'libharfbuzz.so.0': 'libharfbuzz0b',
}
_MISSING_LIBRARY = re.compile(r'^\s*(\S+\.so[.\d]*)\s+=>\s+not found', re.MULTILINE)
_LIBRARY_IN_IMPORT_ERROR = re.compile(r'(lib[\w.+-]+\.so[.\d]*): cannot open shared object file')


def debian_package_for(library):
    """Name of the Debian/Ubuntu package that usually provides a shared library.

    The distribution rule is applied: 'libfoo.so.1' is packaged as 'libfoo1',
    and 'libfoo-2.so.3' as 'libfoo-2-3'.
    """
    if library in _IRREGULAR_PACKAGES:
        return _IRREGULAR_PACKAGES[library]
    match = re.match(r'(.+?)\.so\.?([\d.]*)$', library)
    if not match:
        return library.lower()
    name, version = match.group(1).lower(), match.group(2).split('.')[0]
    if not version:
        return name
    return f"{name}-{version}" if name[-1].isdigit() else f"{name}{version}"


def missing_libraries_from_ldd(output):
    """Shared libraries reported as 'not found' in the output of ldd, sorted, without repeats."""
    return sorted(set(_MISSING_LIBRARY.findall(output)))


def _missing_library_message(libraries):
    packages = ' '.join(sorted({debian_package_for(library) for library in libraries}))
    return "\n".join([
        "LOREON cannot open its window: shared libraries needed by the Qt toolkit are missing",
        "on this system:",
        f"  {', '.join(libraries)}",
        "On Ubuntu and Debian they are usually provided by these packages:",
        f"  sudo apt install {packages}",
        "The conda environment described in the README brings its own copies:",
        f"  {CONDA_ENVIRONMENT_HINT}",
    ])


def _pyqt_binaries():
    """Qt libraries and display plugin inside a PyQt5 installed with pip, when present."""
    try:
        spec = importlib.util.find_spec('PyQt5')
    except (ImportError, ValueError):
        return []
    if spec is None or not spec.submodule_search_locations:
        return []
    root = Path(list(spec.submodule_search_locations)[0]) / 'Qt5'
    candidates = [root / 'plugins' / 'platforms' / 'libqxcb.so']
    candidates += [root / 'lib' / f'libQt5{name}.so.5'
                   for name in ('Core', 'Gui', 'Widgets', 'XcbQpa', 'DBus')]
    return [path for path in candidates if path.exists()]


def pyqt_import_problem(error, binaries=None, linker=None):
    """Message to show when PyQt5 cannot be imported.

    When the cause is a missing system library, the Qt libraries of PyQt5
    are examined so that every missing library is listed at once.
    """
    library = _LIBRARY_IN_IMPORT_ERROR.search(str(error))
    if library:
        linker = _run_ldd if linker is None else linker
        missing = {library.group(1)}
        for binary in (_pyqt_binaries() if binaries is None else binaries):
            missing.update(missing_libraries_from_ldd(linker(binary)))
        return _missing_library_message(sorted(missing))
    return "\n".join([
        f"LOREON cannot open its window: the PyQt5 package cannot be imported ({error}).",
        "Install the Python packages of LOREON in the environment you are using:",
        "  pip install -r requirements.txt",
        "or create the conda environment described in the README:",
        f"  {CONDA_ENVIRONMENT_HINT}",
    ])


def _run_ldd(path):
    try:
        result = subprocess.run(['ldd', str(path)], capture_output=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return ''
    return result.stdout.decode('utf-8', errors='replace')


def qt_library_problem(plugins_dir, linker=_run_ldd):
    """Message listing the system libraries that the Qt display plugin needs and cannot find.

    plugins_dir is the Qt plugin folder; the 'xcb' platform plugin inside it
    is examined with ldd. Returns None when nothing is missing or the plugin
    is not there.
    """
    plugin = Path(plugins_dir) / 'platforms' / 'libqxcb.so'
    if not plugin.exists():
        return None
    libraries = missing_libraries_from_ldd(linker(plugin))
    return _missing_library_message(libraries) if libraries else None


def display_problem(platform=None, environ=None):
    """Message to show on Linux when no graphical session is available, else None."""
    platform = sys.platform if platform is None else platform
    environ = os.environ if environ is None else environ
    if not platform.startswith('linux'):
        return None
    if environ.get('DISPLAY') or environ.get('WAYLAND_DISPLAY') or environ.get('QT_QPA_PLATFORM'):
        return None
    return "\n".join([
        "LOREON cannot open its window: no graphical session was found",
        "(DISPLAY and WAYLAND_DISPLAY are not set).",
        "On a remote or headless computer run the analysis from the command line:",
        "  python metaGenomics_new.py -i <input folder> -o <output folder> -d <database.fasta>",
        "or connect with display forwarding, for example 'ssh -X'.",
    ])


def environment_problems(use_docker, platform=None, which=shutil.which,
                         find_spec=importlib.util.find_spec):
    """Components that are needed to run an analysis and are missing, one message each.

    In Docker mode, and always on Windows, only Docker is needed on this
    computer. Otherwise minimap2, samtools and the Python packages of the
    pipeline must be installed where the interface runs.
    """
    platform = sys.platform if platform is None else platform
    if use_docker or platform == 'win32':
        if which('docker') is None:
            return ["Docker was not found: install Docker Desktop "
                    "(https://www.docker.com/products/docker-desktop) and start it."]
        return []
    problems = []
    for tool in ('minimap2', 'samtools'):
        if which(tool) is None:
            problems.append(
                f"{tool} was not found: install it with "
                f"'conda install -c conda-forge -c bioconda minimap2=2.28 samtools=1.21', "
                f"create the conda environment of the README, or tick 'Use Docker'."
            )
    if not getattr(sys, 'frozen', False):
        missing = [name for name in PIPELINE_MODULES if find_spec(name) is None]
        if missing:
            problems.append(
                f"Python packages not found: {', '.join(missing)}: install them with "
                f"'pip install -r requirements.txt' or create the conda environment of the README."
            )
    return problems


def _data_dir() -> Path:
    """Return platform-specific application data directory."""
    if sys.platform == 'darwin':
        return Path.home() / 'Library' / 'Application Support' / 'LOREON'
    elif sys.platform == 'win32':
        return Path(os.environ.get('LOCALAPPDATA', Path.home())) / 'LOREON'
    else:
        return Path(os.environ.get('XDG_DATA_HOME', Path.home() / '.local' / 'share')) / 'loreon'


VENV_DIR = _data_dir() / 'venv'


def _base_dir() -> Path:
    """Return the directory containing the application files."""
    if getattr(sys, 'frozen', False):
        return Path(sys.executable).parent
    return Path(__file__).parent


def _requirements_path() -> Path:
    return _base_dir() / 'requirements.txt'


def _find_system_python() -> str:
    """Find a usable system python3 interpreter (≥3.9)."""
    candidates = ['python3', 'python']
    if sys.platform == 'win32':
        candidates = ['python', 'python3', 'py -3']

    for cmd in candidates:
        try:
            result = subprocess.run(
                cmd.split() + ['--version'],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0:
                version_str = result.stdout.strip() or result.stderr.strip()
                parts = version_str.split()
                if len(parts) >= 2:
                    major, minor = parts[1].split('.')[:2]
                    if (int(major), int(minor)) >= (3, 9):
                        return cmd
        except (FileNotFoundError, subprocess.TimeoutExpired, ValueError):
            continue

    raise RuntimeError(
        "Cannot find Python 3.9+ on this system.\n"
        "Please install Python 3.9 or later and make sure it is in your PATH."
    )


def _venv_python() -> Path:
    """Return the path to the venv's python interpreter."""
    if sys.platform == 'win32':
        return VENV_DIR / 'Scripts' / 'python.exe'
    return VENV_DIR / 'bin' / 'python3'


def is_venv_ready() -> bool:
    """Check whether the venv exists and has a working interpreter."""
    py = _venv_python()
    if not py.exists():
        return False
    try:
        result = subprocess.run(
            [str(py), '-c', 'import pandas, duckdb, plotly; print("ok")'],
            capture_output=True, text=True, timeout=30
        )
        return result.returncode == 0 and 'ok' in result.stdout
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def create_venv(progress_callback=None):
    """Create the virtualenv and install requirements.

    Parameters
    ----------
    progress_callback : callable, optional
        Called with (message: str) for each step to provide UI feedback.

    Raises
    ------
    RuntimeError
        If any step fails.
    """
    def _log(msg):
        if progress_callback:
            progress_callback(msg)

    _log("Looking for system Python...")
    system_python = _find_system_python()
    _log(f"Found: {system_python}")

    VENV_DIR.parent.mkdir(parents=True, exist_ok=True)
    _log(f"Creating virtual environment in {VENV_DIR}...")
    result = subprocess.run(
        system_python.split() + ['-m', 'venv', str(VENV_DIR), '--clear'],
        capture_output=True, text=True, timeout=120
    )
    if result.returncode != 0:
        raise RuntimeError(f"Failed to create venv:\n{result.stderr}")

    py = str(_venv_python())
    if not Path(py).exists():
        raise RuntimeError(f"Venv python not found at {py} after creation.")

    _log("Upgrading pip...")
    result = subprocess.run(
        [py, '-m', 'pip', 'install', '--upgrade', 'pip'],
        capture_output=True, text=True, timeout=120
    )
    if result.returncode != 0:
        _log(f"WARNING: pip upgrade failed (continuing): {result.stderr}")

    req_path = _requirements_path()
    if not req_path.exists():
        raise RuntimeError(f"requirements.txt not found at {req_path}")

    _log(f"Installing dependencies from {req_path.name}...")
    result = subprocess.run(
        [py, '-m', 'pip', 'install', '-r', str(req_path)],
        capture_output=True, text=True, timeout=600
    )
    if result.returncode != 0:
        raise RuntimeError(f"pip install failed:\n{result.stderr}")

    _log("Virtual environment ready.")


def get_venv_python() -> str:
    """Return the path to the venv python interpreter.

    When running from source (not frozen), returns ``sys.executable``
    (the current interpreter) — no venv needed.

    When frozen (PyInstaller), returns the venv interpreter path.
    Raises RuntimeError if the venv is not ready.
    """
    if not getattr(sys, 'frozen', False):
        return sys.executable

    py = _venv_python()
    if not py.exists():
        raise RuntimeError(
            "LOREON virtual environment not found.\n"
            "Please restart the application to trigger automatic setup."
        )
    return str(py)
