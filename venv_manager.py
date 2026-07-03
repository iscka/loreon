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

import os
import subprocess
import sys
from pathlib import Path


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
    """Find a usable system python3 interpreter (≥3.8)."""
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
                    if int(major) >= 3 and int(minor) >= 8:
                        return cmd
        except (FileNotFoundError, subprocess.TimeoutExpired, ValueError):
            continue

    raise RuntimeError(
        "Cannot find Python 3.8+ on this system.\n"
        "Please install Python 3.8 or later and make sure it is in your PATH."
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
