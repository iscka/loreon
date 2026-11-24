import os
import sys
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parent
SOURCES_DIR = REPO_ROOT / 'sources'
for entry in (str(SOURCES_DIR), str(TESTS_DIR)):
    if entry not in sys.path:
        sys.path.insert(0, entry)

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from helpers import HAVE_TOOLS, SYSTEM_PATH, build_dataset


def pytest_collection_modifyitems(config, items):
    skip_tools = pytest.mark.skip(reason='minimap2 and samtools are required')
    skip_posix = pytest.mark.skip(reason='stand-in tools are shell scripts')
    for item in items:
        if 'real_tools' in item.keywords and not HAVE_TOOLS:
            item.add_marker(skip_tools)
        if 'fake_tools' in item.keywords and os.name == 'nt':
            item.add_marker(skip_posix)


def pytest_configure(config):
    config.addinivalue_line('markers', 'real_tools: needs minimap2 and samtools on PATH')
    config.addinivalue_line('markers', 'fake_tools: replaces external tools with shell scripts')


@pytest.fixture
def dataset(tmp_path):
    return build_dataset(tmp_path / 'data')


@pytest.fixture
def fake_bin(tmp_path, monkeypatch):
    """Directory placed first on PATH; the real minimap2 and samtools are hidden."""
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    monkeypatch.setenv('PATH', f"{bin_dir}:{SYSTEM_PATH}")
    return bin_dir
