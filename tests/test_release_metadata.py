"""Licence and version number must be stated once and agree everywhere."""
import hashlib
import re
import shutil
import subprocess

import pytest

import OtuUtils
from helpers import REPO_ROOT, run_cli

GPL3_SHA256 = '3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986'


def read(name):
    return (REPO_ROOT / name).read_text(encoding='utf-8')


def test_licence_file_is_the_gnu_gpl_version_3():
    data = (REPO_ROOT / 'LICENSE').read_bytes()
    assert hashlib.sha256(data).hexdigest() == GPL3_SHA256


def test_licence_is_declared_consistently():
    assert re.search(r'^license: GPL-3\.0-only$', read('CITATION.cff'), re.MULTILINE)
    for name in ('README.md', 'credits.md'):
        text = read(name)
        assert 'GNU General Public License' in text, name
        assert 'affero' not in text.lower(), name


def test_version_has_three_numbers():
    assert re.fullmatch(r'\d+\.\d+\.\d+', OtuUtils.LOREON_VERSION)


def test_version_is_the_same_in_every_file_that_states_it():
    version = OtuUtils.LOREON_VERSION
    assert f'version: "{version}"' in read('CITATION.cff')
    assert f'#define AppVersion   "{version}"' in read('installer.iss')
    assert f'APP_VERSION="{version}"' in read('build_deb.sh')
    assert f'APP_VERSION="{version}"' in read('build_dmg.sh')
    readme = read('README.md')
    assert f'(v{version})' in readme
    assert 'version   = {' + version + '}' in readme


@pytest.mark.real_tools
def test_run_shows_the_version(tmp_path, dataset):
    returncode, output, hung = run_cli(
        ['-i', dataset['input_dir'], '-o', tmp_path / 'out', '-d', dataset['db_path'], '-T', '2', '-t', '1'],
        timeout=120)
    assert returncode == 0, output
    assert f"(v{OtuUtils.LOREON_VERSION})" in output


def test_readme_tells_loreon_apart_from_lorean():
    assert 'LoReAn' in read('README.md')


def test_no_installer_binaries_are_versioned():
    if shutil.which('git') is None or not (REPO_ROOT / '.git').exists():
        pytest.skip('not a git checkout')
    tracked = subprocess.run(['git', 'ls-files', '--', '*.exe', '*.dmg', '*.deb', '*.msi'],
                             cwd=REPO_ROOT, capture_output=True, text=True, check=True).stdout.split()
    present = [name for name in tracked if (REPO_ROOT / name).exists()]
    assert present == []


def test_debian_package_states_the_same_licence():
    script = read('build_deb.sh')
    assert 'License: MIT' not in script
    assert 'License: GPL-3' in script
    assert '/usr/share/common-licenses/GPL-3' in script


@pytest.mark.parametrize('script', ['build_exe.bat', 'build_mac.sh', 'build_linux.sh',
                                    'build_deb.sh', 'build_dmg.sh'])
def test_every_package_ships_the_licence_text(script):
    lines = [line for line in read(script).splitlines()
             if 'LICENSE' in line and not line.lstrip().startswith(('#', 'REM', '::'))]
    assert lines, script
