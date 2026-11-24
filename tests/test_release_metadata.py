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


def test_windows_installer_does_not_install_python_packages():
    script = read('installer.iss')
    assert 'pip install' not in script
    assert 'Python' not in script.split('[Code]')[1]


def test_macos_bundle_is_signed_after_its_last_change():
    lines = read('build_dmg.sh').splitlines()
    changes = [number for number, line in enumerate(lines)
               if 'PlistBuddy' in line or '"$BUNDLE_DIR/"' in line]
    signing = [number for number, line in enumerate(lines) if line.lstrip().startswith('codesign --force')]
    staging = [number for number, line in enumerate(lines) if '"$DMG_STAGE/"' in line]
    assert changes and signing and staging
    assert max(changes) < signing[0] < min(staging)


def test_debian_package_installs_executable_commands():
    lines = read('build_deb.sh').splitlines()
    reset = max(number for number, line in enumerate(lines) if 'chmod 644' in line)
    for command in ('loreon', 'loreon-cli'):
        granted = [number for number, line in enumerate(lines)
                   if line.strip() == f'chmod 755 "$DEB_ROOT/usr/local/bin/{command}"']
        assert granted and max(granted) > reset, command


def test_debian_command_line_uses_the_environment_of_the_application():
    script = read('build_deb.sh')
    assert 'ln -sf /opt/loreon/metaGenomics_new.py' not in script
    assert "<< 'CLIWRAPPER'" in script
    wrapper = script.split("<< 'CLIWRAPPER'")[1].split('\nCLIWRAPPER')[0]
    assert '/loreon/venv/bin/python3' in wrapper
    assert '/opt/loreon/metaGenomics_new.py "$@"' in wrapper


def test_download_page_links_the_installers_of_this_version():
    version = OtuUtils.LOREON_VERSION
    page = read('installers/README.md')
    for name in (f'LOREON_Setup_{version}.exe', f'LOREON-{version}-macOS.dmg', f'loreon_{version}_amd64.deb'):
        assert f'/releases/download/v{version}/{name}' in page, name


def test_checksums_are_those_of_installers_listed_in_the_download_page():
    version = OtuUtils.LOREON_VERSION
    page = read('installers/README.md')
    entries = [line.split() for line in read('installers/SHA256SUMS').splitlines() if line.strip()]
    assert entries
    for digest, name in entries:
        assert re.fullmatch(r'[0-9a-f]{64}', digest), name
        assert version in name and name in page, name
