"""Spelling of the messages of the mapping step."""
import shutil
import subprocess

import pytest

from helpers import REPO_ROOT, run_cli


@pytest.mark.real_tools
def test_step_4_is_called_mapping_and_tabulation(tmp_path, dataset):
    returncode, output, hung = run_cli(
        ['-i', dataset['input_dir'], '-o', tmp_path / 'out', '-d', dataset['db_path'],
         '-T', '2', '-t', '1', '--debug'], timeout=120)
    assert returncode == 0, output
    assert 'STEP 4: Mapping and Tabulation' in output
    assert 'Tabeling' not in output


def test_progress_label_is_spelled_tabulation():
    pytest.importorskip('PyQt5')
    import pipeline_worker
    worker = pipeline_worker.PipelineWorker()
    labels = []
    worker.progress_signal.connect(lambda value, text: labels.append(text))
    worker.update_progress_from_log('--- STEP 4: Mapping and Tabulation ---')
    worker.update_progress_from_log('  Step 4 progress: 1/2 samples completed (elapsed 00:01:00)')
    assert labels and all('Tabulation' in label and 'Tabeling' not in label for label in labels)


def test_no_tracked_text_still_shows_the_misspelling():
    if shutil.which('git') is None or not (REPO_ROOT / '.git').exists():
        pytest.skip('not a git checkout')
    result = subprocess.run(['git', 'grep', '-l', 'Tabeling', '--', '.', ':!tests'],
                            cwd=REPO_ROOT, capture_output=True, text=True)
    assert result.returncode in (0, 1)
    assert result.stdout.split() == []
