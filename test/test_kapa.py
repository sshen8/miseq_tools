import subprocess
import pytest
import tempfile
import os

@pytest.mark.parametrize(['kapafolder', 'samplesheet'], [
    ('data/SPS246/kapadata', 'data/SPS246/SPS246 miseq - Sheet1.csv'),
    ('data/SPS303/admin_2024-06-30 21-11-35_BR004200 sps303 ngs quant', 'data/SPS303/SPS303 miseq - Sheet1.csv'),
], ids=['SPS246', 'SPS303'])
def test_kapa(kapafolder, samplesheet):
    curdir = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        out = subprocess.run(['python', '-m', 'miseq_tools', 'kapa', os.path.join(curdir, kapafolder), os.path.join(curdir, samplesheet)], capture_output=True)
        assert out.returncode == 0
        txt = out.stderr.decode()
        assert 'Efficiency:' in txt
        assert 'Slope:' in txt
        assert 'R²:' in txt
        assert 'Intercept:' in txt
    os.chdir(curdir)
