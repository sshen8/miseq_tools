from miseq_tools.samplesheet import check_indexes, format_samplesheet
import pandas as pd
import logging
import pytest
import tempfile

LOGGER = logging.getLogger(__name__)

@pytest.mark.parametrize("i7,ok", [
    ("GCTGAGAA", False),
    ("TTCTCAGC", True),
], ids=["bad", "good"])
def test_revcomp_i7(i7, ok, caplog):
    df = pd.DataFrame.from_records([dict(
        Sample_ID="test",
        I7_Index_ID="BC-A06",
        index=i7,
        I5_Index_ID="P5 mix",
        index2="TCTTTCCC",
    )])
    with caplog.at_level(logging.WARNING):
        check_indexes(df)
    assert ('reverse complemented' not in caplog.text) == ok

def test_swapped(caplog):
    df = pd.DataFrame.from_records([dict(
        Sample_ID="test",
        I7_Index_ID="FWD-2a",
        index="ATTACTCG",
        I5_Index_ID="Ad2_1",
        index2="TAAGGCGA",
    )])
    with caplog.at_level(logging.WARNING):
        check_indexes(df)
    assert 'swapped' in caplog.text

@pytest.mark.parametrize("sample_name", ["SPS375 9", "SPS375.9"])
@pytest.mark.parametrize("fix_invalid_characters", [True, False])
@pytest.mark.parametrize("cli", [True, False])
def test_invalid_characters(caplog, fix_invalid_characters, sample_name, cli):
    with tempfile.TemporaryDirectory() as tmpdir:
        fname_in = f"{tmpdir}/in.csv"
        with open(fname_in, 'wt') as f:
            f.write(f"""
Pool label,Sample_ID,I7_Index_ID,index,I5_Index_ID,index2,Reads (million),Amplicon size (bp),READ STRUCTURE,Read 1,278,
SPS375_14,{sample_name},Ad2_1,TAAGGCGA,Ad1_77,TAGGTATG,0.650,558,300-cycle 5M miseq i100,Index 1 (i7),8,
,SPS375_10,Ad2_2,CGTACTAG,Ad1_77,TAGGTATG,0.650,558,,Index 2 (i5),8,
,SPS375_11,Ad2_3,AGGCAGAA,Ad1_77,TAGGTATG,0.650,558,,Read 2,44,
,SPS375_12,Ad2_4,TCCTGAGC,Ad1_77,TAGGTATG,0.650,558,,TOTAL,338,out of 338
,,,,,,,,,,,
,,,,,,,,,,,
PhiX,2.00%,n/a,n/a,n/a,n/a,0.100,500,,,,
TOTAL,,,,,,5.000,out of 5M,,,,
""")
        fname_out = f"{tmpdir}/out.csv"
        if cli:
            import subprocess
            args = ["python", "-m", "miseq_tools", "sheet", fname_in, "-o", fname_out]
            if not fix_invalid_characters:
                args.append("--dont-fix-invalid-characters")
            out = subprocess.run(args, capture_output=True)
            assert out.returncode == 0
            log = out.stderr.decode()
        else:
            with caplog.at_level(logging.WARNING):
                format_samplesheet(fname_in=fname_in, fname_out=fname_out, fix_invalid_characters=fix_invalid_characters)
            log = caplog.text
        df_out = pd.read_csv(fname_out, skiprows=16)
    assert 'Invalid characters' in log
    if fix_invalid_characters:
        assert df_out.loc[0, 'Sample_ID'] == 'SPS375_9'
    else:
        assert df_out.loc[0, 'Sample_ID'] == sample_name
