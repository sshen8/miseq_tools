from miseq_tools.samplesheet import check_indexes
import pandas as pd
import logging
import pytest

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
