import pandas as pd
from miseq_tools.demux_stats import demux
from miseq_tools.utils import parse_samplesheet, pooled_reads

SAMPLESHEET = "test/data/SPS412 miseq - Sheet1.csv"
STATS = "test/data/Demultiplex_Stats.csv"

def test_demux_multilane(tmp_path):
    df_pool = demux(SAMPLESHEET, STATS, out=tmp_path / "demux_stats.pdf")

    # intended reads should not be multiplied by the number of lanes
    expected_intended = pooled_reads(parse_samplesheet(SAMPLESHEET)) * 1e6
    pd.testing.assert_series_equal(df_pool["intended"], expected_intended.reindex(df_pool.index), check_names=False)

    # actual reads should be summed across lanes
    df_stats = pd.read_csv(STATS)
    assert df_stats["Lane"].nunique() > 1
    total = df_stats.loc[df_stats["SampleID"] != "Undetermined", "# Reads"].sum()
    assert df_pool["actual"].sum() == total
