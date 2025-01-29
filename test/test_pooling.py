import pandas as pd
from miseq_tools.pooling import _pools
import pytest

@pytest.mark.parametrize("min_ul_pipettable", [1, 2])
@pytest.mark.parametrize("max_ul_pipettable", [10, 3])
def test_pooling(min_ul_pipettable, max_ul_pipettable):
    num_reads = pd.Series({
        "Sample1": 10000000,
        "Sample2": 20000000,
    })
    concs = pd.Series({ # nM
        "Sample1": 10,
        "Sample2": 8.5,
    })
    pools = _pools(num_reads, concs, min_ul_pipettable=min_ul_pipettable, max_ul_pipettable=max_ul_pipettable)

    for pool_i, pool in enumerate(pools, 1):
        for sample, ul in pool.items():
            # check all pipetting volumes are within range
            assert min_ul_pipettable <= ul, f"Pool {pool_i} asks to pipette {ul} uL for {sample}, which is less than the minimum of {min_ul_pipettable} uL"
            if pool_i < len(pools): # last pool
                assert ul <= max_ul_pipettable, f"Pool {pool_i} asks to pipette {ul} uL for {sample}, which is more than the maximum of {max_ul_pipettable} uL"

            # check water is just called "Water"
            if sample.startswith("Water"):
                assert sample == "Water"
