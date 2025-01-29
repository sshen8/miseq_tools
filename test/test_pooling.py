import pandas as pd
from miseq_tools.pooling import _pools

def test_pooling():
    num_reads = pd.Series({
        "Sample1": 10000000,
        "Sample2": 20000000,
    })
    concs = pd.Series({ # nM
        "Sample1": 10,
        "Sample2": 8.5,
    })
    min_ul_pipettable = 2
    max_ul_pipettable = 10
    pools = _pools(num_reads, concs, min_ul_pipettable=min_ul_pipettable, max_ul_pipettable=max_ul_pipettable)

    # check all pipetting volumes are within range
    for pool_i, pool in enumerate(pools, 1):
        for sample, ul in pool.items():
            assert min_ul_pipettable <= ul <= max_ul_pipettable, f"Pool {pool_i} asks to pipette {ul} uL for {sample}, which is out of range"
