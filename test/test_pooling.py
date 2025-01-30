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

    # check all samples have been added exactly once
    all_samples = set(num_reads.index)
    _check_samples_used_exactly_once(pools, all_samples)

def _check_samples_used_exactly_once(pools: list[dict[str, float]], all_samples: set[str]):
    unused_samples = all_samples.copy()
    for pool in pools:
        samples_in_this_pool = set(pool.keys()) - {'Water', 'Prev Pool'}
        assert samples_in_this_pool.issubset(all_samples), f"Some sample(s) in pool {pool} are not valid samples"
        assert samples_in_this_pool.issubset(unused_samples), f"Some sample(s) in pool {pool} have already been added to a previous pool"
        unused_samples -= samples_in_this_pool
    assert len(unused_samples) == 0, f"Samples {','.join(unused_samples)} have not been added to any pool"

@pytest.mark.parametrize("pools,expected", [
    # unused
    ([], False),
    ([{'Sample1': 1}], False),
    ([{'Sample2': 1}], False),
    # extra
    ([{'Sample1': 1, 'Sample2': 1, 'Sample3': 1}], False),
    # multiple
    ([{'Sample1': 1, 'Sample2': 1}, {'Sample1': 1}], False),
    # good
    ([{'Sample1': 1, 'Sample2': 1}], True),
    ([{'Sample1': 1}, {'Sample2': 1}], True),
])
def test_check_samples_used_exactly_once(pools, expected):
    all_samples = {'Sample1', 'Sample2'}
    if expected:
        _check_samples_used_exactly_once(pools, all_samples)
    else:
        with pytest.raises(AssertionError):
            _check_samples_used_exactly_once(pools, all_samples)
