import pandas as pd
from miseq_tools.pooling import _pools, _check_samples_used_exactly_once, _check_dilution
import pytest

@pytest.mark.parametrize("min_ul_pipettable", [1, 2])
@pytest.mark.parametrize("max_ul_pipettable", [10, 5])
@pytest.mark.parametrize("num_reads,concs", [
    ({
        "Sample1": 10000000,
        "Sample2": 20000000,
    }, { # nM
        "Sample1": 10,
        "Sample2": 8.5,
    }),
    ({
        "SPS303A": 2400000,
        "SPS303B": 480000,
        "II": 910000,
    }, { # nM
        "SPS303A": 21.178463,
        "SPS303B": 14.756816,
        "II": 54.082907,
    })
])
def test_pooling(num_reads, concs, min_ul_pipettable, max_ul_pipettable):
    if min_ul_pipettable == 2 and max_ul_pipettable == 5 and set(num_reads.keys()) == {'SPS303A', 'SPS303B', 'II'}:
        pytest.skip("Impossible to solve")

    num_reads = pd.Series(num_reads)
    concs = pd.Series(concs)
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

    # check dilution is correct
    _check_dilution(pools, num_reads, concs)

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

@pytest.mark.parametrize("kwargs,expected", [
    (dict(
        pools=[{'Sample1': 1, 'Sample2': 1}],
        num_reads=pd.Series({'Sample1': 1, 'Sample2': 1}),
        concs=pd.Series({'Sample1': 4, 'Sample2': 4}),
    ), True),
    (dict(
        pools=[{'Sample1': 2, 'Sample2': 1}],
        num_reads=pd.Series({'Sample1': 2, 'Sample2': 1}),
        concs=pd.Series({'Sample1': 4, 'Sample2': 4}),
    ), True),
    (dict(
        pools=[{'Sample1': 1, 'Sample2': 1}],
        num_reads=pd.Series({'Sample1': 2, 'Sample2': 1}),
        concs=pd.Series({'Sample1': 4, 'Sample2': 4}),
    ), False),
])
def test_check_dilutions(kwargs, expected):
    if expected:
        _check_dilution(**kwargs)
    else:
        with pytest.raises(AssertionError):
            _check_dilution(**kwargs)
