import pandas as pd
from .utils import parse_samplesheet, pooled_reads

def pooling(samplesheet: str, quant_csv: str, **kwargs):
    samples = parse_samplesheet(samplesheet)
    num_reads = (pooled_reads(samples) * 1e6).astype(int)
    concs = pd.read_csv(quant_csv, index_col=0)["nM"]
    pools = _pools(num_reads, concs, **kwargs)
    _check_samples_used_exactly_once(pools, set(num_reads.index))
    _check_dilution(pools, num_reads, concs)
    for i, pool in enumerate(pools, 1):
        print(f"Pool {i}")
        # show water first, then descending order of volume
        for sample, ul in sorted(pool.items(), key=lambda item: (item[0] == 'Water', item[1]), reverse=True):
            print(f"[ ] {sample:<10}:\t{ul:.2f} uL")

def _pools(num_reads: dict[str, int],
           concs: dict[str, float],
           min_ul_pipettable: float = 2,
           max_ul_pipettable: float = 10,
           min_ul_total: float = 10,
           ) -> list[dict[str, float]]:
    num_reads = pd.Series(num_reads)
    concs = pd.Series(concs)
    frac_reads = num_reads / num_reads.sum()
    desired_nm = frac_reads * 4
    ul = min_ul_total * desired_nm / concs
    assert all(not sample_name.startswith('Water') for sample_name in ul.index)
    assert all(not sample_name == 'Prev Pool' for sample_name in ul.index)
    ul['Water'] = min_ul_total - ul.sum()
    assert ul['Water'] >= 0, "Some sample(s) is/are not concentrated enough."

    # initial assign dilution factor groups
    def _assign_dilution_factors(_ul):
        _dilution_factors = pd.Series(None, index=_ul.index, dtype=float)
        while any(_dilution_factors.isna()):
            remaining_samples = _dilution_factors[_dilution_factors.isna()].index
            dilution_factor = min_ul_pipettable / min(_ul[remaining_samples])
            if dilution_factor <= 1:
                _dilution_factors[remaining_samples] = 1
                break
            hypothetical_volumes = _ul[remaining_samples] * dilution_factor
            samples_in_this_pool = hypothetical_volumes.index[hypothetical_volumes <= max_ul_pipettable]
            _dilution_factors[samples_in_this_pool] = dilution_factor
        return _dilution_factors
    
    # check previous pool volumes and dilute with water if necessary
    def _split_water(_ul, _dilution_factors):
        _ul_new = _ul.copy()
        dilution_factors_list = sorted(_dilution_factors.unique(), reverse=True)
        for i, dilution_factor in enumerate(dilution_factors_list, 1):
            samples_in_previous_pool = _dilution_factors[_dilution_factors > dilution_factor].index
            if len(samples_in_previous_pool) == 0:
                continue
            ul_prev_pool = _ul[samples_in_previous_pool].sum() * dilution_factor
            if ul_prev_pool >= min_ul_pipettable:
                continue
            ul_water_split = (max_ul_pipettable - ul_prev_pool) / dilution_factor
            if i >= 2:
                ul_water_split = min(ul_water_split, max_ul_pipettable / dilution_factors_list[i - 2])
            _ul_new['Water'] -= ul_water_split
            _ul_new[f'Water {i}'] = _ul_new.get(f'Water {i}', 0) + ul_water_split
            break
        return _ul_new

    ul_old = None
    ul_new = ul.copy()
    while ul_old is None or not ul_old.equals(ul_new):
        ul_old = ul_new.copy()
        dilution_factors = _assign_dilution_factors(ul_new)
        ul_new = _split_water(ul_new, dilution_factors)

    # calculate volumes per pool
    pools = []
    dilution_factors_list = sorted(dilution_factors.unique(), reverse=True)
    for dilution_factor in dilution_factors_list:
        samples_in_previous_pool = dilution_factors[dilution_factors > dilution_factor].index
        if len(samples_in_previous_pool) == 0:
            ul_prev_pool = dict()
        else:
            ul_prev_pool = {'Prev Pool': ul_new[samples_in_previous_pool].sum() * dilution_factor}
        samples_in_this_pool = dilution_factors[dilution_factors == dilution_factor].index
        ul_diluted = ul_new[samples_in_this_pool] * dilution_factor
        ul_diluted.index = ul_diluted.index.map(lambda x: 'Water' if x.startswith('Water') else x)
        pools.append({
            **ul_diluted.to_dict(),
            **ul_prev_pool,
        })

    return pools

def _check_dilution(pools: list[dict[str, float]], num_reads: pd.Series, concs: pd.Series):
    volume_so_far = pd.Series(0, index=num_reads.index, dtype=float)
    for pool in pools:
        volume_pool = pd.Series(pool, index=num_reads.index.to_list() + ['Water', 'Prev Pool']).fillna(0)
        ul_prev_pool = volume_pool.pop('Prev Pool')
        if ul_prev_pool > 0:
            volume_pool += volume_so_far * ul_prev_pool / volume_so_far.sum()
        volume_so_far = volume_pool
    concs_final = concs * volume_so_far / volume_so_far.sum()
    concs_final.drop('Water', inplace=True)
    fracs_final = concs_final / 4
    pd.testing.assert_series_equal(fracs_final[num_reads.index], num_reads / num_reads.sum(), check_names=False)

def _check_samples_used_exactly_once(pools: list[dict[str, float]], all_samples: set[str]):
    unused_samples = all_samples.copy()
    for pool in pools:
        samples_in_this_pool = set(pool.keys()) - {'Water', 'Prev Pool'}
        assert samples_in_this_pool.issubset(all_samples), f"Some sample(s) in pool {pool} are not valid samples"
        assert samples_in_this_pool.issubset(unused_samples), f"Some sample(s) in pool {pool} have already been added to a previous pool"
        unused_samples -= samples_in_this_pool
    assert len(unused_samples) == 0, f"Samples {','.join(unused_samples)} have not been added to any pool"
