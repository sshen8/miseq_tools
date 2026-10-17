import pandas as pd
from .utils import parse_samplesheet, pooled_reads

def pooling(samplesheet: str, quant_csv: str, desired_nm: float = 4, **kwargs):
    samples = parse_samplesheet(samplesheet)
    num_reads = (pooled_reads(samples) * 1e6).astype(int)
    concs = pd.read_csv(quant_csv, index_col=0)["nM"]
    pools = _pools(num_reads, concs, desired_nm=desired_nm, **kwargs)
    _check_samples_used_exactly_once(pools, set(num_reads.index))
    _check_dilution(pools, num_reads, concs, desired_nm=desired_nm)
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
           desired_nm: float = 4,
           ) -> list[dict[str, float]]:
    num_reads = pd.Series(num_reads)
    concs = pd.Series(concs)
    frac_reads = num_reads / num_reads.sum()
    desired_nm_per_pool = frac_reads * desired_nm
    ul = min_ul_total * desired_nm_per_pool / concs
    assert all(not sample_name.startswith('Water') for sample_name in ul.index)
    assert all(not sample_name == 'Prev Pool' for sample_name in ul.index)
    ul['Water'] = min_ul_total - ul.sum()
    assert ul['Water'] >= 0, "Some sample(s) is/are not concentrated enough."

    # Build a chain of pools, most dilute first. Volumes in `ul` are in units of the final pool;
    # pool i pipettes ul * dilution_factor of each of its components. Each pool is transferred
    # into the next, and pure dilution pools (prev pool + water only) are inserted when needed.
    remaining = ul.drop('Water').sort_values()
    pools = []
    dilution_factor = None
    ul_so_far = 0.  # cumulative volume of everything in previous pools, in final pool units
    while dilution_factor is None or dilution_factor > 1:
        u = remaining.iloc[0] if len(remaining) else None
        if dilution_factor is not None:
            # largest dilution factor allowed for the next pool so that its smallest sample is pipettable
            max_next_dilution_factor = max(max_ul_pipettable / u, 1) if u is not None else 1
            # add water to the current pool so that it can be transferred into the next pool
            ul_needed = min_ul_pipettable / max_next_dilution_factor - ul_so_far
            if ul_needed > 0:
                ul_water = min(max(ul_needed, min_ul_pipettable / dilution_factor), max_ul_pipettable / dilution_factor)
                pools[-1]['Water'] = ul_water * dilution_factor
                ul_so_far += ul_water
                ul['Water'] -= ul_water
                if ul_water < ul_needed:
                    # still not enough volume, so dilute the current pool further before adding more samples
                    dilution_factor = min_ul_pipettable / ul_so_far
                    pools.append({'Prev Pool': min_ul_pipettable})
                    continue
        next_dilution_factor = max(min_ul_pipettable / u if u is not None else 1, 1)
        if dilution_factor is not None:
            next_dilution_factor = max(next_dilution_factor, min_ul_pipettable / ul_so_far)
        if next_dilution_factor <= 1:
            next_dilution_factor = 1
            samples_in_this_pool = remaining.index
        else:
            samples_in_this_pool = remaining.index[remaining * next_dilution_factor <= max_ul_pipettable]
        pool = (ul[samples_in_this_pool] * next_dilution_factor).to_dict()
        if dilution_factor is not None:
            ul_prev_pool = ul_so_far * next_dilution_factor
            assert next_dilution_factor == 1 or ul_prev_pool <= max_ul_pipettable, "Unable to find a pooling scheme."
            pool['Prev Pool'] = ul_prev_pool
        pools.append(pool)
        ul_so_far += ul[samples_in_this_pool].sum()
        remaining = remaining.drop(samples_in_this_pool)
        dilution_factor = next_dilution_factor

    assert ul['Water'] >= 0, "Some sample(s) is/are not concentrated enough."
    if ul['Water'] > 0:
        assert ul['Water'] >= min_ul_pipettable, f"Final pool requires {ul['Water']:.2f} uL water, which is less than the minimum pipettable volume."
        pools[-1]['Water'] = ul['Water']

    return pools

def _check_dilution(pools: list[dict[str, float]], num_reads: pd.Series, concs: pd.Series, desired_nm: float = 4):
    volume_so_far = pd.Series(0, index=num_reads.index, dtype=float)
    for pool in pools:
        volume_pool = pd.Series(pool, index=num_reads.index.to_list() + ['Water', 'Prev Pool']).fillna(0)
        ul_prev_pool = volume_pool.pop('Prev Pool')
        if ul_prev_pool > 0:
            volume_pool += volume_so_far * ul_prev_pool / volume_so_far.sum()
        volume_so_far = volume_pool
    concs_final = concs * volume_so_far / volume_so_far.sum()
    concs_final.drop('Water', inplace=True)
    fracs_final = concs_final / desired_nm
    pd.testing.assert_series_equal(fracs_final[num_reads.index], num_reads / num_reads.sum(), check_names=False)

def _check_samples_used_exactly_once(pools: list[dict[str, float]], all_samples: set[str]):
    unused_samples = all_samples.copy()
    for pool in pools:
        samples_in_this_pool = set(pool.keys()) - {'Water', 'Prev Pool'}
        assert samples_in_this_pool.issubset(all_samples), f"Some sample(s) in pool {pool} are not valid samples"
        assert samples_in_this_pool.issubset(unused_samples), f"Some sample(s) in pool {pool} have already been added to a previous pool"
        unused_samples -= samples_in_this_pool
    assert len(unused_samples) == 0, f"Samples {','.join(unused_samples)} have not been added to any pool"
