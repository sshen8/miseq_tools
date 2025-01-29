import pandas as pd
from .utils import parse_samplesheet, pooled_reads

def pooling(samplesheet, quant_csv, pools: int):
    samples = parse_samplesheet(samplesheet)
    reads = (pooled_reads(samples) * 1e6).astype(int)
    concs = pd.read_csv(quant_csv, index_col=0)["nM"]
    pools = _pools(reads, concs)
    for i, pool in enumerate(pools, 1):
        print(f"Pool {i}")
        for sample, ul in pool.items():
            print(f"{sample}: {ul:.2f} uL")

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
    ul['Water 1'] = min_ul_total - ul.sum()

    # initial assign dilution factor groups
    dilution_factors = pd.Series(None, index=ul.index, dtype=float)
    while any(dilution_factors.isna()):
        remaining_samples = dilution_factors[dilution_factors.isna()].index
        dilution_factor = min_ul_pipettable / min(ul[remaining_samples])
        if dilution_factor <= 1:
            dilution_factors[remaining_samples] = 1
            break
        hypothetical_volumes = ul[remaining_samples] * dilution_factor
        samples_in_this_pool = hypothetical_volumes.index[hypothetical_volumes <= max_ul_pipettable]
        dilution_factors[samples_in_this_pool] = dilution_factor

    # calculate volumes per pool
    pools = []
    dilution_factors_list = sorted(dilution_factors.unique(), reverse=True)
    for dilution_factor in dilution_factors_list:
        samples_in_previous_pool = dilution_factors[dilution_factors > dilution_factor].index
        if len(samples_in_previous_pool) == 0:
            ul_prev_pool = dict()
        else:
            ul_prev_pool = {'Prev Pool': ul[samples_in_previous_pool].sum() * dilution_factor}
        samples_in_this_pool = dilution_factors[dilution_factors == dilution_factor].index
        ul_diluted = ul[samples_in_this_pool] * dilution_factor
        pools.append({
            **ul_diluted.to_dict(),
            **ul_prev_pool,
        })

    return pools
