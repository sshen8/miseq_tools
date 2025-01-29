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
    assert 'Water' not in ul
    ul['Water'] = min_ul_total - ul.sum()

    pools = []
    dilution_factors = []
    ul_prev_pool = None
    while any(ul < min_ul_pipettable):
        dilution_factor = min_ul_pipettable / min(ul)
        if ul_prev_pool is not None:
            dilution_factor = max(dilution_factor, min_ul_pipettable / ul_prev_pool)
        ul_diluted = ul * dilution_factor
        assert all(ul_diluted >= min_ul_pipettable)
        samples_in_this_pool = ul_diluted.index[ul_diluted <= max_ul_pipettable]
        volumes_in_this_pool = ul_diluted[samples_in_this_pool].to_dict()
        if ul_prev_pool is not None:
            volumes_in_this_pool |= {'Prev Pool': ul_prev_pool * dilution_factor}
        pools.append(volumes_in_this_pool)
        dilution_factors.append(dilution_factor)
        ul_prev_pool = ul[samples_in_this_pool].sum()
        ul.drop(samples_in_this_pool, inplace=True)

    volumes_in_this_pool = ul.to_dict()
    if ul_prev_pool is not None:
        volumes_in_this_pool |= {'Prev Pool': ul_prev_pool}
    pools.append(volumes_in_this_pool)
    dilution_factors.append(1)

    return pools
