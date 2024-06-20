import pandas as pd
from .utils import parse_samplesheet, pooled_reads
import sklearn.cluster
import numpy as np

samplesheet = "SPS246 miseq - Sheet1.csv"
quant_csv = "quant_combined.csv"

def pooling(samplesheet, quant_csv, pools: int):
    samples = parse_samplesheet(samplesheet)
    reads = pooled_reads(samples)
    concs = pd.read_csv(quant_csv, index_col=0)["nM"]

    frac_reads = reads / reads.sum()
    desired_nm = frac_reads * 4
    desired_ul = 10
    min_ul_pipettable = 2
    ul = desired_ul * desired_nm / concs
    model = sklearn.cluster.KMeans(n_clusters=pools)
    clusters = pd.Series(model.fit_predict(ul.values.reshape(-1, 1)), index=ul.index)
    ul_water = desired_ul - ul.sum()

    # smallest value per cluster
    small = ul.groupby(clusters).min()
    factors = np.maximum(min_ul_pipettable / small, 1).sort_values(ascending=False)
    factors.rename("factor", inplace=True)
    factors.rename_axis("cluster", inplace=True)

    prev_clusters = []
    for pool_i, (cluster_i, factor) in enumerate(factors.items()):
        cur_ul = ul[clusters == cluster_i] * factor
        if prev_clusters:
            cur_ul[f"Pool {pool_i}"] = ul[clusters.isin(prev_clusters)].sum() * factor
        # anticipate whether the next pool is gonna use less than min_ul_pipettable
        if pool_i < small.size - 1:
            next_factor = factors.iloc[pool_i + 1]
            next_pool = ul[clusters.isin(prev_clusters + [cluster_i])].sum() * next_factor
        if pool_i == small.size - 1:
            cur_ul["Water"] = ul_water * factor
        print(f'=======Pool {pool_i + 1:<2d} (total = {cur_ul.sum():.2f} uL, factor = {factor:.4f})=======')
        for sample, u in cur_ul.items():
            print(f'{sample:14s}: {u:.2f} uL')
        prev_clusters.append(cluster_i)
