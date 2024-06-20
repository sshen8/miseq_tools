import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np
from .utils import parse_samplesheet

def demux(samplesheet, stats):
    df_intended = parse_samplesheet(samplesheet)
    df_intended.set_index(["Pool label", "Sample_ID"], inplace=True)
    df_intended = df_intended["Reads (million)"] * 1e6
    df_intended.rename("intended", inplace=True)

    df_actual = json.load(open(stats, "rt"))
    df_actual = pd.concat([pd.DataFrame(df_actual_lane["DemuxResults"], columns=["SampleId", "NumberReads"]).set_index("SampleId").squeeze("columns") for df_actual_lane in df_actual["ConversionResults"]], names=["Lane"], keys=range(1, 1 + len(df_actual["ConversionResults"])), axis="columns")
    df_actual = df_actual.sum(axis="columns").rename("actual")

    df = pd.merge(left=df_intended, right=df_actual, left_on="Sample_ID", right_index=True, how="outer")
    df_pool = df.groupby("Pool label").sum()
    lim_min = min(df_pool['intended'].min(), df_pool['actual'].min())
    lim_max = max(df_pool['intended'].max(), df_pool['actual'].max())
    lim_range = np.log10(lim_max) - np.log10(lim_min)
    lim = np.log10(lim_min) - 0.1 * lim_range, np.log10(lim_max) + 0.1 * lim_range
    lim = 10 ** np.array(lim)
    fig, ax = plt.subplots(figsize=(3, 3))
    ax.plot(np.linspace(*lim), np.linspace(*lim), color="k", linewidth=0.5)
    ax.scatter(df_pool["intended"], df_pool["actual"])
    for label, row in df_pool.iterrows():
        ax.annotate(xy=(row["intended"], row["actual"]), text=label, ha="left", va="center", xytext=(5, 0), textcoords="offset points", fontsize=6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("Intended num reads")
    ax.set_ylabel("Actual num reads")
    fig.savefig("demux_stats.pdf", bbox_inches="tight")
    plt.close(fig)

    print((df_pool["actual"] / df_pool["intended"]).rename("actual/intended"))
