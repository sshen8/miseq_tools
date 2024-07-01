import pandas as pd
import os
import numpy as np
import scipy.stats
import logging
from .utils import parse_samplesheet, pooled_bp
import matplotlib.pyplot as plt
import seaborn as sns

def kapaquant(kapafolder, samplesheet, dilution, standard_bp: int):
    samples = parse_samplesheet(samplesheet)
    amplicon_sizes = pooled_bp(samples)

    fname = next(filter(lambda x: x.endswith('.csv') and 'Quantification Summary' in x, os.listdir(kapafolder)))
    if not fname:
        logging.error(f'Could not find quantification summary data in {kapafolder}')
        return

    df = pd.read_csv(os.path.join(kapafolder, fname))
    df.drop(columns=df.columns[0], inplace=True)
    df.drop(index=df.index[df.Content == 'Unkn'], inplace=True)
    df.Content = df.Content.where(df.Content != 'Std', other=df.Content + df.Well.str.slice(0, 1))
    df.sort_values(by=['Content', 'Well'], inplace=True)

    # standard curve
    std = df[df.Content.str.startswith('Std')]
    data = std.groupby('Content')[['Cq', 'SQ']].mean()
    if not (std.groupby('Content')['SQ'].std() == 0).all():
        logging.error('Concentration is not constant for all replicates of standards')
    deltaCq = data.Cq.diff()
    logging.info(f'ΔCq: {", ".join(deltaCq[1:].apply("{:.2f}".format))}')
    if not deltaCq[1:].between(3.1, 3.6).all():
        logging.warning('ΔCq is not within 3.1-3.6 for all standards')

    x = np.log10(data.SQ * 1e12) # pM
    y = data.Cq
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    efficiency = 10**(-1/slope) - 1
    logging.info(f'Efficiency: {100 * efficiency:.2f}%')
    if not 0.9 <= efficiency <= 1.1:
        logging.warning(f'Efficiency is {100 * efficiency:.2f}%, expected 90-110%')
    logging.info(f'Slope: {slope:.4f}')
    logging.info(f'R²: {r_value**2:.4f}')
    if r_value**2 < 0.99:
        logging.warning(f'R² is {r_value**2:.4f}, expected >0.99')
    logging.info(f'Intercept: {intercept:.3f}')

    fig, ax = plt.subplots()
    ax.scatter(x, y)
    ax.plot(x, slope * x + intercept, 'r')
    ax.set_xlabel('log10(pM)')
    ax.set_ylabel('Cq')
    ax.set_title('Standard curve')
    fig.savefig('quant_kapa_standards.pdf', bbox_inches='tight')
    plt.close(fig)

    unkn = df[df.Content.str.startswith('Unkn')].iloc[:len(amplicon_sizes) * 3]
    avg_cq = unkn.groupby('Content')['Cq'].mean()
    avg_cq.index = amplicon_sizes.index
    conc = np.power(10, (avg_cq - intercept) / slope) # pM
    conc_size_adjusted = conc * standard_bp / amplicon_sizes # pM
    conc_undiluted = conc_size_adjusted * dilution / 1e3 # nM
    conc_undiluted_mass = conc_undiluted * amplicon_sizes * 617.9 * 1e-6 # ng/uL

    fig, ax = plt.subplots()
    pools = amplicon_sizes.index.to_series(index=[f'Unkn-{i:02d}' for i in range(1, len(amplicon_sizes) + 1)])
    unkn_plot = unkn.join(pools, on="Content")
    std_plot = std.copy()
    std_plot["Pool label"] = "Standards"
    data_plot = pd.concat([std_plot, unkn_plot])
    sns.swarmplot(data=data_plot, x='Pool label', y='Cq', ax=ax, order=["Standards"] + pools.tolist())
    ax.set_xlabel("")
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
        tick.set_ha('right')
    range_cq = std.Cq.min(), std.Cq.max()
    range_cq_diff = range_cq[1] - range_cq[0]
    ax.set_ylim(range_cq[0] - 0.1 * range_cq_diff, range_cq[1] + 0.1 * range_cq_diff)
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    # second y-axis
    ax2 = ax.twinx()
    ax2.set_ylabel('Concentration (pM)')
    range_sq = std.SQ.max(), std.SQ.min()
    range_sq_log = np.log10(range_sq) + 12 # convert M to pM
    range_sq_log_diff = range_sq_log[0] - range_sq_log[1]
    ax2.set_ylim(np.power(10, [range_sq_log[0] + 0.1 * range_sq_log_diff, range_sq_log[1] - 0.1 * range_sq_log_diff]))
    ax2.set_yscale("log")
    ax2.spines['right'].set_visible(True)
    fig.savefig('quant_kapa.pdf', bbox_inches='tight')
    plt.close(fig)

    pd.concat([amplicon_sizes, conc_undiluted_mass, conc_undiluted], axis=1, keys=['bp', 'ng/uL', 'nM']).to_csv('quant_kapa.csv')
