import pandas as pd
import matplotlib.pyplot as plt

def quant_combine(kapa_fname, qubit_fname):

    kapa = pd.read_csv(kapa_fname, index_col=0)
    qubit = pd.read_csv(qubit_fname, index_col=0)

    fig, ax = plt.subplots()
    ax.scatter(qubit['nM'], kapa['nM'])
    ax.set_xlabel('Qubit nM')
    ax.set_ylabel('KAPA nM')
    lim_max = 1.1 * max(kapa['nM'].max(), qubit['nM'].max())
    ax.set_xlim(0, lim_max)
    ax.set_ylim(0, lim_max)
    ax.set_aspect('equal')
    ax.plot([0, lim_max], [0, lim_max], 'k')
    for (i, nM), (j, nM2) in zip(qubit["nM"].items(), kapa["nM"].items()):
        ax.annotate(xy=(nM, nM2), text=i, xytext=(5, 0), textcoords='offset points', ha='left', va='center')
        assert i == j
    fig.savefig('quant_combined.pdf', bbox_inches='tight')
    plt.close(fig)

    df_out = kapa['bp'].to_frame()
    assert (kapa['bp'] == qubit['bp']).all()
    df_out[['ng/uL', 'nM']] = ((kapa[['ng/uL', 'nM']] + qubit[['ng/uL', 'nM']]) / 2)

    df_out.to_csv('quant_combined.csv')
