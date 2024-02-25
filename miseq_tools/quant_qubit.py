import pandas as pd
from .utils import parse_samplesheet, pooled_bp

def qubitquant(samplesheet):
    samples = parse_samplesheet(samplesheet)
    amplicon_sizes = pooled_bp(samples)
    concs = [float(input(f'ng/uL for {sample}: ')) for sample in amplicon_sizes.index]
    concs = pd.Series(concs, index=amplicon_sizes.index)
    concs_molar = 1e6 * concs / (amplicon_sizes * 617.9)
    pd.concat([amplicon_sizes, concs, concs_molar], axis=1, keys=['bp', 'ng/uL', 'nM']).to_csv('quant_qubit.csv')
