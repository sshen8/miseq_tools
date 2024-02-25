import pandas as pd

def parse_samplesheet(f):
    df = pd.read_csv(f, usecols=range(8))
    df.dropna(how='all', inplace=True)
    df.drop(df.index[df['Pool label'].str.upper().isin(('PHIX', 'TOTAL'))], inplace=True)
    df['Pool label'] = df['Pool label'].ffill()
    df['Amplicon size (bp)'] = df['Amplicon size (bp)'].astype(float)
    return df

def pooled_bp(samples):
    def _weighted_average_size(x):
        return (x['Reads (million)'] * x['Amplicon size (bp)']).sum() / x['Reads (million)'].sum()
    return samples.groupby('Pool label', sort=False).apply(_weighted_average_size, include_groups=False)

def pooled_reads(samples):
    return samples.groupby('Pool label', sort=False)['Reads (million)'].sum()
