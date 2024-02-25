import pandas as pd
import warnings
import datetime
from .utils import parse_samplesheet

def format_samplesheet(fname_in, fname_out):
    df = parse_samplesheet(fname_in)

    read_info = pd.read_csv(fname_in, usecols=[9, 10], header=None, nrows=4, index_col=0).squeeze()

    # check validity
    if (errors := df['Sample_ID'].duplicated()).any():
        warnings.warn(f'Duplicate Sample_ID found: {df.loc[errors, "Sample_ID"].tolist()}')
    if (errors := df[['index', 'index2']].apply(tuple, axis=1).duplicated()).any():
        warnings.warn(f'Duplicate index pair found: {df.loc[errors, ["index", "index2"]].values.tolist()}')
    if (errors := df['Sample_ID'].str.len() > 40).any():
        warnings.warn(f'Sample_ID too long: {df.loc[errors, "Sample_ID"].tolist()}')
    if (errors := ~df['Sample_ID'].str.match(r'^[a-zA-Z0-9-_]+$')).any():
        warnings.warn(f'Invalid characters in Sample_ID: {df.loc[errors, "Sample_ID"].tolist()}')
    # make sure index is the right length
    if (errors := df['index'].str.len() != read_info['Index 1 (i7)']).any():
        warnings.warn(f'index is the wrong length for these samples: {df.loc[errors, "Sample_ID"].tolist()}')
    if (errors := df['index2'].str.len() != read_info['Index 2 (i5)']).any():
        warnings.warn(f'index2 is the wrong length: {df.loc[errors, "Sample_ID"].tolist()}')

    # write out
    with open(fname_out, 'wt') as f:
        f.write(f"""[Header],,,,
IEMFileVersion,4,,,
Date,{datetime.date.today().strftime("%-m/%-d/%y")},,,
Workflow,GenerateFASTQ,,,
Application,FASTQ Only,,,
Assay,TruSeq HT,,,
Description,,,,
Chemistry,Amplicon,,,
,,,,
[Reads],,,,
{read_info["Read 1"]},,,,
{read_info["Read 2"]},,,,
,,,,
[Settings],,,,
,,,,
[Data],,,,
Sample_ID,I7_Index_ID,index,I5_Index_ID,index2
""")
        df[['Sample_ID', 'I7_Index_ID','index', 'I5_Index_ID', 'index2']].to_csv(f, index=False, header=False, mode='a')
