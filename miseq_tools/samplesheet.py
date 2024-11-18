import pandas as pd
import logging
import datetime
from .utils import parse_samplesheet
import os
import Bio.Seq

def check_indexes(df):
    dir_known_i7 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'known_barcodes', 'i7')
    known_i7 = pd.concat([pd.read_csv(os.path.join(dir_known_i7, fname), header=None, index_col=0, delimiter="\t" if fname.endswith('.tsv') else ",").squeeze() for fname in os.listdir(dir_known_i7)])
    dir_known_i5 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'known_barcodes', 'i5')
    known_i5 = pd.concat([pd.read_csv(os.path.join(dir_known_i5, fname), header=None, index_col=0, delimiter="\t" if fname.endswith('.tsv') else ",").squeeze() for fname in os.listdir(dir_known_i5)])
    # i7 rev comp
    if (errors := df['index'].isin(known_i7.apply(Bio.Seq.reverse_complement))).any():
        logging.warning(f'i7 might need to be reverse complemented for these samples: {df.loc[errors, "Sample_ID"].tolist()}')
    # i5 rev comp
    if (errors := df['index2'].isin(known_i5.apply(Bio.Seq.reverse_complement))).any():
        logging.warning(f'i5 might need to be reverse complemented for these samples: {df.loc[errors, "Sample_ID"].tolist()}')
    # i7 and i5 swapped
    if (errors := df['index'].isin(known_i5) | df['index2'].isin(known_i7)).any():
        logging.warning(f'i7 and i5 might be swapped for these samples: {df.loc[errors, "Sample_ID"].tolist()}')

def format_samplesheet(fname_in, fname_out, nextseq=False):
    df = parse_samplesheet(fname_in)

    read_info = pd.read_csv(fname_in, usecols=[9, 10], header=None, nrows=4, index_col=0).squeeze()

    # check validity
    if (errors := df['Sample_ID'].duplicated()).any():
        logging.warning(f'Duplicate Sample_ID found: {df.loc[errors, "Sample_ID"].tolist()}')
    if (errors := df[['index', 'index2']].apply(tuple, axis=1).duplicated()).any():
        logging.warning(f'Duplicate index pair found: {df.loc[errors, ["index", "index2"]].values.tolist()}')
    if (errors := df['Sample_ID'].str.len() > 40).any():
        logging.warning(f'Sample_ID too long: {df.loc[errors, "Sample_ID"].tolist()}')
    if (errors := ~df['Sample_ID'].str.match(r'^[a-zA-Z0-9-_]+$')).any():
        logging.warning(f'Invalid characters in Sample_ID: {df.loc[errors, "Sample_ID"].tolist()}')
    # make sure index is the right length
    if (errors := df['index'].str.len() != read_info['Index 1 (i7)']).any():
        logging.warning(f'index is the wrong length for these samples: {df.loc[errors, "Sample_ID"].tolist()}')
    if (errors := df['index2'].str.len() != read_info['Index 2 (i5)']).any():
        logging.warning(f'index2 is the wrong length: {df.loc[errors, "Sample_ID"].tolist()}')
    # check against known indexes
    check_indexes(df)

    # reverse complement i5 for nextseq if necessary
    if nextseq:
        df['index2'] = df['index2'].apply(Bio.Seq.reverse_complement)

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
