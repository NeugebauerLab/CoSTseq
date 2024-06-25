import pandas as pd
import numpy as np
import argparse
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
import re
import pickle

def run_deseq2(counts_table, ctrl, threads, output_name):
    counts_df = pd.read_csv(counts_table, delimiter='\t', skiprows=1, index_col='Geneid')

    # rename columns
    new_col_names = []
    samples = []
    for i in list(counts_df.columns):
        if len(i.split('/')) == 1:
            new_col_names.append(i)
        else:
            new_col_names.append(re.findall(r"\S{4}_\d", i.split('/')[-1])[0])
            samples.append(re.findall(r"\S{4}_\d", i.split('/')[-1])[0])
            
    counts_df.rename(columns=dict(zip(list(counts_df.columns), new_col_names)), inplace=True)

    # run DEseq2 for all samples, compare each to control
    results = {}
    for cond in np.unique([i.split('_')[0] for i in samples]):
        # skip the control
        if cond == ctrl:
            continue
        
        # collect data for ctrl-sample pair
        sample_list = [i for i in samples if ctrl in i] + [i for i in samples if cond in i]
        sample_df  = counts_df[sample_list].T
        conditions = pd.DataFrame({'samples': sample_list, 
                                   'condition': ['A']*len([i for i in samples if ctrl in i]) + ['B']*len([i for i in samples if cond in i])}).set_index('samples')
        
        # run DEseq2
        dds = DeseqDataSet(counts=sample_df, metadata=conditions, design_factors="condition", refit_cooks=True, n_cpus=threads)
        dds.deseq2()
        stat_res = DeseqStats(dds)
        stat_res.summary()
        results[cond] = stat_res.results_df.reset_index()

    with open(output_name, 'wb') as f:
        pickle.dump(results, f)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('counts_table', help='Output from featureCounts package. Sample names must be four characters followed by an underscore and a number.')
    parser.add_argument('-ctrl', help='Wildtype sample ID; sample IDs must be 4 characters.')
    parser.add_argument('-c', help='Number of threads.', type=int)
    parser.add_argument('-o', help='Output file name.')
    args = parser.parse_args()

    run_deseq2(args.counts_table, args.ctrl, args.c, args.o)
