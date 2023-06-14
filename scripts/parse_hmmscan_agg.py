import pandas as pd
import argparse
import os, sys, csv
from glob import glob
import numpy as np
from gene_list import quinones, quinones_full_profile
from quinones_scores_cutoff import cutoffs


def get_genomes_list(dataset_path):
     genomes = [os.path.basename(genome) for genome in glob(dataset_path+"/GC*")]
     return genomes


def read_hmmscan_table(file):
    dataframe = pd.read_csv(file, header=None, delim_whitespace=True,
                            names = ['genome', 'target', '2','tlen','query',
                            '5','qlen','7','8','9','10','11','12','i-evalue',
                            'score','15','hmm_from','hmm_to','ali_from',
                            'ali_to','20','21','22','23'])
    return dataframe


def best_hit(dataframe):
    ind = dataframe.groupby(['genome', 'query'])['i-evalue'].idxmin()
    ind = ind.dropna()
    #print(ind)
    df_besthit = dataframe.loc[ind].sort_index()
    return df_besthit


def coverage(dataframe):
    end = dataframe['hmm_to']
    start = dataframe['hmm_from']
    tlen = dataframe['tlen'].astype(int)
    cov = (end - start + 1) / tlen
    dataframe["coverage"] = cov
    return dataframe


def filter_ieval_cov(dataframe, ieval_thresh, cov_thresh):
    df_ieval = dataframe['i-evalue'] <= ieval_thresh
    df_cov = dataframe['coverage'] >= cov_thresh
    dataframe_filtered = dataframe[(df_ieval) & (df_cov)]
    return dataframe_filtered


def apply_cutga(dataframe, cut_off_dict):
    for key, value in cut_off_dict.items():
        df_target = dataframe['target']==key
        print(key)
        df_score = dataframe['score']<=value
        dataframe = dataframe.drop(dataframe[(df_target) & (df_score)].index)
    return dataframe


def filt_on_genes(dataframe):
    return dataframe[dataframe['target'].isin(quinones_full_profile)]


def get_interm_table(dataframe, dataset):
    dataframe = filt_on_genes(dataframe)
    df_interm = dataframe[['genome', 'target', 'query', 'i-evalue', 'score']]
    df_interm.to_csv('results/'+dataset+'/'+dataset+'_intermediate_results_hmmscan.tsv', sep='\t')


def get_results(dataframe, dataset_path):
    genomes = get_genomes_list(dataset_path)
    df_quin_genes = filt_on_genes(dataframe)
    df_genomes = pd.DataFrame(genomes, columns=['genome'])
    column_names = quinones_full_profile.copy()
    column_names.insert(0, 'genome')
    df_pivoted = pd.pivot_table(df_quin_genes[['genome', 'target']],
                                columns="target", aggfunc="size",
                                index="genome")
    df_pivoted_reindexed = df_pivoted.reset_index()
    df_pivoted_reindexed = df_pivoted_reindexed.rename_axis(None, axis=1)
    df_pivoted_reindexed = df_pivoted_reindexed.reindex(columns = column_names)
    df_final = pd.merge(df_genomes, df_pivoted_reindexed, on="genome", how='outer')
    df_final = df_final.fillna(0)
    return df_final


def extract_id(dataframe, col):
    """
    From the given column (usually FTP ones), extract the id of the genome.
    This id must correspond to the genome directory name.
    """
    return dataframe[col].str.extract(r'(.*/)(GC._.*)')


def get_filenames(df_tax):
    """
    """
    try:
        extracted_id_refseq = extract_id(df_tax, 'RefSeq FTP')
        extracted_id_genbank = extract_id(df_tax, 'GenBank FTP')
        df_tax['genome'] = extracted_id_refseq[1]
        df_tax['genome'] = df_tax['genome'].fillna(extracted_id_genbank[1])
    except KeyError:
        return df_tax
    return df_tax


def add_taxonomy(dataframe, dataframe_tax):
    """
    Merge the dataframe (either count or proteins) &
    the taxonomy dataframe
    """
    dataframe_tax = get_filenames(dataframe_tax)
    #print(dataframe_tax['genome'])
    dataframe_merged = dataframe.merge(dataframe_tax, on = ['genome'], how="inner")
    #print(dataframe_merged)
    return dataframe_merged


def merge_col(dataframe):
    """
    Merge columns:
    - UbiT UbiT2
    - MqnM MqnM2
    - UbiX UbiX2
    - RquA RquA2
    - UbiFHILM_Cyano FMO_Cyano
    """
    dataframe['UbiT'] = dataframe['UbiT']+dataframe['UbiT2']
    dataframe['MqnM'] = dataframe['MqnM']+dataframe['MqnM2']
    dataframe['UbiX'] = dataframe['UbiX']+dataframe['UbiX2']
    dataframe['RquA'] = dataframe['RquA']+dataframe['RquA2']
    dataframe['UbiFHILM'] = dataframe['UbiFHILM_Cyano']+dataframe['FMO_Cyano']
    col_names = quinones.copy()
    col_names.insert(0, 'genome')
    return dataframe[col_names]


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-F', '--file', required=False, help="best_solution results concatenated \
                        into a single best_solutions file")
    parser.add_argument('-O', '--output', required=False)
    parser.add_argument('-T', '--tax', required=False, help="taxonomy file")
    parser.add_argument('-i', '--ieval', required=True, help="i-eval cutoff")
    parser.add_argument('-c', '--cov', required=True, help="coverage cutoff")

    args = parser.parse_args()

    FILE = args.file
    OUT = args.output
    IEVAL = float(args.ieval)
    COV = float(args.cov)
    DS = OUT.split("/")[1]
    DS_PATH = "data/"+OUT

    DF = read_hmmscan_table(FILE)
    DF_BH = best_hit(DF)
    DF_BH = coverage(DF_BH)
    DF_BH_FILT = filter_ieval_cov(DF_BH, IEVAL, COV)
    DF_BH_GA = apply_cutga(DF_BH_FILT, cutoffs)
    get_interm_table(DF_BH_GA, DS)
    DF_QUIN = get_results(DF_BH_GA, DS_PATH)
    DF_QUIN_CLEAN = merge_col(DF_QUIN)
    DF_QUIN_CLEAN.to_csv(OUT, sep='\t')
    ### add taxonomy
    if args.tax:
        OUT_TAX = OUT.replace('.tsv', '_tax.tsv')
        TAX = pd.read_csv(args.tax, low_memory=False)
        print(TAX)
        DF_TAX = add_taxonomy(DF_QUIN_CLEAN, TAX)
        DF_TAX.to_csv(OUT_TAX, sep='\t')

