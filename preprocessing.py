import pandas as pd
import os
"""
This module contains functions for working with many different types of data
"""
def FPKM(df):
   """
   input a dataframe with 3 columns (gene name, count, length)
   :return: gene name, FPKM
   """
   df.columns = ['Gene', 'Count', 'Length']
   total_count = df['Count'].sum()
   df['FPKM'] = (df['Count'] * 10 ** 9) / (total_count * df['Length'])
   output_df = df[['Gene', 'FPKM']]
   return output_df

def CPM(df):
    """
    input a dataframe with 2 columns (peak, count)
    :return: peak, CPM
    """
    df.columns = ['Count']
    total_count = df['Count'].sum()
    df['CPM'] = (df['Count'] * 10 ** 6) / total_count
    output_df = df[['CPM']].copy()
    output_df.index = df.index
    return output_df

def RNA_process_mode_1(sample_name, genome):
    """
    bulk RNA-seq count
    :return: RNA FPKM
    """
    rna_df = pd.read_csv('{}_RNA.txt'.format(sample_name),sep='\t', header=None)
    trans_df = pd.read_csv('../Data/{}_gene_length.txt'.format(genome), sep='\t', header=None)
    merged_df = pd.merge(rna_df, trans_df, on=0)
    rna_FPKM = FPKM(merged_df)
    rna_FPKM = rna_FPKM[rna_FPKM.iloc[:, 1] > 0.1]
    rna_FPKM.to_csv('{}_PSExp.txt'.format(sample_name), sep='\t', header=False, index=False)

def RNA_process_mode_2(sample_name, genome):
    """
    scRNA-seq count without meta
    :return: pseudo bulk RNA FPKM
    """
    scrna_df = pd.read_csv('{}_scRNA.csv'.format(sample_name),index_col=0)
    scrna_df = scrna_df.sum(axis=1).to_frame()
    trans_df = pd.read_csv('../Data/{}_gene_length.txt'.format(genome), sep='\t', header=None)
    merged_df = scrna_df.merge(trans_df, left_index=True, right_on=0)
    merged_df = merged_df.iloc[:, [0, 1, 3]]
    scrna_FPKM = FPKM(merged_df)
    scrna_FPKM = scrna_FPKM[scrna_FPKM.iloc[:, 1] > 0.1]
    scrna_FPKM.to_csv('{}_PSExp.txt'.format(sample_name), sep='\t', header=False, index=False)

def RNA_process_mode_3(sample_name, genome):
    """
    scRNA-seq count with meta
    :return: pseudo bulk RNA FPKM for different cell types
    """
    scrna_df = pd.read_csv('{}_scRNA.csv'.format(sample_name), index_col=0)
    meta = pd.read_csv('{}_scRNA_meta.csv'.format(sample_name), index_col=0)
    meta.columns = ['celltype']
    celltype_name = meta['celltype'].value_counts().reset_index()
    trans_df = pd.read_csv('../Data/{}_gene_length.txt'.format(genome), sep='\t', header=None)
    for i in range(celltype_name.shape[0]):
        barcode_index = meta.loc[meta['celltype'] == celltype_name.iloc[0,0]].index
        scrna_df2 = scrna_df[barcode_index]
        scrna_df2 = scrna_df2.sum(axis=1).to_frame()
        merged_df = scrna_df2.merge(trans_df, left_index=True, right_on=0)
        merged_df = merged_df.iloc[:, [0, 1, 3]]
        scrna_FPKM = FPKM(merged_df)
        scrna_FPKM = scrna_FPKM[scrna_FPKM.iloc[:, 1] > 0.1]
        scrna_FPKM.to_csv('{}_PSExp.txt'.format(str(celltype_name.iloc[0,0])), sep='\t', header=False, index=False)

def ATAC_process_mode_1(sample_name):
    """
    bulk ATAC-seq count
    :return: ATAC CPM
    """
    atac_df = pd.read_csv('{}_ATAC.txt'.format(sample_name),sep='\t', header=None)
    atac_CPM = CPM(atac_df)
    atac_CPM = atac_CPM[atac_CPM.iloc[:, 0] > 1]
    atac_CPM.to_csv('{}_PSOpn.txt'.format(sample_name), sep='\t', header=False, index=False)

def ATAC_process_mode_2(sample_name):
    """
    scRNA-seq count without meta
    :return: pseudo bulk RNA FPKM
    """
    scatac_df = pd.read_csv('{}_scATAC.csv'.format(sample_name),index_col=0)
    scatac_df = scatac_df.sum(axis=1).to_frame()
    scatac_CPM = CPM(scatac_df)
    scatac_CPM = scatac_CPM[scatac_CPM.iloc[:, 0] > 1]
    scatac_CPM.to_csv('{}_PSOpn.txt'.format(sample_name), sep='\t', header=False, index=True)

def ATAC_process_mode_3(sample_name):
    """
    scATAC-seq count with meta
    :return: pseudo bulk ATAC CPM for different cell types
    """
    scatac_df = pd.read_csv('{}_scATAC.csv'.format(sample_name), index_col=0)
    meta = pd.read_csv('{}_scATAC_meta.csv'.format(sample_name), index_col=0)
    meta.columns = ['celltype']
    celltype_name = meta['celltype'].value_counts().reset_index()
    for i in range(celltype_name.shape[0]):
        barcode_index = meta.loc[meta['celltype'] == celltype_name.iloc[0,0]].index
        scatac_df2 = scatac_df[barcode_index]
        scatac_df2 = scatac_df2.sum(axis=1).to_frame()
        scatac_CPM = CPM(scatac_df2)
        scatac_CPM = scatac_CPM[scatac_CPM.iloc[:, 0] > 0.1]
        scatac_CPM.to_csv('{}_PSOpn.txt'.format(str(celltype_name.iloc[0,0])), sep='\t', header=False, index=False)





