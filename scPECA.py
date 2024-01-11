import os
import shutil
import pandas as pd
from PECA_network import GRN
from preprocessing import CPM, FPKM, RNA_process_mode_1, RNA_process_mode_2, RNA_process_mode_3, ATAC_process_mode_1, ATAC_process_mode_2, ATAC_process_mode_3
from homer_motifinding import Split, Homer, Homer_prior, MotifFind
from mf_collect import mf_collect
from mfbs import mfbs, mfbs_c
import pybedtools
import subprocess


class scPECA:
    def __init__(self, path, sample_name, genome):
        self.path = path
        self.sample_name = sample_name
        self.genome = genome
        os.chdir(self.path)
    def RNA_process(self, mode):
        # 在这里根据mode调用其他Python文件中的函数进行数据预处理
        if mode == 1:
            RNA_process_mode_1(self.sample_name, self.genome)
        elif mode == 2:
            RNA_process_mode_2(self.sample_name, self.genome)
        elif mode == 3:
            RNA_process_mode_3(self.sample_name, self.genome)
        else:
            print("Unknown mode")

    def ATAC_process(self, mode):
        # 在这里根据mode调用其他Python文件中的函数进行数据预处理
        if mode == 1:
            ATAC_process_mode_1(self.sample_name)
        elif mode == 2:
            ATAC_process_mode_2(self.sample_name)
        elif mode == 3:
            ATAC_process_mode_3(self.sample_name)
        else:
            print("Unknown mode")

    def network(self, celltype, num_processes=20, prior=0):
        """
        :param num_processes: 线程数
        :param prior: 是否有homer已扫描的文件
        :return: network
        """
        folder_name = celltype
        os.makedirs(folder_name, exist_ok=True)
        os.chdir('./{}'.format(folder_name))
        # openness, region 
        shutil.copy("../{}_PSExp.txt".format(celltype), "{}.txt".format(celltype))
        shutil.copy("../{}_PSOpn.txt".format(celltype), "{}_PSOpn.txt".format(celltype))
        df = pd.read_csv("{}_PSOpn.txt".format(celltype), sep = '\t', header=None, names=['col1', 'col2', 'col3'])
        df[['chr', 'start', 'end']] = df['col1'].str.split('_', expand=True)
        df = df[['chr', 'start', 'end', 'col2']]
        df.to_csv('openness1.bed', sep='\t', header=False, index=False)
        df = pd.read_csv('openness1.bed', sep='\t', header=None)
        sorted_df = df.sort_values(by=[0, 1])
        sorted_df.to_csv('openness1.bed', sep='\t', header=False, index=False)
        sorted_df[3] = sorted_df[0] + '_' + sorted_df[1].astype(str) + '_' + sorted_df[2].astype(str)
        sorted_df[[0, 1, 2, 3]].to_csv('region.txt', sep='\t', header=False, index=False)
        sorted_df[[0, 1, 2]].to_csv('region.bed', sep='\t', header=False, index=False)

        print("step 1: motif binding...")
        MotifFind(self.genome, num_processes, prior)

        print("step2: Opn files...")
        openness1 = pybedtools.BedTool('openness1.bed')
        opn_median = pybedtools.BedTool('../../Prior/Opn_median_{}.bed'.format(self.genome))
        result = openness1.intersect(opn_median, wa=True, wb=True, sorted=True).cut([0, 1, 2, 3, 7]).to_dataframe()
        result.iloc[:, 0] = result.iloc[:, 0].astype(str) + '_' + result.iloc[:, 1].astype(str) + '_' + result.iloc[:,2].astype(str) + '_' + result.iloc[:, 3].astype(str)
        result = result.iloc[:, [0, 4]]
        max_values = result.groupby('chrom')['score'].max().reset_index()
        max_values2 = max_values.iloc[:, 0].str.rsplit('_', 1, expand=True)
        max_values2['score'] = max_values.iloc[:, 1]
        max_values2.to_csv("openness2.bed", sep='\t', header=False, index=False)

        os.mkdir("Enrichment")
        data = pd.read_csv("openness2.bed", sep='\t', header=None)
        data['new_col'] = (data[1] + 0.5) / (data[2] + 0.5)
        sorted_data = data.sort_values('new_col', ascending=False)
        filtered_data = sorted_data[sorted_data.iloc[:, 3] - sorted_data.iloc[:, 2] < 2000].head(10000)
        filtered_data[['chr', 'start', 'end']] = filtered_data.iloc[:, 0].str.split('_', expand=True)
        filtered_data[['chr', 'start', 'end']].to_csv("./Enrichment/region.bed", sep='\t', header=False, index=False)

        print("step3: Motif collect...")
        os.chdir("./Enrichment/")
        if prior==1:
            command = "findMotifsGenome.pl region.bed {} ./. -size given -mask -nomotif -mknown ../../../Data/all_motif_rmdup -preparsedDir /home/zhangjiahao/software/PECA_orig/Homer -p 30".format(self.genome)
            subprocess.run(command, shell=True)
        else:
            command = "findMotifsGenome.pl region.bed {} ./. -size given -mask -nomotif -mknown ../../../Data/all_motif_rmdup -p 30".format(self.genome)
            subprocess.run(command, shell=True)
        mf_collect(self.genome)

        print("step4: Prior...")
        os.chdir('../')
        a = pybedtools.BedTool('region.bed')
        b = pybedtools.BedTool('../../Prior/RE_gene_corr_{}.bed'.format(self.genome))
        result = a.intersect(b, wa=True, wb=True, sorted=True).cut([0, 1, 2, 6, 7, 8]).to_dataframe()
        result.iloc[:, 0] = result.iloc[:, 0].astype(str) + '_' + result.iloc[:, 1].astype(str) + '_' + result.iloc[:,2].astype(str)
        result = result.iloc[:, [0, 3, 4, 5]]
        result.to_csv("peak_gene_100k_corr.bed", sep='\t', header=False, index=False)

        b = pybedtools.BedTool('../../Prior/Enhancer_RE_gene_corr_{}.bed'.format(self.genome))
        result = a.intersect(b, wa=True, wb=True, sorted=True).cut([0, 1, 2, 6, 7, 8]).to_dataframe()
        result.iloc[:, 0] = result.iloc[:, 0].astype(str) + '_' + result.iloc[:, 1].astype(str) + '_' + result.iloc[:,2].astype(str)
        result = result.iloc[:, [0, 3, 4, 5]]
        result.to_csv("peak_gene_100k_corr.bed", mode='a', sep='\t', header=False, index=False)

        print("step5: Network...")
        GRN(celltype,self.genome,num_processes)
        
        # 清理中间文件
        for filename in os.listdir():
            if filename.startswith('.'):
                os.remove(filename)

        print('PECA2 done')



