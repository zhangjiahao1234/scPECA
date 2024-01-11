import numpy as np
import pandas as pd
from ismember import ismember
import scipy.sparse as sparse
import os
import scipy.io as scio
from numpy.matlib import repmat
import numpy_groupies as npg
from mfbs import mfbs, mfbs_c

def GRN(celltype, genome,num_processes):
	if "hg" in genome:
		species = "human"
	elif "mm" in genome:
		species = "mouse"
	C = pd.read_csv('openness2.bed', sep='\t', header=None)
	Element_name = C.iloc[:, 0]
	Opn = C.iloc[:, 1]
	Opn_median = C.iloc[:, 2]

	# 读取mat文件
	# Match2, motifName, motifWeight
	MotifMatch_mouse_rmdup = scio.loadmat('../../Data/MotifMatch_{}_rmdup.mat'.format(species))
	# Exp_median, List, R2, TFExp_median, TFName
	TFTG_corr_mouse = scio.loadmat('../../Prior/TFTG_corr_{}.mat'.format(species))

	Match2 = np.empty((MotifMatch_mouse_rmdup['Match2'].shape[0], 2), dtype='U100')
	for i in range(MotifMatch_mouse_rmdup['Match2'].shape[0]):
		Match2[i, 0] = MotifMatch_mouse_rmdup['Match2'][i, 0].item()
		Match2[i, 1] = MotifMatch_mouse_rmdup['Match2'][i, 1].item()

	motifName = np.empty((MotifMatch_mouse_rmdup['motifName'].shape[0], 1), dtype='U100')
	for i in range(MotifMatch_mouse_rmdup['motifName'].shape[0]):
		motifName[i, 0] = MotifMatch_mouse_rmdup['motifName'][i, 0].item()

	motifWeight = MotifMatch_mouse_rmdup['motifWeight']

	TFName = np.empty((TFTG_corr_mouse['TFName'].shape[0], 1), dtype='U100')
	for i in range(TFTG_corr_mouse['TFName'].shape[0]):
		TFName[i, 0] = TFTG_corr_mouse['TFName'][i, 0].item()

	List = np.empty((TFTG_corr_mouse['List'].shape[0], 1), dtype='U100')
	for i in range(TFTG_corr_mouse['List'].shape[0]):
		List[i, 0] = TFTG_corr_mouse['List'][i, 0].item()

	R2 = TFTG_corr_mouse['R2']
	Exp_median = TFTG_corr_mouse['Exp_median']
	# ---------------------------

	N = num_processes
	TF_binding = mfbs_c(N,TFName, Element_name, motifName, motifWeight, Match2)

	# gene expr
	C = pd.read_csv('{}.txt'.format(celltype), sep='\t', header=None)
	Symbol = C.iloc[:, 0]
	G = C.iloc[:, 1]

	# CR-TF
	CRInfo_mouse = scio.loadmat('../../Data/CRInfo_{}.mat'.format(species))
	C_TFName = np.empty((CRInfo_mouse['C_TFName'].shape[0]), dtype='U100')
	for i in range(CRInfo_mouse['C_TFName'].shape[0]):
		C_TFName[i] = CRInfo_mouse['C_TFName'][i, 0].item()
	CR_TF = CRInfo_mouse['CR_TF']
	CRName = np.empty((CRInfo_mouse['CRName'].shape[0]), dtype='U100')
	for i in range(CRInfo_mouse['CRName'].shape[0]):
		CRName[i] = CRInfo_mouse['CRName'][i, 0].item()
	TFS = CRInfo_mouse['TFS']

	eita0 = -30.4395
	eita1 = 0.8759
	d, f = ismember(C_TFName, TFName)
	C_TFName = C_TFName[d]
	TFS = TFS[d]
	CR_TF = CR_TF[:, d]
	TFB = TF_binding[f].todense()
	C_TFExp = np.zeros(len(C_TFName))
	d, f = ismember(C_TFName, Symbol)
	C_TFExp[d] = np.log2(1 + G[f])

	d1, f1 = ismember(C_TFName, TFName)

	TFBO = np.power(np.multiply(np.multiply(repmat(C_TFExp, len(Opn), 1).T * repmat(C_TFExp / np.squeeze(TFS), len(Opn), 1).T,np.squeeze(TF_binding.todense()[np.where(d1 == True), :])), repmat(Opn, len(C_TFName), 1)), 0.25)
	CRB = eita0 + eita1 * np.dot(CR_TF, TFBO)

	alhfa = 0.5
	Opn_median = np.log2(1 + Opn_median)
	Opn1 = np.log2(1 + Opn)
	Opn = Opn1 * (Opn1 / (Opn_median + 0.5))

	geneName = np.intersect1d(List, Symbol)
	d, f = ismember(geneName, List)
	R2 = R2[:, f]
	Exp_median = Exp_median[f]

	d, f = ismember(geneName, Symbol)
	G = G[f]
	d1 = sorted(G)
	f1 = np.argsort(np.array(G))
	d2 = sorted(Exp_median)
	f2 = np.argsort(np.squeeze(Exp_median))
	G1 = np.empty(len(G))
	for i in range(len(G)):
		G1[f1[i]] = d2[i]
	G = np.multiply((np.power(G1, alhfa)), (G1 / (np.squeeze(Exp_median) + 0.5)))

	d, f = ismember(TFName, geneName)
	TFName = TFName[np.where(d == True)]
	TF_binding = TF_binding[np.where(d == True)[0], :]
	TFExp = G[f]
	R2 = R2[np.where(d == True)[0], :]

	C = pd.read_csv('./Enrichment/knownResults_TFrank.txt', header=None, sep='\t')
	d, f = ismember(TFName, C.iloc[:, 0])
	TF_motif = np.zeros(len(TFName))
	TF_motif[np.where(d == True)[0]] = C.iloc[:, 1][f]
	TFExp = TFExp * TF_motif

	C = pd.read_csv('peak_gene_100k_corr.bed', header=None, sep='\t')
	d, f = ismember(C.iloc[:, 0], Element_name)
	d1, f1 = ismember(C.iloc[:, 1], geneName)
	f_2 = np.zeros(shape=(C.shape[0], 1))
	f1_2 = np.zeros(shape=(C.shape[0], 1))
	f_2[np.where(d == True)[0], 0] = f
	f1_2[np.where(d1 == True)[0], 0] = f1

	f2, ia, ic = np.unique(np.hstack((f_2[(d & d1)], f1_2[(d & d1)])), axis=0, return_index=True, return_inverse=True)
	c3 = npg.aggregate(ic, C.iloc[d & d1, 2], func='min')
	c4 = npg.aggregate(ic, C.iloc[d & d1, 3], func='min')
	c4[np.where(c4 < 0.2)] = 0
	d0 = 500000
	c = np.exp(-1 * c3 / d0) * c4
	Opn[np.isnan(Opn)] = 0
	H1 = sparse.csr_matrix((c, (f2[:, 1], f2[:, 0])), shape=(len(geneName), len(Element_name)))
	TFO = np.multiply(TF_binding.todense(), np.tile(Opn.T, (np.size(TF_binding, 0), 1)))
	H1Tdense = H1.todense().T
	BOH = np.matmul(TFO, H1Tdense)
	Score = np.multiply(np.multiply((np.dot(TFExp.reshape(-1, 1), G.reshape(-1, 1).T)), (2 ** np.abs(R2))), BOH)
	Score[np.isnan(Score)] = 0
	np.savetxt('TFTG_regulationScore.txt', Score, delimiter='\t')
	np.savetxt('TFName.txt', TFName, fmt='%s', delimiter='\n')

	TFTGControl = scio.loadmat('../../Data/TFTG_{}_nagetriveControl.mat'.format(species))
	Back_net = TFTGControl['Back_net']
	for i in range(Back_net.shape[0]):
		for j in range(Back_net.shape[1]):
			Back_net[i, j] = Back_net[i, j].item()
	d, f = ismember(Back_net[:, 0], TFName)
	d1, f1 = ismember(Back_net[:, 1], geneName)
	f_2 = np.zeros(shape=(Back_net.shape[0], 1))
	f1_2 = np.zeros(shape=(Back_net.shape[0], 1))
	f_2[np.where(d == True)[0], 0] = f
	f1_2[np.where(d1 == True)[0], 0] = f1
	f2 = np.hstack((f_2[(d & d1)], f1_2[(d & d1)]))
	Score_T_1col = Score.T.reshape(-1, 1)
	aa = (f2[:, 1]) * Score.shape[0] + f2[:, 0]
	Back_score = Score_T_1col[aa.astype(np.int64)].squeeze().T

	Cut = np.percentile(Back_score, 99)
	[b, a] = np.where(Score.T > Cut)
	c = np.where(Score_T_1col > Cut)[0]
	c1 = Score_T_1col[c]
	Net = np.column_stack((TFName[a], geneName[b]))
	a1 = np.sort(H1.T.todense(), axis=0)[::-1]
	a2 = np.argsort(H1.T.todense(), axis=0)[::-1]
	a1 = a1[:10, :]
	a2 = a2[:10, :]
	TFTG_RE = [';'.join(Element_name[np.asarray(a2[np.where((TFO[a[i], a2[:, b[i]]] > 0) & (a1[:, b[i]] > 0))[0], b[i]]).squeeze()]) for i in range(len(a))]
	for i in range(len(a)):
		kk = np.asarray(a2[np.where((TFO[a[i], a2[:, b[i]]] > 0) & (a1[:, b[i]] > 0))[0], b[i]]).squeeze()
		if kk.size == 1:
			TFTG_RE[i] = TFTG_RE[i].replace(';', '')

	d = np.sort(np.asarray(c1).squeeze())[::-1]
	f = np.argsort(np.asarray(c1).squeeze())[::-1]
	Net = np.column_stack((Net[f], np.asarray(d).squeeze(), np.asarray(TFTG_RE)[f]))
	filename = '{}_network.txt'.format(celltype)
	with open(filename, 'wt') as fid:
		fid.write('\t'.join(['TF', 'TG', 'Score', 'FDR', 'REs']) + '\n')
		for i in range(Net.shape[0]):
			fid.write('\t'.join([Net[i, 0], Net[i, 1], str(Net[i, 2]), str((np.sum(Back_score > d[i]) + 1) / len(Back_score)),Net[i, 3]]) + '\n')
