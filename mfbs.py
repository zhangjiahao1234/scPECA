import numpy as np
import pandas as pd
from ismember import ismember
import scipy.sparse as sparse
import os
import scipy.io as scio
from numpy.matlib import repmat
import numpy_groupies as npg

def mfbs(MCFile, TFName, Element_name, motifName, motifWeight, Match2):
	MC = pd.read_csv(MCFile, sep='\t', header=None)
	Mf3 = MC.iloc[:,2]
	Md1, Mf1 = ismember(MC.iloc[:,0], Element_name)
	Md2, Mf2 = ismember(MC.iloc[:,1], pd.DataFrame(motifName))
	Mf1_2 = np.zeros(shape=(MC.shape[0],1))
	Mf2_2 = np.zeros(shape=(MC.shape[0],1))
	Mf1_2[np.where(Md1 == True)[0], 0] = Mf1
	Mf2_2[np.where(Md2 == True)[0], 0] = Mf2
	Mt1 = np.setdiff1d(np.arange(len(motifName)), np.unique(Mf2))
	Mf2 = np.concatenate((np.squeeze(Mf2_2[(Md1 & Md2)]), Mt1 ))
	Mf1 = np.concatenate((np.squeeze(Mf1_2[(Md1 & Md2)]), np.ones(len(Mt1))))
	Mf3 = np.concatenate((Mf3[np.where(Md1 & Md2)[0]], np.zeros(len(Mt1))))
	Mt1 = np.setdiff1d(np.arange(len(Element_name)), np.unique(Mf1))
	Mf1 = np.concatenate((Mf1, Mt1))
	Mf2 = np.concatenate((Mf2, np.ones(len(Mt1))))
	Mf3 = np.concatenate((Mf3, np.zeros(len(Mt1))))
	MMotif_binding = sparse.csr_matrix((Mf3, (Mf2, Mf1)), shape=(len(motifName), len(Element_name)))
	MMotif_binding = sparse.diags(np.squeeze(np.asarray(1/(motifWeight + 0.1)))) * MMotif_binding
	MMotif_binding = np.log(MMotif_binding.todense()+1)
	MTF_binding = np.zeros((len(TFName), len(Element_name)))
	Mf1_2 = np.zeros(shape=(Match2.shape[0], 1))-1
	Mf2_2 = np.zeros(shape=(Match2.shape[0], 1))-1
	Md1, Mf1 = ismember(pd.DataFrame(Match2[:, 0]), motifName)
	Md2, Mf2 = ismember(pd.DataFrame(Match2[:, 1]), TFName)
	Mf1_2[np.where(Md1 == True)[0], 0] = Mf1
	Mf2_2[np.where(Md2 == True)[0], 0] = Mf2
	Mf1 = np.squeeze(Mf1_2)
	Mf2 = np.squeeze(Mf2_2)
	for i in range(len(TFName)):
		Ma = np.where(Mf2 == i)[0]
		if len(Ma) > 1:
			MTF_binding[i, :] = MMotif_binding[Mf1[Ma].astype(np.int64), :].max(axis=0)
		elif len(Ma) == 1:
			MTF_binding[i, :] = MMotif_binding[Mf1[Ma].astype(np.int64), :]
		else:
			MTF_binding[i, :] = np.zeros((len(Element_name),))
	MTF_binding = sparse.csr_matrix(MTF_binding)
	return MTF_binding

def mfbs_c(N,TFName, Element_name, motifName, motifWeight, Match2, celltype):
	TFB = sparse.csr_matrix(([], ([], [])), shape=(len(TFName), len(Element_name)))
	for i in range(N):
		TFB1 = mfbs("./{}/.MotifTarget".format(celltype)+str(i+1)+".txt",TFName, Element_name, motifName, motifWeight, Match2)
		TFB = TFB + TFB1
	return TFB




