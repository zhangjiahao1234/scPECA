import numpy as np
import os
import scipy.io as scio
from collections import Counter
from ismember import ismember
import numpy_groupies as npg

def mf_collect(genome):
	if "hg" in genome:
		species = "human"
	elif "mm" in genome:
		species = "mouse"
	filename = 'knownResults.txt'
	fileID = open(filename)
	C = np.loadtxt(fileID, dtype=str, skiprows=1)
	fileID.close()

	for i in range(C.shape[0]):
		C[i, 6] = float(C[i, 6].strip('%'))
		C[i, 8] = float(C[i, 8].strip('%'))
	Score = np.zeros(shape=(C.shape[0], 3))
	Score[:, 0] = -np.log10(C[:, 2].astype(float))
	Score[:, 1] = (C[:, 6].astype(float) + 0.1) / (C[:, 8].astype(float) + 0.1)
	Score[Score[:, 0] > 100, 0] = 100
	Score[:, 2] = np.sqrt(Score[:, 0] * Score[:, 1])

	f = np.argsort(Score[:, 2])[::-1]
	Score = Score[f, :]
	Name = C[f, 0]

	filename = 'knownResults_rank.txt'
	fid = open(filename, 'w')
	fid.write('Motif\t-log10(p)\tFoldChange\tScore\n')
	for i in range(Name.shape[0]):
		fid.write('{}\t{}\t{}\t{}\n'.format(Name[i], Score[i, 0], Score[i, 1], Score[i, 2]))
	fid.close()

	Motifmatch = scio.loadmat('../../../Data/MotifMatch_{}_rmdup.mat'.format(species))

	Match2 = np.empty((Motifmatch['Match2'].shape[0], 2), dtype='U100')
	for i in range(Motifmatch['Match2'].shape[0]):
		Match2[i, 0] = Motifmatch['Match2'][i, 0].item()
		Match2[i, 1] = Motifmatch['Match2'][i, 1].item()

	d, f = ismember(Match2[:, 0], Name)
	Score1 = Score[f, :]
	Match2 = Match2[np.where(d == True), :]
	Match2 = np.squeeze(Match2)
	TF, ic = np.unique(Match2[:, 1], return_inverse=True)
	TFScore = npg.aggregate(ic, Score1[:, 2], func='max')
	f = np.argsort(TFScore)[::-1]
	Results = np.column_stack((TF[f], TFScore[f]))
	filename = 'knownResults_TFrank.txt'
	fid = open(filename, 'w')
	for i in range(Results.shape[0]):
		fid.write('{}\t{}\n'.format(Results[i, 0], Results[i, 1]))
	fid.close()
