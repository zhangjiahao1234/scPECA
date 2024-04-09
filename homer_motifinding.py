import numpy as np
import subprocess
from multiprocessing import Process
import os

def Split(N, celltype):
	f = open("./{}/region.txt".format(celltype))
	pn = f.readlines();f.close()
	NN = len(pn);NM = int(NN/N)+1
	for i in range(N):
		g = open("./{}/.region".format(celltype)+str(i+1)+'.txt','w')
		for j in range(len(pn)):
			if int(j/NM) == i:
				g.write(pn[j])
		g.close()

def Homer_prior(peak,genome,celltype, output,pkg_path=None, prior_path=None):
	command="mkdir "+ './{}/'.format(celltype)+output+";findMotifsGenome.pl "+peak+" "+genome+' ./{}/'.format(celltype)+output+f" -size given -find {os.path.join(pkg_path, 'Data/all_motif_rmdup')} -preparsedDir {prior_path} > ./{celltype}/"+output+"/"+output+".bed;cat ./{}/".format(celltype)+output+"/"+output+".bed | awk 'NR>1'|cut -f 1,4,6 > "+'./{}/'.format(celltype)+output+".txt;rm -rf "+'./{}/'.format(celltype)+output
	logfile = open(output + ".log", 'w')
	ph = subprocess.Popen(command, shell=True, stderr=logfile)
	return_code=ph.wait()

def Homer(peak,genome,celltype,output, pkg_path=None):
	command = "mkdir "+'./{}/'.format(celltype)+output+";findMotifsGenome.pl " + peak + " " + genome + ' ./{}/'.format(celltype) + output + " -size given -find {} > ".format(os.path.join(pkg_path, 'Data/all_motif_rmdup')) + './{}/'.format(celltype) + output + "/" + output + ".bed;cat ./{}/".format(celltype) + output + "/" + output + ".bed | awk 'NR>1'|cut -f 1,4,6 > " + './{}/'.format(celltype)+output + ".txt;rm -rf " + './{}/'.format(celltype)+output
	logfile = open('./{}/'.format(celltype)+output + ".log", 'w')
	ph = subprocess.Popen(command, shell=True, stderr=logfile)
	return_code = ph.wait()

def MotifFind(genome, celltype, pkg_path, num_processes=int(os.cpu_count()*0.75), prior=0, prior_path=None):
	"""
	Args:
		genome (str): 
		num_processes (int): 线程数
		prior (int): 是否有已经扫好的homer motif文件
		celltype (str)：
	"""
	processes = []
	Split(num_processes,celltype)
	if prior == 1:
		for p in range(num_processes):
			process = Process(target=Homer_prior, args=("./{}/.region".format(celltype)+str(p+1)+'.txt', str(genome), celltype, ".MotifTarget".format(celltype)+str(p+1), prior_path, pkg_path))
			processes.append(process)
	else:
		for p in range(num_processes):
			process = Process(target=Homer, args=("./{}/.region".format(celltype)+str(p+1)+'.txt', str(genome), celltype, ".MotifTarget".format(celltype)+str(p+1), pkg_path))
			processes.append(process)
	for process in processes:
		process.start()
	for process in processes:
		process.join()

