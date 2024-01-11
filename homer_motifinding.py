import numpy as np
import subprocess
from multiprocessing import Process
import os

def Split(N):
	f = open("region.txt")
	pn = f.readlines();f.close()
	NN = len(pn);NM = int(NN/N)+1
	for i in range(N):
		g = open(".region"+str(i+1)+'.txt','w')
		for j in range(len(pn)):
			if int(j/NM) == i:
				g.write(pn[j])
		g.close()

def Homer_prior(peak,genome,output):
	command="mkdir "+output+";findMotifsGenome.pl "+peak+" "+genome+" ./"+output+" -size given -find ../../Data/all_motif_rmdup -preparsedDir ../../Homer/ > ./"+output+"/"+output+".bed;cat ./"+output+"/"+output+".bed | awk 'NR>1'|cut -f 1,4,6 > "+output+".txt;rm -rf "+output
	ph=subprocess.Popen(command,shell=True)
	return_code=ph.wait()

def Homer(peak,genome,output):
	command = "mkdir "+output+";findMotifsGenome.pl " + peak + " " + genome + " ./" + output + " -size given -find ../../Data/all_motif_rmdup > ./" + output + "/" + output + ".bed;cat ./" + output + "/" + output + ".bed | awk 'NR>1'|cut -f 1,4,6 > " + output + ".txt;rm -rf " + output
	ph = subprocess.Popen(command, shell=True)
	return_code = ph.wait()

def MotifFind(genome, num_processes, prior):
	"""
	Args:
		genome (str): 
		num_processes (int): 线程数
		prior (int): 是否有已经扫好的homer motif文件
	"""
	processes = []
	Split(num_processes)
	if prior == 1:
		for p in range(num_processes):
			process = Process(target=Homer_prior, args=(".region"+str(p+1)+'.txt',str(genome),".MotifTarget"+str(p+1)))
			processes.append(process)
	else:
		for p in range(num_processes):
			process = Process(target=Homer, args=(".region"+str(p+1)+'.txt',str(genome),".MotifTarget"+str(p+1)))
			processes.append(process)
	for process in processes:
		process.start()
	for process in processes:
		process.join()

