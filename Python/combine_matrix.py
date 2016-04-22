# -*-coding:Utf-8 -*

#######################################
## combine_matrix.py 
##
## functions that allow you to take in 
## two gene expression matrices
## (.txt format), make them consistent 
## with each other, and merge them
#######################################
## Developed by:
## Sisi Chen 
## sisi.chen1@gmail.com
#######################################
import glob
import copy

def readdata(fname):
	data = []
	with open(fname, 'r') as f:
	  for line in f:
	    line = line.strip()
	    if len(line) > 0:
	    	data.append(line.split('\t'))
	return(data)

def find(lst, a):
    return [i for i, x in enumerate(lst) if x==a]

def match_matrix(m1,m2):
	genelist1 = [x[0] for x in m1[1:]]
	genelist2 = [x[0] for x in m2[1:]]
	if genelist1==genelist2:
		return(m1,m2)
	else:
		# match the gene lists
		newgenes1 = set(genelist1).difference(set(genelist2))
		newgenes2 = set(genelist2).difference(set(genelist1))
		# add newgenes1 to m2, with all entries = 0
		for i in list(newgenes1):
			newline = [0 for x in m2[0]]
			newline[0] = i
			m2.append(newline)
		# add newgenes2 to m1, with all entries = 0
		for j in list(newgenes2):
			newline = [0 for x in m1[0]]
			newline[0] = j
			m1.append(newline)
		# reorder m2 to match m1: 
		finlist = [x[0] for x in m1[0:]]
		m2list = [x[0] for x in m2[0:]]
		m2_new = copy.deepcopy(m2)
		for k in range(1, len(m2list)):
			m2gene = m2list[k]
			fingene = finlist[k]
			if (m2gene != fingene): 
				ind = find(m2list, fingene)[0]
				m2_new[k] = m2[ind]

		return(m1, m2_new)

def combine_matrix(m1,m2):
	matched = match_matrix(m1,m2)
	m3 = matched[0]
	m4 = matched[1]
	mfin =  copy.deepcopy(m3)
	# for each element in the list, we need to combine the two lists
	for i in range(0, len(m3)):
		mfin[i]=m3[i]+m4[i][1:]
	return(mfin)