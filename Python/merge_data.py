# -*-coding:Utf-8 -*

#######################################
## merge_data.py
##
## script that merges the *_GEM.txt
## and *_samples.txt files in a given
## folder
#######################################
## Developed by:
## Sisi Chen 
## sisi.chen1@gmail.com
#######################################

#########################################
##  Load in data and combine data
#########################################

import sys
import os
from combine_matrix import *

#datapath = sys.argv[0]
#newbasename = sys.argv[1]

datapath = '/home/iamcam/Documents/RNAseqdata/Separate/combine/'
newbasename = 'ndiff_pooled'

gemfiles = [filename for filename in glob.glob(datapath + '*GEM.txt') if newbasename not in filename]
gemfiles.sort()
samplefiles = [filename for filename in glob.glob(datapath + '*samples.txt') if newbasename not in filename]
samplefiles.sort()

n = len(gemfiles)
m = len(samplefiles)

if (n!=m):
	print 'Unmatched GEMfiles or samplefiles...'
else: 
	# read in all matrices in folder: 
	mfile = gemfiles[0]
	sfile = samplefiles[0]
	m = readdata(mfile)
	s = readdata(sfile)

	for i in range(1,n):
		mfile = gemfiles[i]
		sfile = samplefiles[i]

		curr_m = readdata(mfile)
		curr_s = readdata(sfile)

		m = combine_matrix(m, curr_m)
		s.extend(curr_s)

	matrix_file = open(datapath+newbasename+'_GEM.txt', 'w+')
	for item in m:
		matrix_file.write('\t'.join([str(i) for i in item])+'\n')
	matrix_file.close()

	# make sample list into a single list
	s = [x[0] for x in s]

	sample_file = open(datapath+newbasename+'_samples.txt', 'w+')
	sample_file.write('\n'.join(s))
	sample_file.close()
	