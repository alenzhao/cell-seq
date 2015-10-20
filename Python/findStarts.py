# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 16:47:14 2015

@author: sisichen


#######################################
## findStarts.py is designed to report 
## a matrix containing all the start sites 
## within a SAM file
#######################################
## Developed by:
## Sisi Chen
## sisi.chen1@gmail.com
## Summer 2015
#######################################

"""
from __future__ import division
from itertools import izip
from collections import defaultdict
from drop_params import *
import re
import math
import os
import os.path
import time
import gzip
import json


#Getting a file path for the fasta file
#starts_path = '../../RNAseqdata/DS4-IH/50bpstarts/'
#alignment_path = '../../RNAseqdata/DS4-IH/50bpalignments/'
#samfilename = 'IH50'
#startsfilename = 'IH50_refset'

starts_path = '../../RNAseqdata/DS-4/starts/'
alignment_path = '../../RNAseqdata/DS-4/alignments/'
samfilename = 'DS4_all'
startsfilename = 'DS4_refset'
ASfilter=0

gene_list = [
        'NM_008084', #GAPDH
        'NM_007393', #	ACTB
        'NM_008907',#PPIA
        'NM_013633', #POU5F1
        'NM_011443', #SOX2
        'NM_001177354', #MYC
        'NM_009556', #ZFP42
        'NM_028016' #NANOG
        
#        'NM_008899', #POU3F2
#        'NM_001165982', #DCAF17
#        'NM_001081154', #MARF1
#        'NM_001039483', #TMCO1
#        'NM_008553', #ASCL1
#        'NM_023279' #TUBB3
        ]


if os.path.isfile(starts_path+ startsfilename + '_' + 'startsmatrix.txt'):
	print "alignment_starts.txt already exists.......................................... 100 %"

else:

    gene_counter = 0;
    sam_file = gzip.open(alignment_path+samfilename+'.sam.gz','rb')
    dict_gene_starts = defaultdict(int)
    
    while True:
        line=sam_file.readline()
        if not line: 
            break
        else:
            columns = line.split("\t")
            gene = columns[2]
            startsite = columns[3]
            if (gene != '*') and (gene in gene_list):
                AS_score = int(columns[11][5:])
                if ASfilter==0: 
                    AS_score=0
                if (AS_score>-3):
                    if gene in dict_gene_starts: 
                        currlist = dict_gene_starts[gene]
                        currlist.append(startsite)
                        dict_gene_starts[gene]=currlist
                    else:
                        dict_gene_starts[gene] = [startsite]
                        gene_counter +=1
                    
                    
    sam_file.close()
    print "Creating list of start sites...\n"
    maxstarts=1
    genes=dict_gene_starts.keys()
    
    for x in genes:
        currstarts = len(dict_gene_starts[x])
        maxstarts=max(currstarts,maxstarts)

    print "maxstarts is: " 
    print maxstarts    
    
    matrix = [[0 for x in range(gene_counter)] for x in range(maxstarts)]
    for gene_key in dict_gene_starts:
        row_num = genes.index(gene_key)
        sitelist = dict_gene_starts[gene_key]
        numsites = len(sitelist)
        print numsites
        for i in range(0,numsites):
            matrix[i][row_num] = sitelist[i]

    print "Start sites tabulated"
    matrix_file = open(starts_path+ startsfilename + '_' + 'startsmatrix.txt', 'w+')
    for item in matrix:
        matrix_file.write('\t'.join([str(i) for i in item])+'\n')
    matrix_file.close()
    
    
    gene_names_file = open(starts_path + startsfilename + '_' + 'startsgenes.txt','w+')
    gene_names_file.write('\n'.join(genes))
    gene_names_file.close()
    
    
print "\n"
print "*****************************************"
print "***********End Counting Starts***********"
print "*****************************************"
print "\n"
