# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 16:47:14 2015

@author: sisi


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


sum_path = os.path.commonprefix([dir_path_fastqs, dir_path_alignment])
samfilename = '1_TAGATCGC_L001_R1_001_noTA';

if os.path.isfile(sum_path+'alignment_starts.txt'):
	print "alignment_starts.txt already exists.......................................... 100 %"
else:

    gene_counter = 1;
    sam_file = gzip.open(dir_path_alignment+samfilename+'.sam.gz','rb')
    dict_gene_starts = defaultdict(int)
    
    while True:
        line=sam_file.readline()
        if not line: 
            break
        else:
            columns = line.split("\t")
            gene = columns[2]
            startsite = columns[3]
            if gene != '*':
                if gene in dict_gene_starts: 
                    currlist = dict_gene_starts[gene];
                    currlist.append(startsite);
                    dict_gene_starts=currlist;
                else:
                    dict_gene_counter[gene] = [startsite];
                    gene_counter +=1;
                    
                    
    sam_file.close()
    print "Creating list of start sites...\n"
    maxstarts=1;
    genes=dict_genes.keys();
    
    for x in genes:
        currstarts = size(dict_genes[x]);
        maxstarts=max(currstarts,maxstarts);
    
    matrix = [[0 for x in range(gene_counter)] for x in range(maxstarts)]
    for gene_key in dict_genes:
        col_num = genes.index(gene_key);
        sitelist = dict_gene_starts[gene_key]
        for startsite in sitelist:
            row_num = sitelist.index(startsite);
            matrix[rownum]
        

### everything is left off here. need to continue sorting the code above.         
        
        
     while True:
		line=sam_file.readline()
		if not line:
			break
		else:
			columns = line.split("\t")
			gene = columns[2]
			barcode = barcode_file2.readline()
			barcode = barcode.replace('\n','')
			#If read aligned, columns[2] is different from '*'
			#print gene, barcode
			if gene != '*' and barcode in dict_barcode_occurences:
				if gene not in dict_gene_counter:
					dict_gene_counter[gene] = gene_counter
					gene_counter+=1
				if barcode not in dict_barcode_counter:
					dict_barcode_counter[barcode] = barcode_counter
					barcode_counter+=1
				if gene in dict_genes_barcode:
					if barcode in dict_genes_barcode[gene].keys():
						dict_genes_barcode[gene][barcode] +=1
					else:
						dict_genes_barcode[gene][barcode] = 1
				else:
					dict_genes_barcode[gene] = {barcode : 1}
	barcode_file2.close()
	sam_file.close()
	print "Data stored in dictionaries........................................",percent,"%"
	print "Creating genes-cells matrix...\n"
	print gene_counter-1, "genes"
	print barcode_counter-1, "cells"
	matrix = [[0 for x in range(barcode_counter)] for x in range(gene_counter)]
	for key_barcode in dict_barcode_counter:
		col_num = dict_barcode_counter[key_barcode]
		matrix[0][col_num] = key_barcode
	for key_gene in dict_genes_barcode:
		row_num = dict_gene_counter[key_gene]
		matrix[row_num][0] = key_gene
		for key_barcode in dict_genes_barcode[key_gene]:
			col_num = dict_barcode_counter[key_barcode]
			matrix[row_num][col_num]=dict_genes_barcode[key_gene][key_barcode]
	print "Genes-cells matrix created.........................................",percent,"%"
	matrix_file = open(sum_path+'matrix.txt', 'w+')
	for item in matrix:
		matrix_file.write('\t'.join([str(i) for i in item])+'\n')
	matrix_file.close()
print "\n"
print "**********************************"
print "***********End Counting Starts***********"
print "**********************************"
print "\n"