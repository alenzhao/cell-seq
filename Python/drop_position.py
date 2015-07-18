# -*-coding:Utf-8 -*

#######################################
## drop.py is designed to preprocess
## the fastq files from a plate exper 
## iment, and align data to a referen
## ce genome.
#######################################
## Developed by:
## Paul Rivaud
## paulrivaud.info@gmail.com
## Summer 2015
#######################################

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

def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)

def round_figures(x, n):
	return round(x, int(n - math.ceil(math.log10(abs(x)))))

#Getting a files list from the directory
files = os.listdir(dir_path_fastqs)
#Cleaning list of files in the directory(removing *noTA.fastq and *umi.txt files):
files_noNoTA = [f for f in files if 'noTA' not in f]
files_noNoTA_noUMI = [f for f in files_noNoTA if 'umi' not in f]
#Removing files that do not have a .fastq file extension:
fastq_files = [f for f in files_noNoTA_noUMI if '.fastq.gz' in f]
#Sorting the list alphabetically in order to get R1 and R2 pairs) 
fastq_files.sort()
############################### PREPROCESSING ###############################
print "\n"
print "**********************************"
print "**    Starting preprocessing    **"
print "**********************************"
print "\n"

nb_fastqs = len(fastq_files)
curr_fastq = 2

for f1, f2 in grouped(fastq_files, 2):
	percent = int((curr_fastq/nb_fastqs)*100)
	curr_fastq+=2
	if os.path.isfile(dir_path_fastqs+f1.split(os.extsep)[0]+'_noTA.fastq.gz'):
		print f1
		print f2
		print "\tFiles already preprocessed. Moving on to other files...",percent,"%"
	else:
		print "Begin..."
		print "Opening fastq files..."
		fastq1_path = dir_path_fastqs+f1.split(os.extsep)[0]
		fastq2_path = dir_path_fastqs+f2.split(os.extsep)[0]
		#file1_fastq = open(fastq1_path+'.fastq','r')
		#file2_fastq = open(fastq2_path+'.fastq','r')
		file1_fastq = gzip.open(fastq1_path+'.fastq.gz','rb')
		file2_fastq = gzip.open(fastq2_path+'.fastq.gz','rb')
		#Dictionary containing (seq)-(umis list) pairs
		seq_dictionary = defaultdict(list)
		print f1
		print f2
		print "\tReading files..."
		print "\tWriting files..."
		#file_noTA = open(fastq1_path+'_noTA.fastq', 'w+')
		file_noTA = gzip.open(fastq1_path+'_noTA.fastq.gz', 'wb')
		#file_umi = open(fastq1_path+'_umi.txt', 'w+')
		file_barcode = open(fastq1_path+'_barcode.txt', 'w+')
		#Stats about the trimming process
		total_reads = 0
		#complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
		while True:
			f1_line1 = file1_fastq.readline()
			f1_line2 = file1_fastq.readline()
			f1_line3 = file1_fastq.readline()
			f1_line4 = file1_fastq.readline()
			f2_line1 = file2_fastq.readline()
			#f2_line2 = "".join(complement.get(base, base) for base in reversed(file2_fastq.readline()))
			f2_line2 = file2_fastq.readline()
			f2_line3 = file2_fastq.readline()
			#f2_line4 = "".join(reversed(file2_fastq.readline()))
			f2_line4 = file2_fastq.readline()
			if not f1_line1:
				break
			else:
				total_reads+=1
			if f1_line2[:3] == 'TAC' and f1_line2[:6] != 'TACGGG':
				if tso not in f2_line2:
					barcode = f1_line2[tac_length:tac_length+barcode_length]
					umi = f1_line2[tac_length+barcode_length:tac_length+barcode_length+umi_length]
					#Checking trimmed sequence length
					if f2_line2 in seq_dictionary:
						if umi not in seq_dictionary[f2_line2]:
							seq_dictionary[f2_line2].append(umi)
							file_barcode.write(barcode+'\n')
							#file_umi.write(umi+'\n')
							file_noTA.write(f1_line1)
							file_noTA.write(f2_line2)
							file_noTA.write(f1_line3)
							file_noTA.write(f2_line4)
					else:
						seq_dictionary[f2_line2].append(umi)
						file_barcode.write(barcode+'\n')
						#file_umi.write(umi+'\n')
						file_noTA.write(f1_line1)
						file_noTA.write(f2_line2)
						file_noTA.write(f1_line3)
						file_noTA.write(f2_line4)
		#file_umi.close()
		file_barcode.close()
		file1_fastq.close()
		file2_fastq.close()
		file_noTA.close()
		print "\tTotal reads: ",total_reads
		print "...................................................................",percent,"%"


############################### ALIGNMENT ###############################
start_time = time.time()
#Listing *_noTA.fastq files:
preprocessed_files = os.listdir(dir_path_fastqs)
#Cleaning list of files in the directory (Only *noTA.fastq files):
files_noTA = [f for f in preprocessed_files if 'noTA.fastq.gz' in f]
files_noTA.sort()
bowtie_opt = ' '.join(bowtie2_options)
print "\n"
print "**********************************"
print "**      Starting alignment      **"
print "**********************************"
print "\n"

nb_noTA = len(files_noTA)
curr_noTA = 1

for file_noTA in files_noTA:
	percent = int((curr_noTA/nb_noTA)*100)
	curr_noTA+=1
	if os.path.isfile(dir_path_alignment+file_noTA.split(os.extsep)[0]+'.sam.gz'):
		print "Skipping alignment: "+file_noTA.split(os.extsep)[0]+'.sam.gz already exists...',percent,"%"
	else:
		print "Starting alignment with bowtie2 for", file_noTA
		os.system("nice "+bowtie2_dir+" bowtie2 "+bowtie_opt+\
			" -x "+reference_genome+" -U "+dir_path_fastqs+file_noTA+\
			" -S "+dir_path_alignment+file_noTA.split(os.extsep)[0]+'.sam')
		print file_noTA+" aligned..."
		print file_noTA.split(os.extsep)[0]+'.sam created...'
		print "...................................................................",percent,"%"
		os.system("gzip "+dir_path_alignment+file_noTA.split(os.extsep)[0]+'.sam')
		print "Sam file compressed..."
		print "...................................................................",percent,"%"
total_time = time.time() - start_time
print "Reads alignment time:"
print int(total_time/60),"min",int(total_time%60),"sec"

############################### Data gathering ###############################
print "\n"
print "**********************************"
print "**       Retrieving data        **"
print "**********************************"
print "\n"

sum_path = os.path.commonprefix([dir_path_fastqs, dir_path_alignment])
if os.path.isfile(sum_path+'matrix.txt'):
	print "matrix.txt already exists.......................................... 100 %"
else:
	gene_counter = 1
	barcode_counter = 1
	sam_file = gzip.open(dir_path_alignment+file_noTA.split(os.extsep)[0]+'.sam.gz','rb')
	#umi_file = open(fastq1_path+'_umi.txt','r')
	list_f = os.listdir(dir_path_fastqs)
	barcode_list = [f for f in list_f if '_barcode.txt' in f]
	barcode_file = open(dir_path_fastqs+barcode_list[0],'r')
	dict_genes_barcode = defaultdict(dict)
	dict_gene_counter = defaultdict(int)
	dict_barcode_counter = defaultdict(int)
	dict_barcode_occurences = defaultdict(str)
	print "Creating barcode occurence dictionary..."
	while True:
		barcode = barcode_file.readline()
		barcode = barcode.replace('\n','')
		if not barcode:
			break
		else:
			if barcode not in dict_barcode_occurences:
				dict_barcode_occurences[barcode] = 1
			else:
				dict_barcode_occurences[barcode] += 1
	print "Trimming barcode occurence dictionary..."
	for barc in dict_barcode_occurences.keys():
		if dict_barcode_occurences[barc] < occ_threshold:
			del dict_barcode_occurences[barc]
	print "Barcode occurence dictionary trimmed..."
	barcode_file.close()
	barcode_file2 = open(dir_path_fastqs+barcode_list[0],'r')
	print "Storing data in dictionaries..."
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
print "Ready for PCA..."
print "\n"
print "**********************************"
print "***********Pipeline end***********"
print "**********************************"
print "\n"