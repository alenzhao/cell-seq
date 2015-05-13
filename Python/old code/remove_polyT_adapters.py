# -*- coding: utf-8 -*-
"""
Created on Mon Apr 06 18:37:46 2015

@author: Sisi
"""
from Bio import SeqIO
import itertools
import re
import sys

filename1 = sys.argv[1];
filename2 = sys.argv[2];

fastq_iter1 = SeqIO.parse(open(filename1),"fastq")
fastq_iter2 = SeqIO.parse(open(filename2),"fastq")

output1_fn = filename1[0:-6]+"_trimmed.fastq";
output2_fn = filename2[0:-6]+"_trimmed.fastq";

output1_handle = open(output1_fn, "w")
output2_handle = open(output2_fn, "w")

allreadcount = 0;
count = 0; #count the number of records excluded for polyT sequences

for rec1, rec2 in itertools.izip(fastq_iter1, fastq_iter2):
   
    allreadcount = allreadcount + 1;
    sequence1 = rec1.seq;
    
    trimmedSeq1 = sequence1[4:].lstrip('G'); #removes TAC, polyG 
    trimmedStr1 = trimmedSeq1.__str__();
    
    trimmedLen1 = len(trimmedStr1);
    origLength = rec1.seq.__len__();
    newStart1 = origLength - trimmedLen1;
    newRec1 = rec1[newStart1:]

    # run regular expression parsing to find all polyT and polyG strings
    polyTmatch = re.match(r'T{5,100}', trimmedStr1);    
    
    if polyTmatch:
        
        count=count+1; #update count

    else:
    
        SeqIO.write(newRec1, output1_handle, "fastq")
        SeqIO.write(rec2, output2_handle, "fastq")

print("records removed", count)
print('out of ', allreadcount)

output1_handle.close()
output2_handle.close()
