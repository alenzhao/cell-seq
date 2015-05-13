# -*- coding: utf-8 -*-
"""
Created on Mon Apr 06 18:37:46 2015

@author: Sisi
"""
from Bio import SeqIO
import re
import sys

filename = "BPD1_S5_L001_R1_001.fastq";
output_filename = "BPD_R1_trimmed.fastq";
output_handle = open(output_filename, "w");

count = 0; #count the number of records excluded for polyT sequences

for record in SeqIO.parse(filename,"fastq"):
   
    sequence = record.seq;
    trimmedSeq = sequence[4:].lstrip('G'); #removes TAC and polyG
    trimmedStr = trimmedSeq.__str__();
    trimmedLength = len(trimmedStr);
    origLength = record.seq.__len__();
    newStart = origLength - trimmedLength;
    newRecord = record[newStart:]

    # run regular expression parsing to find all polyT and polyG strings
    polyTmatch = re.search(r'T{5,100}', trimmedStr);    
    
    if polyTmatch:
        
        count=count+1; #update count

    else:
    
        SeqIO.write(newRecord, output_handle, "fastq")

print(count)

output_handle.close()
