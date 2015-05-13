# -*- coding: utf-8 -*-
"""
Created on Mon Apr 06 18:37:46 2015

argv[1]: directory with fastq files
only for paired end reads right now

@author: Sisi
"""
from Bio import SeqIO
import itertools
import re
import sys
import os
import csv

currPath = sys.argv[1];

fileList = os.listdir(currPath);
#find fastq files only
fastqList = [s for s in fileList if "001.fastq" in s]

#find number of fastq files
nSamples = len(fastqList)/2;
countList = [0]*nSamples;
allreadcountList=[0]*nSamples;

for i in range(1,nSamples+1):

    currFileList = [s for s in fastqList if "_S"+str(i)+"_" in s]
    filename1 = [s for s in currFileList if "R1" in s][0]
    filename2 =  [s for s in currFileList if "R2" in s][0]
    
    fastq_iter1 = SeqIO.parse(open(os.path.join(currPath,filename1)),"fastq");
    fastq_iter2 = SeqIO.parse(open(os.path.join(currPath,filename2)),"fastq");
    
    output1_fn = os.path.join(currPath, filename1[0:-9]+"trimmed.fastq");    
    output2_fn = os.path.join(currPath, filename2[0:-9]+"trimmed.fastq");
    
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
        newRec1 = rec1[newStart1:];
    
        # run regular expression parsing to find all polyT and polyG strings
        polyTmatch = re.match(r'T{5,100}', trimmedStr1);    
        
        if polyTmatch:
            
            count=count+1; #update count
    
        else:
        
            SeqIO.write(newRec1, output1_handle, "fastq")
            SeqIO.write(rec2, output2_handle, "fastq")
    
    print("records removed", count)
    print('out of ', allreadcount)
    countList[i]=count;
    allreadcountList[i]=allreadcount;
    
    output1_handle.close()
    output2_handle.close()
    

#open csv file to write records removed data to: 
#with open(currPath+"\\records.csv", 'rb') as csvfile:
    
csvfilename=os.path.join(currPath,"records.csv");    
    
with open(csvfilename, 'wb') as csvfile:
    fieldnames=['sample','removed', 'total'];
    writer=csv.writer(csvfile,fieldnames);
    writer.writerow(fieldnames);
    
    for i in range(0,nSamples):
        writer.writerow([i, countList[i], allreadcountList[i]]);
    