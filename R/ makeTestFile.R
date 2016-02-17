####################################################
##
## author: sisichen
## date: 10/29/15
##
## This script is for generating a large test file
## to test the pipeline
##
####################################################

library(stringi)
library(dplyr)
library(gdata)
library(gtools)
library(tidyr)

###############################
##   generate main barcodes 
###############################

setwd('~/Documents/cell-seq/R/testfiles/')

numBCs<-20
fractRand <- 0.5
pmrate<-1e-3
minreads<-10000
meanreads<-50000

mainBCs <- stri_rand_strings(numBCs, 12, pattern='[ATCG]')
mainNum <- round(meanreads*(rexp(numBCs,rate=1)))
mainNum[which(mainNum<minreads)] <- minreads # set the minimum number of reads
fragNum <-rpois(numBCs,1.5)

BCfrag<-function(BC,num) {
  
  if (num!=0) {
    baseBC<-substring(BC,1,nchar(BC)-num)
    x<-data.frame(permutations(4, num, c('A','T','C','G'), repeats.allowed=TRUE))
    suffixes<-unlist(unite(x,'new', 1:num, sep=''))
    newBCs<-paste(baseBC,suffixes,sep='')
    return(unname(newBCs))
  } else {
    return(BC)
  }
  
}

numfrag <- function(readnum, num){
  newreads<-rep(round(readnum/(4^num)), 4^num)
  return(newreads)
}

newBCs<-unname(mapply(BCfrag,mainBCs,fragNum))
newBCs<-unlist(newBCs)
newNum<-mapply(numfrag,mainNum,fragNum)
newNum<-unlist(newNum)

## save list of barcodes, and readnums/barcode as text files
write.table(cbind(newBCs, newNum),file='test4_BCs_afterfrag.txt',sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(cbind(mainBCs, mainNum),file='test4_BCs.txt',sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)

cellBCs <- unname(unlist(mapply(rep,newBCs,newNum)))
numExtra<-sum(newNum)

## make random mutations in the newBCs (sequencing error)

mutator<-function(BC){
  chance<-rbinom(nchar(BC),1,pmrate)
  inds<-which(chance>0)
  n<-length(inds)
  
  if (n>0) {
    for (i in 1:n) {
      char<-substring(BC,inds[i],inds[i])
      newset<-setdiff(c('A','T','C','G'),char)
      newchar<-sample(newset,1)
      BC<-paste(substring(BC,1,inds[i]-1),newchar,substring(BC,inds[i]+1, nchar(BC)), sep='')
    }
  } 
  return(BC)
}

mutBCs<-sapply(cellBCs,mutator)
extraBCs <- stri_rand_strings(numExtra, 12, pattern='[ATCG]')
allBCs <- c(mutBCs,extraBCs) 
numReads <- length(allBCs)

## generate R2 sequences
seq1='AGGCTGCAGTACAAGGATGATGCCCTTATTCAGGAGAGGCTGGAATATGA' #xkr4
seq2='AAATATGGTCAGGACAAGAAATATATCCTTAAAATATTGTCTAGTACATA' #RP1
seq3='CACAAACAGCGGGGTGTTCAAATGTGATATAATTGTTCTGAGAAAAATTC' #SOX17
seq4='CCTTAAAAAGCAATAAGTCTTGAGTGTAGGTCTGTCCTTCCTGGAATGTG' #MRPL15
seq5='CAACATCTGGGCTTAAAGTCAGGGCAAAGCCAGGTTCCTTCCTTCTTCCA' #NANOG
seq6='TATGTGAAGTTTAGAAGCCTCAAGCTGTGAGGCCCAGGGCTGAGGAATAA' #DPPA3
seq7='GTATGACAATGAATACGGCTACAGCAACAGGGTGGTGGACCTCATGGCCT' #GAPDH

seqs<-c(seq3,seq4,seq5,seq6,seq7)

equiSeq<-function(num, seqlist) {
  num_each<-floor(num/length(seqlist))
  newlists<-lapply(seqlist,rep,num_each)
  newlist<-unlist(newlists)
  
  diff <- num - length(newlist)
  if (diff!=0) {
   newlist<-c(newlist,rep(tail(seqlist,1), diff))
  }
  return(newlist)
}

allR2s<-unlist(lapply(newNum, equiSeq, seqs))
extraR2s<-stri_rand_strings(numExtra,50, pattern='[ATCG]')
allR2s<-c(allR2s, extraR2s)

## generate UMIs
UMIs <- stri_rand_strings(length(seqs), 8, pattern='[ATCG]')
allUMIs <- unlist(lapply(newNum, equiSeq, UMIs))
extraUMIs<-stri_rand_strings(numExtra, 50, pattern='[ATCG]')
allUMIs<-c(allUMIs, extraUMIs)

## concatenate strings 
TACs<-rep('TAC',numReads) 
allR1s<-paste(TACs,allBCs,allUMIs,sep='')

## Convert all sequences to fastq format: 
## need: 
#  @M02038:45:000000000-AJ5AB:1:1102:7054:12812 1:N:0:1
#  TACGGGATAGCGAGACATTTCTTNN
#  +
#  CCCCCGGGGGGGGGGGGGGGGGGGG

n <- length(allR1s)
quallength1<-nchar(allR1s[1])
quallength2<-nchar(allR2s[1])
ats<-rep('@M',n)
R1quals<-stri_rand_strings(n,quallength1,pattern='[CG]')
R2quals<-stri_rand_strings(n,quallength2,pattern='[CG]')
pluses<-rep('+',n)

R1all<-rbind(ats,allR1s,pluses,R1quals)
R1all<-c(R1all)

R2all<-rbind(ats,allR2s,pluses,R2quals)
R2all<-c(R2all)

write.table(R1all,file='test4_R1.fastq',sep='\n', quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(R2all,file='test4_R2.fastq',sep='\n', quote=FALSE, col.names=FALSE, row.names=FALSE)


