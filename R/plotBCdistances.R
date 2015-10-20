####################################################
##
## author: sisichen
## date: 09/29/15
##
## This script plots the pairwise distance histogram
## of a set of barcodes
##
####################################################

library(ggplot2)
library(plyr)
library(R.matlab)
library(stringdist)
library(gdata)

setwd("/home/iamcam/Documents/RNAseqdata/Barcodes/051415_BeadQC/fastqs/")
barcodes_05 = read.table("R1_pre_barcode.txt",stringsAsFactors = FALSE)
setwd("/home/iamcam/Documents/RNAseqdata/Barcodes/072815_BeadQC/fastqs/")
barcodes_07 = read.table("R1_pre_barcode.txt",stringsAsFactors = FALSE)
setwd("/home/iamcam/Documents/RNAseqdata/Barcodes/McCarroll_BeadQC/fastqs/")
barcodes_McC = read.table("R1_pre_barcode.txt",stringsAsFactors = FALSE)

barcodedf_05 = data.frame(table(barcodes_05[,1]))
barcodedf_07 = data.frame(table(barcodes_07[,1]))
barcodedf_McC = data.frame(table(barcodes_McC[,1]))
colnames(barcodedf_05) <- c('barcode', 'Freq')
colnames(barcodedf_07) <- c('barcode', 'Freq')
colnames(barcodedf_McC) <- c('barcode', 'Freq')

result<-stringdistmatrix(barcodes_05[1:20000,1], barcodes_05[1:20000,1],method="hamming")
resultlist<-upperTriangle(result)
distdf_05<-data.frame(table(resultlist))

result<-stringdistmatrix(barcodes_07[1:20000,1], barcodes_07[1:20000,1],method="hamming")
resultlist<-upperTriangle(result)
distdf_07<-data.frame(table(resultlist))

result<-stringdistmatrix(barcodes_McC[1:20000,1], barcodes_McC[1:20000,1],method="hamming")
resultlist<-upperTriangle(result)
distdf_McC<-data.frame(table(resultlist))

randBC<-function(x){
  BC<-paste(sample(c('A','T','C','G'),12,replace=TRUE), collapse="")
  return(BC)
}
zeroArray<-rep(0,20000)
randomBCs<-mapply(randBC, zeroArray)

randresult<-stringdistmatrix(randomBCs, randomBCs,method="hamming")
randresultlist<-upperTriangle(randresult)
randdistdf<-data.frame(table(randresultlist))

##############################

u_sample<-sample(bc_unique,20000, replace=TRUE)
u_result<-stringdistmatrix(u_sample, u_sample,method="hamming")
u_resultlist<-upperTriangle(u_result)
u_distdf<-data.frame(table(u_resultlist))
