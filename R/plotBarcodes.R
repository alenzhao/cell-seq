####################################################
##
## author: sisichen
## date: 07/22/15
##
## This script plots a histogram of the barcodes distribution
## given the list of barcodes from the pipeline, as well as
## the cumsum of barcodes ordered by barcode 
##
####################################################

library(ggplot2)
library(plyr)
library(R.matlab)

setwd("~/Documents/RNAseqdata/DS-3")

barcodes = read.table("50-cells_S1_L001_R1_001_barcode.txt",stringsAsFactors = FALSE)
#barcodes = read.table("50-cells-barcode-10bp.txt",stringsAsFactors = FALSE)
 
barcodedf = data.frame(table(barcodes[,1]))
colnames(barcodedf) <- c('barcode', 'Freq')
bc_50reads = barcodedf[which(barcodedf$Freq>50),]

barcodeplot<-ggplot(barcodedf[which(barcodedf$Freq>1),],aes(Freq)) + 
  geom_histogram() + 
  labs('DS3 barcode frequences, xlin-ylin');
plot(barcodeplot);

barcodeplot<-ggplot(barcodedf[which(barcodedf$Freq>1),],aes(Freq)) + 
  geom_histogram() + 
  scale_x_log10() + 
  labs('DS3 barcode frequences, xlog-ylin');
plot(barcodeplot);

barcodeplot<-ggplot(barcodedf,aes(Freq)) + 
  geom_histogram() + 
  scale_x_log10() + scale_y_log10() + 
  labs('DS3 barcode frequences, xlog-ylog');
plot(barcodeplot);