####################################################
##
## author: sisichen
## date: 07/27/15
##
## This script plots the distribution of reads, and 
## number of genes expressed for a mixed species
## experiment
##
####################################################

library(ggplot2)
library(plyr)
library(R.matlab)
library(grep)
library(reshape2)
library(MASS)

setwd("~/Documents/RNAseqdata/DS-3")

s1df<-read.table('sample2_truncBarcode_matrix_stats.txt', header=TRUE, nrows=44097)
s1df_summ<-read.table('sample2_truncBarcode_matrix_stats.txt', skip=44098, nrows=3)
#s1df<-read.table('sample2_truncBarcode_matrix_stats.txt', header=TRUE, nrows=46225)
#s1df_summ<-read.table('sample2_truncBarcode_matrix_stats.txt', skip=46226, nrows=3)
barcodelist<-colnames(s1df)
gene_info<-s1df[,1:2];
numcells<-dim(s1df)[2] - 2;

# transpose data matrix:
m<-s1df[,2:(numcells+2)]
t_matrix<-t(m)
new_s1df<-data.frame(t_matrix)
colnames(new_s1df)<-gene_info$X0
new_s1df<-tail(new_s1df,-1)

## make s1df_summ variable

s1df_summ$V1=NULL
s1df_summ$V3=NULL
row.names(s1df_summ)<-s1df_summ$V2
s1df_summ$V2=NULL
ishuman<-round(s1df_summ[1,]/s1df_summ[3,])
rownames(ishuman)<-c('ishuman')
s1df_summ<-rbind(s1df_summ,ishuman)
colnames(s1df_summ)<-rownames(new_s1df)

# add ishuman column to s1df
new_s1df<-cbind(t(s1df_summ["ishuman",]),new_s1df)

# make numreads/numgenes summary data frame
y_dim=dim(new_s1df)[2];
numreads<-rowSums(data.matrix(new_s1df[,2:y_dim]))
numgenes<-apply(new_s1df,1,function(c)sum(c!=0))

meangenes<-apply(new_s1df,1,function(c)mean(c))

# plot num reads separated by human vs mouse
s1df_a<-data.frame(numreads, numgenes)
s1df_a<-cbind(factor(new_s1df$ishuman),s1df_a)
colnames(s1df_a)<-c('ishuman','numreads', 'numgenes')

nrplot<-ggplot(data=s1df_a, aes(x=ishuman,y=numreads)) + 
  geom_violin() +
  scale_y_continuous(limits=c(5000,60000));
print(nrplot)

ngplot<-ggplot(data=s1df_a, aes(x=ishuman,y=numgenes)) + 
  geom_violin() +
  scale_y_continuous(limits=c2000,12000));
print(ngplot)

#############################################
#  Try to cluster data
#############################################

s1_dataonly<-new_s1df
s1_dataonly$ishuman<-NULL

#normalize data: 
new<-t(apply(s1_dataonly, 1, function(x)(x/sum(x))))

fit<-kmeans(new,2,iter.max=100)
compare<-cbind(fit$cl,new_s1df$ishuman)

d <- dist(s1_dataonly, method = "euclidean")
fit <- hclust(d, method="ward.D2")
groups <- cutree(fit, k=2)

compare<-cbind(groups,new_s1df$ishuman)
numwrong=sum(compare[,1]-compare[,2]-1)
