####################################################
##
## author: sisichen
## date: 07/22/15
##
## This script plots statistics on the polyT tails in 
## R1 of a drop-seq run.
##
####################################################
library(ggplot2)
library(plyr)
library(R.matlab)
library(grep)
library(reshape2)
library(MASS)

setwd("~/Documents/RNAseqdata/DS-3")
read1s = read.table("50cells_R1_trimmed.fastq",stringsAsFactors = FALSE)

setwd("~/Documents/RNAseqdata/072815_BeadQC/")
read1s = read.table("output.fastq",stringsAsFactors = FALSE)


read1subset= read1s[1:50000,]

i_list<-gregexpr('TTTTTTTTTT',read1subset)
i<-sapply(i_list, "[[", 1);
i=i-3; #to account for TAC
idf<-data.frame(i)
colnames(idf)<-c('indexofpolyT')
head(idf)

# summarize the number of occurrences of indices
i_summ <- data.frame(table(idf))
colnames(i_summ)<- c('indexofpolyT','Freq')

# calculate expected number of occurrences if BCs are random
f<-function(x) {
  if (x<=21)
    (1/4)^(21-x)
  else
    0
}
# make new index_summary df
prob_dist = sapply(1:dim(i_summ)[1],f)
i_summ$E_Freq <- i_summ$Freq*prob_dist
i_summ<-melt(i_summ, "indexofpolyT")
i_summ$indexofpolyT = as.numeric(i_summ$indexofpolyT);

T_hist<-ggplot(i_summ,aes(indexofpolyT,value,fill=variable)) +
  geom_bar(stat="identity", position="dodge") +
  labs(title="index of start of polyT", x="index of polyT start", y='frequency') +
  scale_x_continuous(limits=c(15,25),breaks=15:25)
plot(T_hist)

##################################################
###  Calculate difference between Freq and E_Freq
##################################################

Freq = i_summ[which(i_summ$variable=='Freq'),]$value
E_Freq = i_summ[which(i_summ$variable=='E_Freq'),]$value

diff_Freq <- Freq-E_Freq
tofit <- rev(diff_Freq[1:20])/sum(Freq)

f <-fitdistr(tofit, "poisson")
