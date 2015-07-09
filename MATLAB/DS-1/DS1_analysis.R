
library(R.matlab)
library(data.table)
library(ggplot2)
library(plyr)
library(MultinomialCI)
library(RColorBrewer)
theme_set(theme_gray(24))

## read in files

setwd("~/Documents/Thomson Lab/cell-seq/MATLAB/DS-1");
fname1 <- 'matrix1.txt';
fname2 <- 'matrix2.txt';
df1 <- read.table(fname1, header = TRUE, row.names=1);
df2 <- read.table(fname2, header = TRUE, row.names=1);

## create a summary data frame with: 
## number of reads per barcode
## number of genes expressed per barcode

sumdf1 = colSums(data.matrix(df1))
genesdf1 = colSums(data.matrix(df1)!=0)
sumdf2 = colSums(data.matrix(df2))
genesdf2 = colSums(data.matrix(df2)!=0)

n1 = length(sumdf1);
n2 = length(sumdf2);

# make groupings
newcol = t(rbind(sumdf1, genesdf1, rep(1,length(sumdf1))));
newcol2 = t(rbind(sumdf2, genesdf2, rep(2,length(sumdf2))));
summarydf = data.frame(rbind(newcol,newcol2));
colnames(summarydf) <- c('numreads', 'numgenes', 'group');

# plot violin plots in ggplot2
vplot <- ggplot(summarydf, aes(x=factor(group), y=numreads));
vplot <- vplot+geom_violin()+geom_boxplot();
print(vplot);

# plot violin plots in ggplot2
vgplot <- ggplot(summarydf, aes(x=factor(group), y=numgenes));
vgplot <- vgplot+geom_violin()+geom_boxplot();
print(vgplot);
