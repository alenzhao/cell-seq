
library(R.matlab)
library(data.table)
library(ggplot2)
library(plyr)
library(MultinomialCI)
library(RColorBrewer)
theme_set(theme_gray(24))

setwd('~/Documents/Thomson Lab/cell-seq/MATLAB/drop-test/');
filename = 'HundredSTAMPs_MOUSE.csv'

df1 <- read.table(filename, header = TRUE, row.names=1);