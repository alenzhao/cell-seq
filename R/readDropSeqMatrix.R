####################################################
##
## author: sisichen
## date: 07/27/15
##
## This script reads in a matrix from the drop-seq 
## pipeline so that the genes are columns and the 
## cell barcodes are the rows. 
##
####################################################
readDSmatrix<-function(filename,numrows,skiprows){
  
  temp_df <- read.table(filename, header=TRUE, nrows=numrows, skip=skiprows)
  gene_info=temp_df[,1:2]
  numcells <- dim(temp_df)[2];
  m <- temp_df[,2:(numcells+2)]
  t_matrix <- t(m)
  df <- data.frame(t_matrix)
  colnames(df)<-gene_info$X0
  df<-tail(new_s1df,-1)

  }