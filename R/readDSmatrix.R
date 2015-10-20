readDSmatrix <-
function(filename,...){
  
  temp_df <- read.table(filename, header=TRUE, ...)
  gene_info=temp_df[,1]
  numcells <- dim(temp_df)[2]
  t_matrix <- t(temp_df[,2:numcells])
  df <- data.frame(t_matrix)
  colnames(df)<-gene_info
  #df<-tail(df,-1)
  df<-data.frame(df)
  colnames(df)<-toupper(colnames(df))
  return(df)
  }

#rm(list=setdiff(ls(), "readDSmatrix"))