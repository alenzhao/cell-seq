write.df<-function(matrix,filename){
  
  write.table('', filename, sep=',', eol = ',', row.names = FALSE, col.names = FALSE)
  write.table(matrix, filename, append=TRUE, sep=',', quote=FALSE)
  
}