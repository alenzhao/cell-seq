
diffexp <- function(df1,df2,genes){
  # returns the differntial level of expression of genes in genelist
  avg1 <- apply(df1, 2, mean)
  avg2 <- apply(df2, 2, mean)
  diffexp <- avg2/avg1
  return(diffexp[genes])
}