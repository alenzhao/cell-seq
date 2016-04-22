findEndpts <- function (data){
  
  colnames(data)<-c('x','y','z')
  meanX = colMeans(data)
  
  result<-prcomp(data)
  dirVect = result$rotation[,1]
  t1 = c(min(result$x[,1]-0.2), max(result$x[,1]+0.2))
  endpts = rbind(meanX + t1[1]*t(dirVect), meanX + t1[2]*t(dirVect))
  
  return(endpts)
}