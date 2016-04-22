

remZero<-function(currdf){
  if is.data.frame(currdf){
    vector = colSums(currdf)
  }
  # remove zero and Nan genes from a dataframe: 
  bad_inds<-which(vector==0) | is.nan(vector))
  currdf<-currdf[,-bad_inds]
  return(currdf)
}