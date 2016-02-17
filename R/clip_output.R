clip_output<-function(data,type, minperc, maxperc){
  
#clip top and bottom 10% of data, gene-wise
clipgene<- function(x, minperc, maxperc) {
  xmax<-max(x)
  xmin<-min(x)
  
  xmaxthresh<-max(x)-maxperc*(xmax-xmin)
  xminthresh<-min(x)+minperc*(xmax-xmin)
  
  clip<-function(y) {
    if (y>xmaxthresh) {
      xmaxthresh 
    } else if (y<xminthresh) {
      xminthresh
    } else y
  }
  
  listout<-sapply(x, clip) 
}

#whole matrix max and mins
mmax<-max(max(data))
mmin<-min(min(data))
mmax_thresh<-mmax - maxperc*(mmax-mmin)
mmin_thresh<-mmin + minperc*(mmax-mmin)

clipmatrix<-function(x){
  #x is a column vector
  
  x[x>mmax_thresh]<-mmax_thresh
  x[x<mmin_thresh]<-mmin_thresh
  return(x)
  
}

output<-switch(type,
       bygene=data.frame(apply(data,2,clipgene, minperc,maxperc)),
       bymatrix=data.frame(apply(data,2,clipmatrix)),
       none=data)

return(output)  

}