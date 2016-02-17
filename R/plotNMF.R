plotNMF<- function(figpath, df, df.nmf, nmfBool,r,n, inds){

  W<-df.nmf@fit@W
  H<-df.nmf@fit@H
  
  Wnames<-apply(W, 2, function(x) (names(sort(x, decreasing=TRUE))[1:n]))
  Wreshape<-reshape(Wnames,n*r,1)
  
  jpeg(paste(figpath,'NMF_df_by_W.jpg',sep=''), width=1000, height=1000, units="px")
  imagesc(t(as.matrix(df[inds,Wreshape])))
  dev.off()
  
  jpeg(paste(figpath,'NMF_H.jpg',sep=''), width=1000, height=1000, units="px")
  imagesc(H[,inds])
  #imagesc(rbind(H[,clusterInds], clusters[clusterInds]*25))
  dev.off()
  
  return(Wnames)
}

