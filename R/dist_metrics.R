#distribution divergence metrics

kldiv<-function(v1,v2){
  # return symmetric KL divergence  
  r<-seq(from = 0, to = max(c(v1,v2))+1, length.out=50)
  h1<-hist(v1,r,plot=FALSE)
  f1<-h1$counts/(sum(h1$counts))
  h2<-hist(v2,r,plot=FALSE)
  f2<-h2$counts/(sum(h2$counts))
  
  kld<-kl.dist(f1,f2)
  return(kld$D)
}

JSD<- function(v1,v2) {
  
  # return symmetric KL divergence  
  r<-seq(from = 0, to = max(c(v1,v2))+1, length.out=50)
  h1<-hist(v1,r,plot=FALSE)
  f1<-h1$counts/(sum(h1$counts))
  h2<-hist(v2,r,plot=FALSE)
  f2<-h2$counts/(sum(h2$counts))
  
  A<-kl.dist(f1,(f1+f2)/2)$D1
  B<-kl.dist(f2,(f1+f2)/2)$D1
  return(sqrt(0.5 * A + 0.5 * B))
}

kstest<-function(v1,v2){
  result<- ks.test(v1,v2)
  return(result$statistic)
}

emd_v<-function(v1,v2){
  # return earth mover distance
  len<-50
  r<-seq(from = 0, to = max(c(v1,v2))+1, length.out=len)
  h1<-hist(v1,r,plot=FALSE)
  f1<-h1$counts/(sum(h1$counts))
  h2<-hist(v2,r,plot=FALSE)
  f2<-h2$counts/(sum(h2$counts))
  
  dist<-emdw(r[1:len-1],f1,r[1:len-1],f2)
  return(dist)
}