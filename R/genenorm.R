
quantile_normalize <- function(x){
  q <- quantile(x,seq(0,1,0.01))
  names(q)<-seq(0,1,0.01)
  q_norm<-sapply(x,function(y) names(q)[min(which(q>=y))])
  q_norm<-as.numeric(q_norm)
  return(q_norm)
}

genenorm <- function(x_d, type){
  
  switch(type, 
         standard = data.frame(apply(x_d, 2,function(x)(x-mean(x))/(sd(x)))),
         sd = data.frame(apply(x_d, 2,function(x) x/sd(x))),
         mean = data.frame(apply(x_d, 2, function(x) (x/mean(x)))),
         max = data.frame(apply(x_d, 2, function(x) (x/max(x)))),
         quantile = data.frame(apply(x_d, 2, quantile_normalize)),
         scaling = data.frame(apply(x_d, 2, rescale)),
         rootms = data.frame(apply(x_d,2,function(x) (x)/sqrt(sum(x^2)))), # divide by root mean square
         none = x_d)
  
}