####################################################
##
## author: sisichen
## date: 08/28/15
##
## This script is for playing around with dimensionality
## reduction algorithms and hierarchical clustering
##
####################################################

library(ggplot2)
library(plyr)
library(R.matlab)
library(grep)
library(reshape2)
library(MASS)
library(RColorBrewer)
library(gridExtra)
library(dplyr)

x=matrix(c(exp(-0.2*(-(1:300)/10))*cos(-(1:300)/10),
           exp(-0.2*(-(1:300)/10))*sin(-(1:300)/10)),
         ncol=2)

plot(x)

fit.all = prcomp(x)
approx.all=fit.all$x[,1] %*% t(fit.all$rotation[,1])
plot(x,xlab=expression(x[1]),ylab=expression(x[2]))
points(approx.all,pch=4)

fit = prcomp(x[270:280,])
pca.approx = fit$x[,1]%*%t(fit$rotation[,1])+colMeans(x[270:280,])
plot(x[270:280,], xlab=expression(x[1]),ylab=expression(x[2]))
points(pca.approx,pch=4)

##########################
## Code from Lecture Notes
##########################

# Find the k smallest entries in each row of an array
# Inputs: n*p array, p>=k, number of smallest entries to find
# Output: n*k array of column indices for smallest entries per row 
smallest.by.rows<-function(m,k) {
  stopifnot(ncol(m)>=k) # otherwise k smallest is meaningless
  row.orders = t(apply(m,1,order))
  k.smallest = row.orders[,1:k]
  return(k.smallest)
}









# Find multiple nearest neighbors in a data frame
# Inputs: n*p matrix of data vectors, number of neighbors to find, # optional arguments to dist function
# Calls: smallest.by.rows
# Output: n*k matrix of the indices of nearest neighbors
find.kNNs<-function(x,k,...){
  x.distances = dist(x,...)
  x.distances=as.matrix(x.distances)
  kNNs = smallest.by.rows(x.distances,k+1)
  return(kNNs[,-1])
}



# Least-squares weights for linear approx. of data from neighbors
# Inputs: n*p matrix of vectors, n*k matrix of neighbor indices, # scalar regularization setting
# Calls: local.weights
# Outputs n*n matrix of weights