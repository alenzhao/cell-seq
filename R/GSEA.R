####################################################
##
## author: sisichen
## date: 10/13/15
##
## This script implements gene set enrichments using 
## the hypergeometric probaility distribution
####################################################

# genelist is a ranked list of genes that are the top genes in a NMF part or PCA vector
# totgenelist: total number of genes in transcriptome
# j: the number of top gene sets to display
# 
# library(R.matlab)
# 
# GSEA<-readMat('/home/iamcam/Documents/cell-seq/R/GSEA_lists_trunc.mat')
# listnames<-unlist(GSEA[[2]])
# listdesc<-unlist(GSEA[[3]])
# lists<-lapply(GSEA[[1]], unlist)

runGSEA <-function(genelist, totgenelist, j) {

  #############################
  ## Read in GSEA lists
  #############################
  
  load('/home/iamcam/Documents/cell-seq/R/GSEA_variables.R')
  
  #############################
  ## Calculate number of genes 
  ## represented in GSEA lists
  #############################
  
  allgenes <- unique(unlist(lists))
  numtot <- length(intersect(allgenes,totgenelist))
  
  # implement checks on numbers
  
  #######################################################
  ## 12 Find overlaps for all genes
  #######################################################
  
  overlap<-sapply(lists, intersect, genelist)
  
  x <- sapply(overlap, length)
  m <- sapply(lists, length)
  k <- length(genelist)
  
  hyperwrapper<-function(x, m, ntot, k){
    
    n <- ntot - m;
    prob <- dhyper(x,m,n,k)
    return(prob)
    
  }
  
  probs <- mapply(hyperwrapper, x, m, numtot, k)
  
  sort_inds<-order(probs, decreasing=FALSE)
  col1<-listnames[sort_inds][1:j]
  col2<-probs[sort_inds][1:j]
  
  result<-cbind(col1,col2)
  
  return(result)

}
