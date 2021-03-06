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
# GSEA<-readMat('~/Documents/Thomson Lab/cell-seq/R/GSEA_lists_trunc.mat')
# GO<-readMat('~/Documents/Thomson Lab/cell-seq/R/GOSets.mat')
# GOnames<-unlist(GO$path.desc.all)
# GOlists<-lapply(GO[[2]],unlist)
# GSEAnames<-unlist(GSEA[[2]])
# GSEAlists<-lapply(GSEA[[1]], unlist)
# listnames<-c(GOnames,GSEAnames)
# lists<-c(GOlists,GSEAlists)
# save(listnames,lists,file='GS_variables.R')

runGSEA <-function(genelist, totgenelist, j) {

  #############################
  ## Read in GSEA lists
  #############################
  
  load('~/Documents/Thomson Lab/cell-seq/R/GS_variables.R')
  
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
