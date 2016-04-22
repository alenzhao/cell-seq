preprocessing<-function(df, df_samples,
                        readmin=5000,
                        readmax=200000,
                        geneset='none',
                        minTotReads=20,
                        noisethresh=5,
                        percneeded=0.05,
                        CVratio=1.5,
                        CVplot=FALSE,
                        knowngenes=c(),
                        riboRemove=TRUE){
  
  # subset by reads: 
  keep<-(rowSums(df)>readmin) & (rowSums(df) < readmax)
  df_cellsubset<-df[keep,]
  if (dim(df_samples)[1]==dim(df)[1]){
    df_samples2 <- df_samples[keep,,drop=FALSE]
  }
  
  # cell-normalize data by total number of reads
  normdf<-data.frame(t(apply(df_cellsubset, 1, function(x)(x/sum(x)*1e6))))
  
  # remove noisy genes (i.e. < minReads)
  enough_keep<-apply(df_cellsubset, 2, function(x) if (sum(x)<minTotReads) FALSE else TRUE) #inds to keep based on noisiness
  
  # noise threshold
  signal_inds<-apply(df_cellsubset,2,max)>noisethresh
  
  # keep genes based on gene set:
  keepgeneset<-function(x_d,type){
    df_symbols <- colnames(x_d)
    # make list of TF genes to keep
    TFkeep <- (df_symbols %in% TFsymbols) 
    # make list of ST genes to keep
    STkeep <- (df_symbols %in% STsymbols)
    # make list of cytokine genes to keep
    CYkeep <- (df_symbols %in% CYsymbols) 
    switch(type,
           TF = TFkeep,
           ST = STkeep,
           CY = CYkeep,
           none = (df_symbols %in% df_symbols))
  }
  
  GS_keep<-keepgeneset(df_cellsubset, geneset) # inds to keep based on gene-set
  
  # remove genes that are not often expressed 
  numneeded <- floor(percneeded*dim(df_cellsubset)[1])
  numexpressing <- apply(df_cellsubset,2,function(x)sum(x>0))
  oftenkeep<- (numexpressing >= numneeded)
  
  # # keep genes based on CV:
  if (is.numeric(CVratio)) {
    source('~/Dropbox (Thomson Lab)/analysis/cell-seq/R/plotcv.R')
    CVkeep <- plotcv(figpath,normdf,df_samples,10,CVratio,CVplot) #last parameter is whether or not to save plots
  } else CVkeep <- rep(TRUE,dim(normdf)[2])
   
  # add back genes that were specified, only if they have reads
  hasReads <- colSums(df_cellsubset[,knowngenes]) > minTotReads
  knowngenes <- knowngenes[hasReads]
  knownkeep <- colnames(df) %in% knowngenes
  
  # remove ribosomal genes
  no_r_inds<-!(colnames(df) %in% as.character(rgenes[,1]))
  
  # make a list of all indices to keep 
  all_inds <- (GS_keep & enough_keep & oftenkeep & CVkeep & signal_inds & no_r_inds) | knownkeep
  
  # set df2
  df2<-normdf[,all_inds]
  
  return(list(df2,df_samples2,normdf,df_cellsubset))
}