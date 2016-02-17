plotcv<-function(figpath, normdf, df_samples, meanthresh, thresh,plotBool){

  new_df <- normdf
  numentries<-apply(new_df,2,function(x)sum(x>0))
  means <- colMeans(new_df)
  sds <- apply(new_df,2,sd)
  vars <- apply(new_df,2,var)
  cvs <- sds/means
  fano <- sds^2/means
  gene <- as.character(colnames(normdf))
  allcvdf <- data.frame(gene,numentries,means,sds,vars, cvs, fano)
  rownames(allcvdf)<-NULL
  
  #remove all data with fewer than 1% entries
  numthresh<-round(0.05*dim(new_df)[1])
  enough_inds<-(allcvdf$numentries>numthresh)
  allcvdf<-allcvdf[enough_inds,]
  
  #remove all data that has zero mean
  allcvdf<-allcvdf[complete.cases(allcvdf),]
  ols <- lm(log(vars)~log(means), data=allcvdf)
  
  #pick genes whose log(cvs2) > upper blue line
  decider <- function(m,c) {
    if (m < meanthresh) {
      return(FALSE)
    } else {
      bound <- ols$coefficient[2]*log(m)+ols$coefficient[1]+log(thresh)
      if(log(c)>bound){
        return(TRUE)  
      } else {
        return(FALSE)
      }
    }
  }
  
  cv_inds<-mapply(decider,allcvdf$means,allcvdf$vars)
  allcvdf<-data.frame(allcvdf,cv_inds)
  hicv<-as.character(allcvdf[cv_inds,'gene'])
  
  if (plotBool) {
    
    jpeg(paste(figpath,'CV_all.jpg',sep=''), width=500, height=500, units="px")
    cvplot2<-ggplot(allcvdf,aes(log(means),log(vars),colour=cv_inds)) + 
      geom_point(alpha=0.2) +
      #geom_abline(slope=ols$coefficient[2], intercept = ols$coefficient[1], colour="red") +f
      geom_abline(slope=ols$coefficient[2], intercept = (ols$coefficient[1] + log(thresh)) , colour='blue') +
      #geom_abline(slope=ols$coefficient[2], intercept = (ols$coefficient[1] - log(thresh)) , colour='blue)' + 
      geom_vline(xintercept=log(meanthresh), colour='blue') + 
      theme(legend.position="none")+
      scale_colour_manual(values=c('black', 'red')) + 
      xlab('log(mean)') + 
      ylab('log(variance)')
    print(cvplot2)
    dev.off()
    
    # plot a spread of genes from high cv genes
    r_x <- round(runif(49, min=0, max=length(hicv)))
    r_genes <- hicv[r_x]
    test <- normdf[,r_genes]
    melted <- melt(test)
    
    jpeg(paste(figpath,'rand_gene_hists.jpg',sep=''), width=1000, height=800, units="px")
    geneplots<-ggplot(melted,aes(log(value)))+ geom_histogram(binwidth=0.1) + facet_wrap(~variable, scales='free')
    print(geneplots)
    dev.off()
    
    # plot a spread of genes from random genes
    r_x <- round(runif(49, min=0, max=dim(allcvdf)[1]))
    r_genes <- as.character(allcvdf[r_x,'gene'])
    test <- normdf[,r_genes]
    melted <- melt(test)
    
    jpeg(paste(figpath,'rand_gene_hists.jpg',sep=''), width=1000, height=800, units="px")
    geneplots<-ggplot(melted,aes(log(value)))+ geom_histogram(binwidth=0.1) + facet_wrap(~variable, scales='free')
    print(geneplots)
    dev.off()
    
    #############################
    ## Make plots of CV by sample
    #############################
    
    cvdf = data.frame()
    for (i in unique(df_samples)){
      inds <- find(df_samples==i)
      new_df <- normdf[inds,]
      means <- colMeans(new_df)
      numentries<-apply(new_df,2,function(x)sum(x>0))
      sds <- apply(new_df,2,sd)
      vars <- apply(new_df,2,var)
      cvs <- sds/means
      fano <- vars/means
      label <- rep(i,length(means))
      gene <- colnames(normdf)
      curr_df <- data.frame(gene,label,numentries,means,sds,vars,cvs, fano)
      
      #remove all data with fewer than 5% entries
      numthresh<-round(0.05*dim(new_df)[1])
      enough_inds<-(curr_df$numentries>numthresh)
      curr_df<-curr_df[enough_inds,]
      
      curr_df <- curr_df[complete.cases(curr_df),]
      rownames(curr_df)<-NULL
      cvdf <- rbind(cvdf, curr_df)
    }
    
  
    jpeg(paste(figpath,'mean_var_by_sample.jpg',sep=''), width=500, height=500, units="px")
    cvplot1<-ggplot(cvdf,aes(log(means),log(vars))) + 
      geom_point(alpha=0.2) + 
      facet_wrap(~label) +
      xlab('log(mean)') + 
      ylab('log(variance)')
    print(cvplot1)
    dev.off()
    
  }
  
  CVkeep<-(colnames(normdf) %in% hicv)
  return(CVkeep)

}