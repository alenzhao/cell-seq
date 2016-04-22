plotPCA <- function(figpath, scores, genenames, colspec, contcolors=FALSE, sampname='sample'){
  
  if (contcolors){
    colorscale<- scale_color_continuous(low='black', high='magenta', name=sampname)
  } else {
    colorscale<-scale_colour_manual(values=colspec, name=sampname)
  }
  
  pca_opts <- ggplot(scores) +
    theme_minimal(base_size=14) + 
    colorscale
  
  colvar<-colnames(scores)[1]
  
  pc1.2 <- pca_opts + geom_point(aes_string('PC1','PC2',color=colvar), size=3)
  pc1.3 <- pca_opts + geom_point(aes_string('PC1','PC3',color=colvar), size=3)
  pc2.3 <- pca_opts + geom_point(aes_string('PC2','PC3',color=colvar), size=3)

  jpeg(paste(figpath,'PC1-2.jpg',sep=''), width=1600, height=1400, pointsize=18, units="px", res=300, quality=100)
  print(pc1.2)
  dev.off()
  
  jpeg(paste(figpath,'PC2-3.jpg',sep=''), width=1600, height=1400, pointsize=18, units="px", res=300, quality=100)
  print(pc2.3)
  dev.off()
  
  jpeg(paste(figpath,'PC1-3.jpg',sep=''), width=1600, height=1400, pointsize=18, units="px", res=300, quality=100)
  print(pc1.3)
  dev.off()

}

