ggheatmap<-function(H, xname='cell', yname='module', cols=jet.colors(7)) {
  #plots a heatmap with jet colors as default
  
  Hm <- melt(H)
  colnames(Hm) <- c('y','x','value')
  Hm$y<-factor(Hm$y)
  n <- length(unique(Hm$y)) # number of modules
  
  p <- ggplot(Hm, aes(x,y)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradientn(colours = cols) +
    scale_x_discrete(breaks=NULL, name=xname)+
    scale_y_discrete(name=yname)
  
  print(p)
}