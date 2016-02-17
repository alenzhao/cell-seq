plotPCA <- function(figpath, pcadata, genenames, cols){

melted <- cbind(genenames, melt(pcadata$rotation[,1:3]))
scores <- data.frame(cols, pcadata$x[,1:3])
pc1.2 <- qplot(x=PC1, y=PC2, data=scores, colour=cols, size=5) #+ theme(legend.title = element_blank()) + scale_color_manual(values=kcolspec)
pc2.3 <- qplot(x=PC2, y=PC3, data=scores, colour=cols, size=5) #+ theme(legend.title = element_blank()) + scale_color_manual(values=kcolspec)
pc1.3 <- qplot(x=PC1, y=PC3, data=scores, colour=cols, size=5) #+ theme(legend.title = element_blank()) + scale_color_manual(values=kcolspec)

# if (contcolors){
#   pc1.2 = pc1.2 + scale_color_continuous(low=colorscale[1], high=colorscale[2])
#   pc2.3 = pc2.3 + scale_color_continuous(low=colorscale[1], high=colorscale[2])
#   pc1.3 = pc1.3 + scale_color_continuous(low=colorscale[1], high=colorscale[2])
# }

jpeg(paste(figpath,'PC1-2.jpg',sep=''), width=500, height=500, units="px")
print(pc1.2)
dev.off()

jpeg(paste(figpath,'PC2-3.jpg',sep=''), width=500, height=500, units="px")
print(pc2.3)
dev.off()

jpeg(paste(figpath,'PC1-3.jpg',sep=''), width=500, height=500, units="px")
print(pc1.3)
dev.off()

# 
# #make heatmap of cells and the highest/lowest ranked genes in PC1/PC2
# PC1_inds<-sort(df4.pca$rotation[,1], decreasing=TRUE)
# PC2_inds<-sort(df4.pca$rotation[,2], decreasing=TRUE)
# PC3_inds<-sort(df4.pca$rotation[,3], decreasing=TRUE)
# 
# PCdf<-cbind(as.numeric(unname(PC1_inds)),names(PC1_inds), rep('PC1', length(PC1_inds)))
# PCdf<-rbind(PCdf,cbind(as.numeric(unname(PC2_inds)),names(PC2_inds), rep('PC2', length(PC2_inds))))
# PCdf<-rbind(PCdf,cbind(as.numeric(unname(PC3_inds)),names(PC3_inds), rep('PC3', length(PC3_inds))))
# PCdf<-data.frame(PCdf)
# PCdf[,1]<-as.numeric(as.character(PCdf[,1]))
# 
# PC_hist<-ggplot(PCdf,aes(X1)) + facet_grid(X3~.) + geom_histogram();
# 
# jpeg(paste(figpath,'PChist.jpg',sep=''), width=500, height=500, units="px")
# print(PC_hist)
# dev.off()
# 
# return(PCdf)

}

