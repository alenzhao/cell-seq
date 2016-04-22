plot3DPCA<-function(pcadata,cols,ptsize, pttype,viewpt, bbox){
  
  th=viewpt[1]
  ph=viewpt[2]
  fo=viewpt[3]
  zo=viewpt[4]
  
  rgl.open()
  rgl.bg(color = "white")
  rgl.viewpoint(theta=th, phi=ph, fov=fo, zoom=zo)
  x<-pcadata[,1]
  y<-pcadata[,2]
  z<-pcadata[,3]
  
  if (bbox==FALSE) {
    
    plot3d(x,y,z, col=cols, size=ptsize, type=pttype, 
           xlab='PC1', ylab='PC2', zlab='PC3')
  
    } else {
    
    plot3d(x,y,z, col=cols, size=ptsize, type=pttype, 
           xlab='PC1', ylab='PC2', zlab='PC3', 
           xlim=bbox[1:2], ylim=bbox[3:4], zlim=bbox[5:6])
    }

}