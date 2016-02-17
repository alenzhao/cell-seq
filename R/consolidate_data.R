consolidate_data <- function(fname, sname, pattern){
  
currfiles<-list.files();
file_inds<-grep(pattern, currfiles)
keepfiles<-sort(currfiles[file_inds])

comb_df<-data.frame();
comb_samples<-c();
for (i in keepfiles) {
  df<-readDSmatrix(i);
  numcells<-dim(df)[1]
  f_label<-strsplit(i, '_')[[1]][1]
  sample<-rep(f_label,numcells)
  comb_samples<-c(comb_samples,sample)
  comb_df<-rbind(comb_df,df)
}

write.table(0, file=fname, sep='\t', eol='\t',quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(t(comb_df), file=fname, append = TRUE, col.names = TRUE, sep='\t', quote=FALSE) 
write.table(comb_samples, file=sname, sep='\t',quote=FALSE, col.names=FALSE, row.names=FALSE)

}