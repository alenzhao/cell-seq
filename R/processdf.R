processdf <- function(df,
                      df_samples,
                      inds=TRUE,
                      cliptype='none',
                      clipmin=0,
                      clipmax=0.05,
                      normtype='mean',
                      readmin = 5000,
                      readmax = 200000,
                      geneset = 'none',
                      minTotReads = 20,
                      noisethresh = 5,
                      percneeded = 0.05,
                      CVratio = 1.5,
                      CVplot = FALSE,
                      knowngenes = c(),
                      riboRemove = TRUE) {
  
  source('~/Documents/Thomson Lab/cell-seq/R/preprocessing.R')
  source('~/Documents/Thomson Lab/cell-seq/R/clip_output.R')
  source('~/Documents/Thomson Lab/cell-seq/R/genenorm.R')
  source('~/Documents/Thomson Lab/cell-seq/R/remZero.R')

  if (length(inds)==1) {
    n <- dim(df_samples)[1]
    inds <- 1:n
  }
  
  result <- preprocessing(
    df[inds,],
    df_samples[inds,,drop=FALSE],
    readmin = readmin,
    readmax = readmax,
    geneset = geneset,
    minTotReads = minTotReads,
    noisethresh = noisethresh,
    percneeded = percneeded,
    CVratio = CVratio,
    CVplot = CVplot,
    knowngenes = knowngenes,
    riboRemove = riboRemove)
  
  df2 <- result[[1]]
  df_samples2 <- result[[2]]
  normdf <- result[[3]]
  df_cellsubset <- result[[4]]
  
  df3 <- clip_output(df2, cliptype, clipmin, clipmax)
  
  df4 <- genenorm(df3, normtype)
  
  varlist<-c('df_samples','df_cellsubset','normdf','df2','df4')
  answer <- vector(mode="list", length=length(varlist))
  names(answer) <- varlist
  answer$df_samples<-df_samples2
  answer$df_cellsubset<-df_cellsubset
  answer$normdf<-normdf
  answer$df2<-df2
  answer$df4<-df4
  
  # answer <- vector(mode="list", length=5)
  # names(answer) <- c('df','df2','df_samples','df3','df4')
  # answer$df<-df
  # answer$df2<-df2
  # answer$df_samples<-df_samples2
  # answer$df3<-df3
  # answer$df4<-df4

  return(answer)
  
}
                      
                      
