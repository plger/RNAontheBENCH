#' Comparison quantification with the expected spike-in concentrations.
#'
#' Loads and processes the RNAseq data before comparing it with the spike-in concentrations.
#' Will produce a number of benchmarking plots and files in the current working directory.
#'
#' @param ANALYSIS_NAME A string indicating the name of the analysis/pipeline. Will be used in filenames, plot titles, etc.
#' @param rnaseq The path to the gene-level RNAseq expression matrix. If not given, will look for relevant files in the working directory. The expression matrix should have gene symbols in the first column/row.names, and sample names (e.g. 'AJ80' for the 12-samples dataset, "A_1" for SEQC, etc) as column headers.
#' @param qt A string indicating the unit of the expression matrix (either "FPKM", "TPM" or "COUNTS").
#' @param fc.undetected The foldchange to assign to undetected spike-ins (should be either 1 or NA, default 1)
#'
#' @return Nothing, but produces many files in the working directory...
#'
#' @examples
#' # first we create a directory and put the example quantification file in it:
#' data(exampledata)
#' dir.create("example")
#' write.table(exampleGeneLevel,"w12.genes.quant",sep="\t",quote=F)
#' # then we run the function, giving a name to the analysis, 
#' # specifying the file and type of quantification:
#' analyzeSpikein("tophat.featureCount", "w12.genes.quant", qt="COUNTS")
#'
#' @export
analyzeSpikein <- function(ANALYSIS_NAME, rnaseq=NULL, qt, fc.undetected=1){
  sl <- list(analysis=ANALYSIS_NAME)
  qt <- match.arg(toupper(qt),c("FPKM","TPM","COUNTS"))
  
  if(is.null(rnaseq)){
    sp <- .getSpikeinQt()
  }else{
    sp <- read.delim(rnaseq,header=T,row.names=1)
    if("DQ854994" %in% row.names(sp))   sp <- .getSpikeinsFromTranscripts(sp)
  }
  dinf <- .checkDataset(sp, ANALYSIS_NAME)
  if(dinf[["dataset"]] == "sim")    stop("There are no spike-ins to analyze in the simulated data.")
  dinf$ptitle <- paste(dinf$dataset,":",dinf$analysis," (spike-ins)",sep="")
  dinf$fname <- .plfilename(dinf$dataset,dinf$analysis,"spikein")
  
  data("ercc")
  if(dinf$level=="gene"){
    row.names(ercc) <- ercc$ERCC.ID
  }else{
    row.names(ercc) <- ercc$TX.ID
  }

  if(qt=="COUNTS"){
    # values are counts - we first convert to FPKM
    sp <- counts2fpkmWrapper(sp,dinf$level,uniquelyMappableLengths=F)
    qt <- "FPKM"
  }  
  if(qt != "TPM"){
    for(i in 1:ncol(sp))	sp[,i] <- fpkm2tpm(sp[,i])
    qt <- "TPM"
  }

  # ensure that all spike-ins are there, set to 0 the missing ones
  sp <- merge(sp,ercc,by="row.names",all.y=T)[,1:(ncol(sp)+1)]
  row.names(sp) <- sp$Row.names
  sp$Row.names <- NULL
  sp <- as.matrix(sp)
  sp[which(is.na(sp) | is.nan(sp) | is.infinite(sp) | sp<0)] <- 0

  # put everything in the same order
  genes <- row.names(sp)[which(row.names(sp) %in% row.names(ercc))]
  sp <- as.data.frame(sp[genes,dinf$samples])
  ercc  <- ercc[genes,]

  # in case it is unknown, determine which sample has which mix
  df <- data.frame(row.names=names(sp))
  df$mix1.cor <- apply(sp,2,FUN=function(x){ cor(as.numeric(x),as.numeric(ercc$concentration.in.Mix.1..attomoles.ul.))})
  df$mix2.cor <- apply(sp,2,FUN=function(x){ cor(as.numeric(x),as.numeric(ercc$concentration.in.Mix.2..attomoles.ul.))})
  df$best.matching.mix <- apply(df,1,FUN=function(x){ order(as.numeric(x),decreasing=T)[1]})

  # renormalize samples according to their expected spike-in amounts, 
  # thereby avoiding discrepancies due to the different amounts of spike-ins put in the cells.
  # We first make a copy of the spike-in mesurements, which we'll multiply by FC between mixes to make an artificial situation were the samples have the same mix:
  sp2 <- sp
  for(i in which(df$best.matching.mix==2))	sp2[,i] <- as.numeric(sp2[,i])*as.numeric(ercc$expected.fold.change.ratio)
  # we then use this to get linear normalizers
  nf <- getLinearNormalizers(sp2)
  # which we can then apply to the original measurements
  spn <- t(t(sp)*nf)
  # to check that everything worked as planned, we can compare the sum of euclidean distance between the normalized and un-normalized matrices:
  # sum(dist(t(sp))) > sum(dist(t(spn)))
  # alternatively, we check that mix assignment worked:
  if(dinf$dataset=="w12"){
    mix1 <- c("AJ86","AJ90","AJ91")
  }else{
    if(dinf$dataset=="w6"){
        mix1 <- c("CFG.sh2","c3391S.sh2","w306o.shCTR","w306o.sh2")
    }else{
        mix1 <- paste("A",1:5,sep="_")
    }
  }
  df$expected.mix <- ifelse(row.names(df) %in% mix1, 1, 2)
  if(any(df$best.matching.mix != df$expected.mix)){
    warning(paste(ANALYSIS_NAME,"There appears to be a problem with the spike-in quantification or normalization, because samples do not correlate with the correct spike-in mixes. This can happen with very poor quantifications, can be due to a mis-labeling of the columns, or wrong specification of the samples having mix1. Analysis will proceed, but you should double-check the samples' mix assignment:",sep=" - "),immediate.=T)
    print(df)
  }
  
  # next, we must create the x values (expected spike-in concentrations), according to the mix injected in each sample
  xdf <- ercc[,df$expected.mix+4]
  op <- options(warn=-1)
  on.exit(options(op))

  
  #calculate by-sample correlation
  sc <- data.frame(analysis=ANALYSIS_NAME, sample=colnames(spn))
  sc$cor <- sapply(1:ncol(spn),FUN=function(i){ 
    x <- as.numeric(xdf[,i])
    y <- as.numeric(spn[,i])
    w <- which(x>0 & y>0)
    cor(x[w],y[w])
  })
  write.table(sc,paste(dinf$fname,".perSampleCorrelation.tab",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
  
  x <- as.numeric(t(xdf))
  
  # to include also undetected genes in the plot, put them at half of the minimum detected FPKM
  med <- median(x)
  finc <- 0.001
  y <- as.numeric(t(spn))
  w <- which(y>0)	# w indicates the genes that were detected by RNA-seq
  
  # calculate correlations
  sl[["cor.detected"]] <- cor(x[w],y[w],use="pairwise")
  sl[["log.cor.detected"]] <- cor(log10(x[w]),log10(y[w]),use="pairwise")
  sl[["cor.all"]] <- cor(x,y,use="pairwise")
  sl[["cor.detected.abovemedian"]] <- try(cor(x[which(y>med)], y[which(y>med)],use="pairwise"),silent=T)
  sl[["cor.detected.belowmedian"]] <- try(cor(x[which(y>0 & y<=med)], y[which(y>0 & y<=med)],use="pairwise"),silent=T)
  sl[["prop.undetected"]] <- 1-sum(spn>0)/ncol(spn)/92
  ww <- which(!as.logical(lapply(sl,FUN=is.numeric)))
  for(wo in ww) sl[[wo]] <- NA
  simStats(x,y,dinf$ptitle,paste(dinf$fname,"expr",sep="."),norm=T)
  
  # We can now plot the correlation between (normalized) FPKM and expected values
  samples_colors <- sapply(dinf$sample_colors, alpha=150, FUN=maketrans)
  png(paste(dinf$fname,"cor.png",sep="."),width=700,height=600)
  plot(x,y+finc,xlab="Expected spike-in concentration (attomoles/ul)", ylab=qt, col=samples_colors,pch=20,cex=1.5,log="xy",main=dinf$ptitle)
  legend("topleft",bty="n",legend=paste(names(sl[2:6]),lapply(sl[2:6],digits=3,FUN=round),sep=": "))
  legend("bottomright",bty="n",legend=c(paste(round(100*sl[["prop.undetected"]],2),"% undetected",sep="")))
  dev.off()

  # we average the samples according to their mix, and calculate the observed foldchange between mixes
  d <- data.frame(expected=ercc$expected.fold.change.ratio)
  d$avg.mix1=apply(spn[,df$expected.mix==1],1,na.rm=T,FUN=mean)
  d$avg.mix2=apply(spn[,df$expected.mix==2],1,na.rm=T,FUN=mean)
  d$observed <- foldchange(d$avg.mix2,d$avg.mix1,fc.undetected)
  
  w <- isCleanData(d$observed,d$expected)
  sl[["cor.foldchange"]] <- cor(d$observed[w], d$expected[w])
  w <- isCleanData(log2(d$observed),log2(d$expected))
  sl[["cor.log2.foldchange"]] <- cor(log2(d$observed[w]), log2(d$expected[w]), use="pairwise")

  simStats(d$observed,d$expected,dinf$ptitle,paste(dinf$fname,"foldchanges",sep="."),norm=F)
  
  # we calculate statistical significance of the difference between mixes
  pvals <- apply(spn, 1, FUN=function(x){ t.test(as.numeric(x[df$expected.mix==1]),as.numeric(x[df$expected.mix==2]))$p.value})
  sl[["significant"]] <- sum(pvals < 0.05 & ercc$expected.fold.change.ratio != 1, na.rm=T)
  sl[["false.positives"]] <- sum(pvals < 0.05 & ercc$expected.fold.change.ratio == 1, na.rm=T)
  
  # for plotting, we will color the transcripts according to their real concentration
  genes_colors <- plColorMap(log10(rowMeans(ercc[,5:6],na.rm=T))) 
  
  # plot observed vs expected foldchanges
  png(paste(dinf$fname,"corFC.png",sep="."),width=500,height=500)
  plot(log2(d$expected), log2(d$observed), pch=20, cex=1.2, col=genes_colors, xlab="Exepected log2 foldchange (mix1/mix2)", ylab="Observed log2 foldchange",main=dinf$ptitle)
  abline(a=0,b=1,lty="dashed",col="grey")
  legend("bottomright",bty="n",legend=c( paste("foldchanges cor:",round(sl[["cor.foldchange"]],3)), 
					paste("log2(foldchanges) cor:",round(sl[["cor.log2.foldchange"]],3)),
					"",
					paste(round(100*sl[["significant"]]/sum(ercc$expected.fold.change.ratio!=1)),"% of real differences have p<0.05",sep="")
				  ))
  dev.off()

  ld <- log2(d$expected)-log2(d$observed)
  png(paste(dinf$fname,"FC_deviations.png",sep="."),width=500,height=500)
  layout(matrix(1:2,nrow=2))
  boxplot(ld,horizontal=T,xlab="Deviation in log2(Foldchange)",main=dinf$ptitle)
  hist(ld[which(ld >-3 & ld < 3)],xlim=c(-3,3),breaks=seq(from=-3,to=3,by=0.2), main="", xlab="Deviation in log2(foldchange), excluding outliers")
  dev.off()

  write.table(as.data.frame(sl),paste(dinf$fname,"stats",sep="."),col.names=T,row.names=F,sep="\t",quote=F)

}
 
