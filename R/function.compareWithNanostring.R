#' Comparison with Nanostring
#'
#' Loads and processes the RNAseq data before comparing it with the Nanostring quantification. 
#' Will produce a number of benchmarking plots and files in the current working directory.
#'
#' @param ANALYSIS_NAME A string indicating the name of the analysis/pipeline. Will be used in filenames, plot titles, etc.
#' @param rnaseq The path to the gene-level RNAseq expression matrix. If not given, will look for relevant files in the working directory. The expression matrix should have refseq id or gene symbols in the first column/row.names, and sample names (e.g. 'AJ80') as column headers.
#' @param qt A string indicating the unit of the expression matrix (either "FPKM", "TPM" or "COUNTS").
#' @param normMethod The normalization method to use (see \code{\link{donorm}}). Defaults to 'housekeeping' for gene-level, and 'TMM' for transcript-level.
#'
#' @return Nothing, but produces many files in the working directory...
#'
#' @examples
#' # first we create a directory and put the example quantification file in it:
#' data(exampledata)
#' dir.create("example")
#' write.table(exampleGeneLevel,"w12.genes.quant",sep="\t",quote=FALSE)
#' # then we run the function, giving a name to the analysis, 
#' # specifying the file and type of quantification:
#' compareWithNanostring("tophat.featureCount", "w12.genes.quant", qt="COUNTS")
#'
#' @export
compareWithNanostring <- function(ANALYSIS_NAME, rnaseq=NULL, qt, normMethod=NULL){
  qt <- match.arg(toupper(qt),c("FPKM","TPM","COUNTS"))

  if(is.null(rnaseq)){
    ra <- .getGeneQt()
  }else{
    if(class(rnaseq)=="character"){
      ra <- read.delim(rnaseq,header=T,row.names=1)
    }else{
      ra <- rnaseq
      rm(rnaseq)
    }
  }
  dinf <- .checkDataset(ra, ANALYSIS_NAME)
  ra <- ra[,dinf$samples]

  if(dinf$dataset %in% c("sim","seqc"))	stop("This appears to be either the simulated data, or the SEQC data, both of which cannot be compared to Nanostring...")
  
  if(is.null(normMethod)){
    if(dinf$level=="gene"){
      normMethod <- "housekeeping"
    }else{
      normMethod <- "TMM"
    }
  }
  message(paste("Normalization method:",normMethod))
  
  if(dinf$dataset == "w12"){
    data("w12_nanostring")
  }else{
    data("w6_nanostring")
  }
  
  if(dinf$level=="gene"){
    # create a nanostring table that aggregates (sums) isoforms of the same gene
    na <- aggregate(nanostring[,3:ncol(nanostring)],by=list(gene=nanostring$gene),FUN=sum)
    row.names(na) <- na$gene
    na$gene <- NULL
  }else{
    na <- nanostring[,3:ncol(nanostring)]
    row.names(na) <- nanostring$probe
    row.names(nanoprobes) <- nanoprobes$probe
    nanoprobes <- nanoprobes[row.names(na),]    
  }
  na <- na[,dinf$samples]

  if(qt=="COUNTS"){
    # values are counts - we first convert to FPKM
    ra <- counts2fpkmWrapper(ra,dinf$level)
    qt <- "FPKM"
  }  
  if(qt!="TPM"){
    for(i in 1:ncol(ra))	ra[,i] <- fpkm2tpm(ra[,i])
    qt <- "TPM"
  }
  ra <- donorm(ra[,unique(names(ra))], normMethod)

  if(dinf$level=="transcript"){
    # summarize transcripts by probe
    rn <- t(sapply(strsplit(nanoprobes$isoforms.hit,";",fixed=T),FUN=function(x){ colSums(ra[as.character(x),],na.rm=T) }))
    row.names(rn) <- nanoprobes$probe
  }else{
    rn <- ra
  }

  rn <- checkMatIntegrity(rn,na)  

  # put everything in the same order
  genes <- row.names(na)[which(row.names(na) %in% row.names(rn))]
  na <- na[genes,]
  rn <- rn[genes,]

  if(any(colSums(rn)==0))   stop("Some of the samples appear to have no non-null quantifications!")
  
  op <- options(warn=-1)
  on.exit(options(op))
  
  sl <- list(analysis=ANALYSIS_NAME, level=dinf$level)

  gmeds <- apply(na,1,na.rm=T,FUN=median)

  samples_colors <- sapply(dinf$sample_colors, alpha=150, FUN=maketrans)
  
  x <- as.numeric(t(rn))
  y <- as.numeric(t(na))
  w <- isCleanData(x,y)

  sl[["prop.undetected"]] <- 1-sum(x>0)/length(x)
  
  # for the purpose of plotting (i.e. to avoid weird plots), values are set to a minimum
  x2plot <- x
  x2plot[x<0.001] <- 0.001

  med <- median(y,na.rm=T)
  sl[["cor.detected"]] <- cor(x[w],y[w],use="pairwise")
  sl[["cor.log.detected"]] <- cor(log10(x[x>0]),log10(y[x>0]),use="pairwise")
  sl[["cor.all"]] <- cor(x,y,use="pairwise")
  sl[["cor.abovemedian"]] <- cor(x[which(y>med)], y[which(y>med)],use="pairwise")
  if(sum(x <= med & x>0)>0){
    sl[["cor.belowmedian"]] <- try(cor(x[which(y<=med & x>0)], y[which(y<=med & x>0)],use="pairwise"), silent=T)
  }else{
    sl[["cor.belowmedian"]] <- NA
  }
  
  png(paste(dinf$fname,"nano2seq.png",sep="."),width=800,height=600)
  plot(x2plot,y,xlab=qt,ylab="Nanostring intensity",col=samples_colors,pch=20,cex=1.5,log="xy",main=dinf$ptitle)
  legend("topleft",bty="n",legend=c("Pearson correlation",
      paste("Cor detected, before log-transform:",round(sl[["cor.detected"]],4)),
      paste("Cor detected, log-transformed:",round(sl[["cor.log.detected"]],4)),
      paste("Cor all:",round(sl[["cor.all"]],4)),
      paste("Cor below median:",round(sl[["cor.belowmedian"]],4)),
      paste("Cor above median:",round(sl[["cor.abovemedian"]],4))
    ))
  dev.off()

  #calculate by-sample correlation
  sc <- data.frame(analysis=ANALYSIS_NAME, sample=names(na))
  sc$cor <- sapply(1:ncol(na),FUN=function(i){ 
    x <- as.numeric(na[,i])
    y <- as.numeric(rn[,i])
    w <- isCleanData(x,y)
    cor(x[w],y[w])
  })
  write.table(sc,paste(dinf$fname,"perSampleCorrelation.tab",sep="."),col.names=T,row.names=T,sep="\t",quote=F)
  
  simStats(x,y,dinf$ptitle,paste(dinf$fname,"expr",sep="."),norm=T)  

  # get gene colors by mean expression level
  genes_colors <- plColorMap(log10(rowMeans(na)+0.01))  
  
  if(dinf$level == "transcript"){
    # summarize txInfos by probe
    data("annotation")
    tstats <- t(sapply(strsplit(nanoprobes$isoforms.hit,";",fixed=T),FUN=function(x){ c(mean(txInfo[as.character(x),"Length"],na.rm=T),median(txInfo[as.character(x),"nbExons"],na.rm=T),max(txInfo[as.character(x),"nbTxInGene"],na.rm=T))}))
    row.names(tstats) <- nanoprobes$probe
    colnames(tstats) <- c("mLength","mdExons","nbTxInGene")
    tstats <- as.data.frame(tstats[row.names(na),])
    tstats$medianExpression <- apply(na,1,na.rm=T,FUN=median)
    tstats$cor <- NA
    
    for(i in 1:nrow(rn))	tstats$cor[i] <- cor(as.numeric(rn[i,]),as.numeric(na[i ,]),use="pairwise")
    write.table(tstats, paste(dinf$fname,"txInfos",sep="."),col.names=T,row.names=T,sep="\t",quote=F)

    png(paste(dinf$fname,"cor_vs_nbExons.png",sep="."),width=500,height=500)
    plot(tstats$mdExons,tstats$cor,col=genes_colors,main=dinf$ptitle,pch=20,cex=1.5,xlab="Number of exons",ylab="Cross-sample correlation between RNA-seq and Nanostring")
    legend("bottomright",bty="n",legend=c("Correlation:",round(cor(tstats$mdExons,tstats$cor),2)))
    dev.off()
    png(paste(dinf$fname,"cor_vs_length.png",sep="."),width=500,height=500)
    plot(tstats$mLength,tstats$cor,col=genes_colors,main=dinf$ptitle,pch=20,cex=1.5,xlab="Exonic length",ylab="Cross-sample correlation between RNA-seq and Nanostring")
    legend("bottomright",bty="n",legend=c("Correlation:",round(cor(tstats$mLength,tstats$cor),2)))
    dev.off()
    png(paste(dinf$fname,"cor_vs_expression.png",sep="."),width=500,height=500)
    plot(tstats$medianExpression,tstats$cor,main=dinf$ptitle,col=genes_colors,pch=20,cex=1.5,xlab="Median Nanostring expression",ylab="Cross-sample correlation between RNA-seq and Nanostring")
    legend("bottomright",bty="n",legend=c("Correlation:",round(cor(tstats$medianExpression,tstats$cor),2)))
    dev.off()
    png(paste(dinf$fname,"cor_vs_nbtx.png",sep="."),width=500,height=500)
    plot(tstats$nbTxInGene,tstats$cor,main=dinf$ptitle,col=genes_colors,cex=1.5,pch=20,xlim=c(0,35),xlab="Number of transcripts in gene",ylab="Cross-sample correlation between RNA-seq and Nanostring")
    dev.off()
  }
  
  
  # plot foldchange-to-the-mean, color by expression level
  png(paste(dinf$fname,"foldchange2mean.png",sep="."),width=800,height=600)
  x <- as.numeric(as.matrix(fc2mean(rn)))
  y <- as.numeric(as.matrix(fc2mean(na)))
  plot(log2(x),log2(y),xlab="log2(RNA-seq foldchange to the mean)",ylab="log2(Nanostring foldchange to the mean)",col=genes_colors,pch=20,cex=1.5,main=dinf$ptitle)
  w2 <- which(x!=0 & !is.infinite(x))
  w <- isCleanData(x,y)
  abline(a=0,b=1,lty="dashed",col="darkgrey")
  mod <- plGetModel(log2(x[w2]),log2(y[w2]),no.intercept=T)
  if("lm" %in% class(mod))   abline(mod, col="red", lty="dashed")
  legend("topleft",bty="n",legend=c("Pearson correlation",paste("Before log-transform:",round(cor(x[w],y[w],use="pairwise"),4)),paste("After log:",round(cor(log2(x[w2]),log2(y[w2]),use="pairwise"),4))))
  dev.off()
  sl[["fc.cor"]] <- cor(x[w],y[w],use="pairwise")
  sl[["logfc.cor"]] <- cor(log2(x[w2]),log2(y[w2]),use="pairwise")
  
  simStats(x,y,dinf$ptitle,paste(dinf$fname,"fc2mean",sep="."),norm=F)  

  x <- as.numeric(as.matrix(fc2mean(rn[which(gmeds>median(gmeds)),])))
  y <- as.numeric(as.matrix(fc2mean(na[which(gmeds>median(gmeds)),])))
  sl[["fc.cor.abovemedian"]] <- cor(x[which(!is.infinite(x))], y[which(!is.infinite(x))], use="pairwise")
  if(sum(gmeds > 0 & gmeds <= median(gmeds))>0){
    x <- as.numeric(as.matrix(fc2mean(rn[which(gmeds > 0 & gmeds<=median(gmeds)),])))
    y <- as.numeric(as.matrix(fc2mean(na[which(gmeds > 0 & gmeds<=median(gmeds)),])))
    sl[["fc.cor.belowmedian"]] <- cor(x[which(!is.infinite(x))], y[which(!is.infinite(x))], use="pairwise")
  }else{
    sl[["fc.cor.belowmedian"]] <- NA
  }

  x <- as.numeric(as.matrix(zscore(rn)))
  y <- as.numeric(as.matrix(zscore(na)))

  # plot z-scores deviations
  ld <- x-y
  png(paste(dinf$fname,".z_deviations.png",sep=""),width=500,height=500)
  layout(matrix(1:2,nrow=2))
  boxplot(ld,horizontal=T,xlab="Deviation in z-score",main=dinf$ptitle)
  hist(ld[which(ld >-3 & ld < 3)],xlim=c(-3,3),breaks=seq(from=-3,to=3,by=0.2), main="", xlab="Deviation in z-score, excluding outliers")
  dev.off()

  simStats(x,y,dinf$ptitle,paste(dinf$fname,"zscores",sep="."),norm=F)  

  # plot z scores
  png(paste(dinf$fname,"zscores.png",sep="."),width=800,height=600)
  plot(x,y,xlab="RNA-seq Z-score",ylab="Nanostring Z-score",col=genes_colors,pch=20,cex=1.5,main=dinf$ptitle)
  abline(a=0,b=1,lty="dashed",col="darkgrey")
  mod <- plGetModel(x,y,no.intercept=T)
  if("lm" %in% class(mod))   abline(mod, col="red", lty="dashed")
  legend("topleft",bty="n",legend=c("Pearson correlation",round(cor(x,y,use="pairwise"),4)))
  dev.off()
  sl[["z.cor"]] <- cor(x,y,use="pairwise")
  x <- as.numeric(as.matrix(zscore(rn[which(gmeds>median(gmeds)),])))
  y <- as.numeric(as.matrix(zscore(na[which(gmeds>median(gmeds)),])))
  sl[["z.cor.abovemedian"]] <- cor(x[which(!is.infinite(x))], y[which(!is.infinite(x))], use="pairwise")
  if(sum(gmeds > 0 & gmeds <= median(gmeds))>0){  
    x <- as.numeric(as.matrix(zscore(rn[which(gmeds > 0 & gmeds<=median(gmeds)),])))
    y <- as.numeric(as.matrix(zscore(na[which(gmeds > 0 & gmeds<=median(gmeds)),])))
    sl[["z.cor.belowmedian"]] <- cor(x[which(!is.infinite(x))], y[which(!is.infinite(x))], use="pairwise")
  }else{
    sl[["z.cor.belowmedian"]] <- NA
  }
    
  if(dinf$dataset=="w12"){
    # calculate median correlation with copy-number of WBS genes
    data("copynumbers")
    if(dinf$level=="transcript"){
      cns <- tx_copynumbers
    }else{
      cns <- copynumbers
    }
    ra <- ra[which(row.names(ra) %in% row.names(cns)),dinf$samples]
    cns <- cns[row.names(ra),dinf$samples]
    cc <- vector("numeric",nrow(cns))
    for(i in 1:nrow(ra)){ cc[i] <- cor(as.numeric(cns[i,]),as.numeric(ra[i,]))}
    png(paste(dinf$fname,"corWithCopynumber.png",sep="."),width=500,height=500)
    hist(cc,xlab=paste(dinf$level,"s' correlation with copy-number",sep=""),breaks=20,col="grey",main=dinf$ptitle)
    legend("topleft",bty="n",legend=c(paste("median",round(median(cc, na.rm=T),3)), paste("mean",round(mean(cc, na.rm=T),3))))
    dev.off()
    sl[["median.cor.with.copynumber"]] <- median(cc, na.rm=T)
    sl[["avg.cor.with.copynumber"]] <- mean(cc, na.rm=T)
  }
  
  # check for relative abundances of isoforms
  eif4h <- c("NM_022170_var1.1","NM_031992.1")
  if(dinf$level == "transcript" & all(eif4h %in% row.names(rn)) & all(as.numeric(as.matrix(rn[eif4h,]))>0)){
    png(paste(dinf$fname,"EIF4H.relativeAbundance.png",sep="."),width=400,height=400)
    rr <- as.numeric(rn[eif4h[1],]/colSums(rn[eif4h,]))
    nr <- as.numeric(na[eif4h[1],]/colSums(na[eif4h,]))
    sl[["relAbundance.EIF4H"]] <- median(as.numeric(rn[eif4h[1],]/colSums(rn[eif4h,])))
    sl[["cor.rel.eif4h.isoforms"]] <- cor(as.numeric(rn[eif4h[1],]/colSums(rn[eif4h,])),as.numeric(na[eif4h[1],]/colSums(na[eif4h,])))
    plot(rr,nr,main=dinf$ptitle,ylab="Nanostring relative abundance of EIF4H isoforms",xlab="RNA-seq relative abundance of EIF4H isoforms",pch=20,cex=1.2)
    abline(a=0,b=1,lty="dashed",col="grey")
    legend("bottomright",bty="n",legend=c(paste("R^2:", round(sl[["cor.rel.eif4h.isoforms"]],3))))
    dev.off()
  }  
  
  write.table(as.data.frame(sl),paste(dinf$fname,"stats",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
    
}
 
