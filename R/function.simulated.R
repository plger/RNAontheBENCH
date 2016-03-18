#' Transcript-level benchmark of the simulated data
#'
#' Loads and processes the RNAseq data before comparing it with the real foldchanges. 
#' Will produce a number of benchmarking plots and files in the current working directory.
#'
#' @param ANALYSIS_NAME A string indicating the name of the analysis/pipeline. Will be used in filenames, plot titles, etc.
#' @param txfile The path to the transcript-level RNAseq expression matrix. If not given, will look for relevant files in the working directory. The expression matrix should have gene symbols in the first column/row.names, and sample names (either 's1', 's2', etc., or '1' '2', etc..) as column headers.
#' @param qt A string indicating the unit of the expression matrix (either "FPKM", "TPM" or "COUNTS").
#' @param normMethod The normalization method to use (see \code{\link{donorm}}). Defaults to 'TMM'.
#' @param uniquelyMappableLengths Logical, whether to use uniquely mappable length, rather than full length, for FPKM calculation. Default TRUE.
#' @param requireAll logical, whether all samples are required to proceed.
#'
#' @return Nothing, but produces many files in the working directory...
#'
#' @export
compareSimulated <- function(ANALYSIS_NAME, txfile="simulated.quant"){
  if(!file.exists(txfile)){
    message("Could not find expression matrix for simulated data. Skipping comparison.");
  }else{
    e <- read.delim(txfile,header=T,row.names=1)
    data("simulated")
    if(sum(names(e) %in% names(simulated_foldchanges))!=8) names(e) <- paste("s",sub("X","",names(e),fixed=T),sep="")
    if(sum(names(e) %in% names(simulated_foldchanges))!=8)	stop("The column names do not match the simulated samples' names.")

    op <- options(warn=-1)
    on.exit(options(op))
    
    g <- row.names(e)
    g <- g[g %in% row.names(simulated_foldchanges)]
    e <- e[g,names(simulated_foldchanges)]
    fc <- simulated_foldchanges[g,]
    fc2 <- fc
    efc <- as.matrix(e/rowMeans(e[,1:4]))
    ld <- as.numeric(log2(efc)-log2(as.matrix(fc2)))
    png("sim.fc_deviations.png",width=500,height=500)
    layout(matrix(1:2,nrow=2))
    boxplot(as.numeric(ld),horizontal=T,xlab="Deviation in log2(Foldchange)",main="")
    hist(ld[which(ld >-3 & ld < 3)],xlim=c(-3,3),breaks=seq(from=-3,to=3,by=0.1), main="", xlab="Deviation in log2(foldchange), excluding outliers")
    dev.off()
    for(i in 5:8){ fc2[,i] <- 1/fc2[,i] }
    x <- as.numeric(as.matrix(fc2))
    y <- as.numeric(as.matrix(efc))
    w <- which(!is.na(log2(y)) & !is.nan(log2(y)) & !is.infinite(log2(y)))
    x <- x[w]
    y <- y[w]
    png("sim.cor_foldchanges.png",width=500,height=500)
    plot(log2(x),log2(y),pch=4,col=maketrans("black"),xlim=c(-3.5,3.5),ylim=c(-5,5),ylab="RNA-seq log2(foldchange)",xlab="Real log2(foldchange)",main=ANALYSIS_NAME)
    legend("topleft",bty="n",legend=c(paste("Rho:",round(cor(log2(x),log2(y),use="pairwise"),3)),paste("MdAE:",round(MdAE(log2(x),log2(y)),3))))
    abline(a=0,b=1,lty="dashed",col="grey")
    dev.off()
    write(paste(ANALYSIS_NAME,round(cor(log2(x),log2(y),use="pairwise"),4),sep="\t"),"sim.foldchangesCor.txt")
    write(paste(ANALYSIS_NAME,round(MdAE(log2(x),log2(y)),4),sep="\t"),"sim.foldchangesMdAE.txt")
    
    data("annotation")
    g <- g[g %in% row.names(txInfo)]
    txInfo <- txInfo[g,]
    fc2 <- fc2[g,]
    efc <- efc[g,]
    txInfo$cor <- NA
    for(i in 1:nrow(fc2))	txInfo$cor[i] <- suppressWarnings(cor(as.numeric(fc2[i,]),as.numeric(efc[i ,]),use="pairwise"))

    write.table(txInfo, "simulated.stats", col.names=T,row.names=T,sep="\t",quote=F)
    
    genes_colors <- plColorMap(log10(txInfo$Length))
    
    png("sim.cor_vs_nbExons.png",width=500,height=500)
    plot(txInfo$nbExons,txInfo$cor,col=genes_colors,pch=20,xlab="Number of exons",ylab="Cross-sample correlation",xlim=c(1,100))
    dev.off()
    png("sim.cor_vs_length.png",width=500,height=500)
    plot(log10(txInfo$Length),txInfo$cor,col=maketrans("black"),pch=20,xlab="log10(Exonic length in bp)",ylab="Cross-sample correlation")
    dev.off()
    ll <- list()
    for(a in unique(txInfo$nbTxInGene)){
      if(a < 35 & sum(txInfo$nbTxInGene>2))	ll[[a]] <- txInfo$cor[which(txInfo$nbTxInGene==a & !is.na(txInfo$cor))]
    }
    png("sim.cor_vs_nbtx.png",width=500,height=500)
    boxplot(ll,xlab="Number of transcripts in gene",ylab="Cross-sample correlation")
    dev.off()
    
  }
} 
