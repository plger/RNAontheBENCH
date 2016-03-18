#' Gene-level comparison with RT-qPCR data
#'
#' Loads and processes the RNAseq data before comparing it with the RT-qPCR quantification. 
#' Will produce some benchmarking plots and files in the current working directory.
#'
#' @param ANALYSIS_NAME A string indicating the name of the analysis/pipeline. Will be used in filenames, plot titles, etc.
#' @param rnaseq The path to the gene-level RNAseq expression matrix. If not given, will look for relevant files in the working directory. The expression matrix should have gene symbols in the first column/row.names, and sample names (e.g. 'AJ80', etc) as column headers.
#' @param qt A string indicating the unit of the expression matrix (either "FPKM", "TPM" or "COUNTS").
#'
#' @return Nothing, but produces many files in the working directory...
#'
#' @export
compareWithPCR <- function(ANALYSIS_NAME,rnaseq=NULL,qt){
  sl <- list(analysis=ANALYSIS_NAME)
  qt <- match.arg(toupper(qt),c("FPKM","TPM","COUNTS"))
    
  if(is.null(rnaseq)){
    r <- .getGeneQt()
  }else{
    if(class(rnaseq)=="character"){
        r <- read.delim(rnaseq,header=T,row.names=1)
    }else{
        r <- rnaseq
        remove(rnaseq)
    }
  }
  
  dinf <- .checkDataset(r)
  r <- r[,dinf$samples]

  if(dinf$dataset != "w12"){
    warning("RT-qPCR measurements are available only for the 12-samples dataset.")
    return(NULL)
  }
  
  
  if(qt=="COUNTS"){
    # values are counts - we first convert to FPKM
    r <- counts2fpkmWrapper(r,"gene")
    qt <- "FPKM"
  }  
  if(qt!="TPM"){
    for(i in 1:ncol(r))	r[,i] <- fpkm2tpm(r[,i])
    qt <- "TPM"
  }
 
  r <- donorm(r[,unique(names(r))], "TMM")
  
  data("qpcr")  
  
  # remove genes and samples that were not quantified by both methods
  samples <- names(qpcr)[which(names(qpcr) %in% names(r))]
  genes <- row.names(qpcr)[which(row.names(qpcr) %in% row.names(r))]
  # put everything in the same order
  r <- r[,samples]
  rn <- r[genes,]
  q <- qpcr[genes,samples]

  if(!("GAPDH" %in% row.names(r)) | any(as.numeric(r["GAPDH",])==0)){
    message("Could not find GAPDH expression levels. Skipping comparison with qPCR...")
    return(NULL)
  }
  gapdh <- as.numeric(r["GAPDH",])
  r1 <- rn
  for(i in 1:ncol(rn)) r1[,i] <- r1[,i]/gapdh[i]
 
op <- options(warn=-1)
  on.exit(options(op))

 
  a <- function(x){ as.numeric(as.matrix(x))}
 
  cols <- c("black","blue","red","yellow","pink")
  png("w12.correlation_with_qPCR.png",height=500,width=500)
  plot(a(q),a(r1),xlab="qPCR (% of GAPDH)",ylab="RNA-seq (% of GAPDH)",log="xy",col=cols,pch=20,cex=1.5,main=ANALYSIS_NAME)
  legend("topleft",bty="n",legend=c("Pearson correlation","(before log-transform)",round(cor(as.numeric(as.matrix(q)),as.numeric(as.matrix(r1))),3)))
  dev.off()
  sl[["cor.qPCR"]] <- cor(as.numeric(as.matrix(q)),as.numeric(as.matrix(r1)))
  
  q2 <- fc2mean(q)
  r2 <- fc2mean(rn)
  png("w12.fc2mean_qPCR.png",height=500,width=500)
  plot(log2(a(q2)),log2(a(r2)),xlab="qPCR log2(Foldchange to the mean)",ylab="RNA-seq log2(Foldchange to the mean)",main=ANALYSIS_NAME,pch=20,cex=1.5,col=cols)
  legend("topleft",bty="n",legend=c("Pearson correlation","(before log-transform)",round(cor(a(q2),a(r2)),3)))
  dev.off()
  sl[["cor.qPCR.fc2mean"]] <- cor(as.numeric(as.matrix(q2)),as.numeric(as.matrix(r2)))
  
  write.table(as.data.frame(sl),"qPCR.stats",col.names=T,row.names=F,sep="\t",quote=F)

}
