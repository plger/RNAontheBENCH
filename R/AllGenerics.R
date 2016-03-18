#' isCleanData
#'
#' Returns which couples of x and y values are clean (non-missing, non-infinite and non-NaN)
#'
#' @param x a vector
#' @param y a vector of same length as x
#'
#' @return indices that are 'clean' (non-missing, non-infinite and non-NaN) in both x and y.
#'
#' @examples
#' isCleanData(c(2,NA,3),c(0,5,NA))
#'
#' @export
isCleanData <- function(x,y){
  if(length(x) != length(y))	stop("x and y do not have the same length!")
  return(which(!is.infinite(x) & !is.infinite(y) & !is.nan(x) & !is.nan(y) & !is.na(y) & !is.na(x)))
}


#' checkMatIntegrity
#'
#' Replaces missing and invalid values with 0, and makes sure the matrix has all rows and samples in 'against'.
#'
#' @param x a matrix or data.frame
#' @param against (optional) another matrix or data.frame, of which all rows and columns should be present in x.
#'
#' @return a clean version of x, or an error if there are missing columns.
#'
#' @examples
#' checkMatIntegrity(matrix(c(1:8,NA),nrow=3))
#'
#' @export
checkMatIntegrity <- function(x,against=NULL){
  x <- as.matrix(x[which(!is.na(row.names(x))),])
  if(!is.null(against)){
    if(ncol(x)<ncol(against) | !all(colnames(against) %in% colnames(x)))	stop("The expression matrix misses some samples!")
    missingFeatures <- sum(!(row.names(against) %in% row.names(x)))
    if(missingFeatures>0){
      message(paste(missingFeatures," features are missing from the expression matrix. They will be added with zero expression.",sep=""))
      d2 <- matrix(rep(0,ncol(x)*missingFeatures),nrow=missingFeatures)
      row.names(d2) <- row.names(against)[which( !(row.names(against) %in% row.names(x)) )]
      x <- rbind(x,d2)
    }
  }
  if(any(is.na(x)) | any(is.nan(x)) | any(is.infinite(x)) | any(x<0)){
    message("The expression matrix contains missing/invalid values. These will be replaced by zero.")
    x[which(is.na(x) | is.nan(x) | is.infinite(x) | x<0)] <- 0
  }
  if(any(is.infinite(x)))	stop("The expression matrix contains infinite values!")
  return(as.data.frame(x))
}

#' Median absolute error
#'
#' Returns the median absolute error between vectors x and y.
#'
#' @param x a numeric vector
#' @param y a numeric vector of same length as x.
#'
#' @return The median absolute error between x and y.
#'
#' @examples
#' MdAE(1:3,jitter(1:3))
#'
#' @export
MdAE <- function(x,y){
  return(median(abs(y-x)))
}

#' Make transparent
#'
#' Adds transparency to a color.
#'
#' @param x a color specification (see \code{\link[grDevices]{col2rgb}}
#' @param alpha a numeric value between 0 and 255 indicating the degree of transparency.
#'
#' @return A color code.
#'
#' @examples
#' maketrans("blue")
#'
#' @export
maketrans <- function(tcol,alpha=100){
  c <- col2rgb(tcol)
  rgb(c["red",1][[1]],c["green",1][[1]],c["blue",1][[1]],alpha,maxColorValue=255)
}

#' normalizes a dataset
#'
#' Wrapper to normalize a dataset, using the given method.
#'
#' @param dataset is a data.frame
#' @param method is a character string of either 'linear', 'housekeeping', 'quantile', or any method supported by edgeR's \code{\link[edgeR]{calcNormFactors}}.
#'
#' @return The normalized data.frame.
#'
#' @examples
#' donorm(matrix(1:12,nrow=4),"linear")
#'
#' @export
donorm <- function(dataset, method="TMM"){
  method <- match.arg(method, c("TMM","RLE","quantile","linear","housekeeping","none"))
  if(method=="none")	return(dataset)
  if(method=="housekeeping"){
    # defines housekeeping genes
    housekeepers <- data.frame(gene=c("GAPDH","TUBB","TBP"), transcript=c("NM_002046","NM_001172085","NM_178014"))
    # check whether this is tx- or gene-level data
    w1 <- which(row.names(dataset) %in% as.character(housekeepers$gene))
    w2 <- which(row.names(dataset) %in% as.character(housekeepers$transcript))
    if(length(w1)==0 & length(w2)==0)	stop("No housekeeping gene found for normalization!")
    if(length(w1) > length(w2)){
      normFactors <- getLinearNormalizers(log2(dataset[w1,]))
    }else{
      normFactors <- getLinearNormalizers(log2(dataset[w2,]))
    }
    return(as.data.frame(t(2^(log2(t(dataset))*normFactors))))
  }
  if(method=="linear"){
    nf <- getLinearNormalizers(dataset)
    return(t(t(dataset)/nf))
  }
  if(method=="quantile"){
    return(preprocessCore::normalize.quantiles(dataset))
  }
  # call edgeR
  library(edgeR)
  d <- DGEList(counts=dataset,group=colnames(dataset))
  d <- calcNormFactors(d, method=method)
  return(as.data.frame(cpm(d, normalized.lib.sizes=T)*mean(d$samples$lib.size)/1000000))
}


#' z-score
#'
#' Wrapper to calculate the z-score
#'
#' @param x a numeric matrix or data.frame.
#'
#' @return The matrix of corresponding row z-scores.
#'
#' @examples
#' zscore(matrix(1:12,nrow=3))
#'
#' @export
zscore <- function(x){
  t(scale(t(x)))
}


#' foldchange-to-the-mean
#'
#' Returns the foldchange to the row's mean for each value.
#'
#' @param x a numeric matrix or data.frame.
#'
#' @return The matrix of corresponding foldchanges to the row's mean.
#'
#' @examples
#' fc2mean(matrix(1:12,nrow=3))
#'
#' @export
fc2mean <- function(x){ 
  for(i in 1:nrow(x))    x[i,] <- as.numeric(x[i,])/mean(as.numeric(x[i,]),na.rm=T)
  x
}

#' plColorMap
#'
#' Maps colors onto numeric values.
#'
#' @param x a numeric vector.
#'
#' @return A vector of colors of the same length as x.
#'
#' @examples
#' plColorMap(1:3)
#'
#' @export
plColorMap <- function(x){
  pal <- colorRampPalette(c("blue","red"),0.5)(100)
  xmin <- min(x,na.rm=T)
  xmax <- max(x,na.rm=T)-xmin
  sapply(x,FUN=function(y){pal[1+floor((99)*(y-xmin)/xmax)]})
}

#' getLinearNormalizers
#'
#' Returns linear scale factors for each column of the dataset, using the first column as a reference and fitting a linear model between each other colum and the first.
#'
#' @param dataset a numeric data.frame or matrix.
#' @param tryrobust logical, whether to try to fit a robust linear model instead of a normal one.
#'
#' @return A vector of length ncol(dataset) containing the scale factors of each column of dataset.
#'
#' @examples
#' getLinearNormalizers(matrix(1:12,nrow=3))
#'
#' @export
getLinearNormalizers <- function(dataset, tryrobust=F){
  normBy <- as.numeric(dataset[,1])
  co <- c(1,apply(dataset[,2:ncol(dataset)],2,FUN=function(x){ mod <- plGetModel(x,normBy,tryrobust=tryrobust); suppressWarnings(if(is.na(mod)) return(NA)); as.numeric(coef(mod)[[1]])}))
  names(co) <- colnames(dataset)
  co
}

#' plTryRead
#'
#' Tries to find and read a file in the working directory, using an ordered list of regular expression and returning the first matching file.
#'
#' @param patterns a vector of regular expressions.
#'
#' @return The content of a tab-delimited file, if any matches the set of expressions, otherwise NULL.
#'
#' @examples
#' plTryRead(c("test"))
#'
#' @export
plTryRead <- function(patterns){
  for(p in patterns){
    tmp <- list.files(pattern=p)
    if(length(tmp)>0) return(read.delim(tmp[[1]],header=T,row.names=1))
  }
  return(NULL)
}

# load the transcript-level quantification file
.getTxQt <- function(){
  tmp <- plTryRead(c("transcripts.quant", "isoforms.fpkm_table", "transcript.counts", "transcript.fpkm", "transcript.tpm", "transcripts.tpm", "transcripts.fpkm","tpm$"));
  if(!is.null(tmp))	return(tmp)
  stop("Could not find RNA-seq quantification file.")
}

# load the gene-level quantification file
.getGeneQt <- function(){
  tmp <- plTryRead(c("genes.quant", "genes.fpkm_table", "gene.counts", "gene.fpkm", "gene.tpm", "genes.tpm"));
  if(!is.null(tmp))	return(tmp)
  message("Could not find the gene-level RNA-seq quantification file. Trying to build it by summing transcripts.")
  return(convertTx2Genes(.getTxQt()))
}

# load the spike-in quantification
.getSpikeinQt <- function(){
  tmp <- plTryRead(c("spikein.quant","spikein.fpkm", "spikein.counts", "spikein.tpm"));
  if(!is.null(tmp))	return(tmp)
  if(file.exists("genes.quant")){
	system('head -1 genes.quant > spikein.quant; grep "^ERCC-" genes.quant >> spikein.quant')
	return(getSpikeinQt())
  }
  return(.getSpikeinsFromTranscripts())
}

# load the spike-in quantification using the tx-level quantification matrix
.getSpikeinsFromTranscripts <- function(fname=NULL){
  if(is.null(fname)){
    e <- .getTxQt()
  }else{
    if(class(fname)=="character"){
        e <- read.delim(fname,header=T,row.names=1)
    }else{
        e <- fname
    }
  }
  data("ercc")
  m <- merge(ercc[,1:2],e,by.x="TX.ID",by.y="row.names")
  m$TX.ID <- NULL
  row.names(m) <- m$ERCC.ID
  m$ERCC.ID <- NULL
  return(m)
}

#' Summarizes transcript quantifications to the gene level
#'
#' Summarizes transcript quantifications to the gene level by summing transcript values belonging to the same gene.
#'
#' @param tx A data.frame of transcript values, with refseq IDs as row.names.
#'
#' @return A data.frame of gene values.
#'
#' @export
convertTx2Genes <- function(tx){
  if(class(tx)=="character") tx <- read.delim(tx,header=T,row.names=1)
  data("annotation")
  t <- intersect(row.names(tx),t2g$V2)
  if(length(t) < nrow(tx))	message(paste(nrow(tx)-length(t),"transcripts were not associated to any gene."))
  tx <- tx[t,]
  row.names(t2g) <- t2g$V2
  t2g <- t2g[row.names(tx),]
  ag <- aggregate(tx,by=list(gene=t2g$V1),na.rm=T,FUN=sum)
  row.names(ag) <- ag$gene
  ag$gene <- NULL
  return(ag)
}

#' Convert FPKM values to TPM values
#'
#' Convert FPKM values to transcripts per million (TPM).
#'
#' @param fpkm A numeric vector of fpkm values.
#'
#' @return A numeric vector of the corresponding TPM values.
#'
#' @export
fpkm2tpm <- function(fpkm){
  tpm <- exp(log(fpkm) - log(sum(fpkm,na.rm=T)) + log(1e6))
  tpm[which(is.na(tpm))] <- 0
  return(tpm)
}

#' Convert counts to FPKM values
#'
#' Convert fragment counts to FPKM values, using the sum as library size. By default, effective length is used, set mean.frag.size to 0 to use total (given) length.
#'
#' @param counts A numeric vector of counts (for one sample).
#' @param lengths A numeric vector of same length as 'counts', indicating the length of each transcript.
#' @param mean.frag.size A numeric value indicating the mean fragment size (default 220) used for effective length calculation.
#'
#' @return A numeric vector of the corresponding FPKM values.
#'
#' @export
counts2fpkm <- function(counts, lengths, mean.frag.size=220){
  lengths <- lengths-mean.frag.size
  counts[which(lengths < 0)] <- 0
  lengths[which(lengths < mean.frag.size)] <- mean.frag.size
  10^9*counts/(lengths*sum(counts))
}


#' Retrieves a linear model, if possible.
#'
#' Retrieves a linear model, if possible.
#'
#' @param x a numeric vector containing the independent variable.
#' @param y a numeric vector of same length as x, containing the dependent variable.
#' @param no.intercept logical, indicating whether to disable the intersect in the linear model (defaults TRUE, i.e. models must go through the origin)
#' @param tryrobust logical, indicating whether to try to fit a robust linear model (see \code{\link[MASS]{rlm}}) instead of a normal one.
#'
#' @return Returns the linear model, if the fit was possible, NA otherwise.
#'
#' @export
plGetModel <- function(x,y, no.intercept=T, tryrobust=F){
  if(tryrobust)	library(MASS)
  x <- as.numeric(x)
  y <- as.numeric(y)
  mod <- NA
  if(no.intercept){
    if(tryrobust) mod <- try(rlm(y~x+0), silent=T)
    if(!tryrobust | !("lm" %in% class(mod))) mod <- try(lm(y~x+0), silent=T)
  }else{
    if(tryrobust) mod <- try(rlm(y~x), silent=T)
    if(!tryrobust | !("lm" %in% class(mod))) mod <- try(lm(y~x), silent=T)
  }
  if("lm" %in% class(mod))	return(mod)
  return(NA)
}

#' Retrieves the log2(foldchange)
#'
#' Retrieves the log2(foldchange).
#'
#' @param x a numeric vector.
#' @param groups a vector of same length as x, indicating to which of two groups each value belongs. If more than two groups, only the first two are used.
#'
#' @return Returns the log2 foldchange, or the log2 of the largest mean if one of the mean is 0, or 0 if both means are 0.
#'
#' @export
getLogFC <- function(x, groups){
  m1 <- mean(as.numeric(x[which(groups==unique(groups)[[1]])]),na.rm=T)
  m2 <- mean(as.numeric(x[which(groups==unique(groups)[[2]])]),na.rm=T)
  if(m1==0 & m2==0)	return(0)
  if(m2==0)	return(log2(1/m1))
  if(m1==0)	return(log2(m2))
  log2(m2/m1)
}

#' Convert counts matrix to FPKM matrix
#'
#' Convert fragment counts to FPKM values, using the column sums as library sizes and effective length (see \code{\link{counts2fpkm}}).
#'
#' @param r A numeric matrix or data.frame with samples as columns and genes/transcripts as rows.
#' @param level Character, either "gene" or "tx", used to fetch the right lengths.
#' @param uniquelyMappableLengths Whether to use the uniquely mappable lengths for transcripts (default TRUE).
#'
#' @return A matrix/data.frame of the corresponding FPKM values.
#'
#' @export
counts2fpkmWrapper <- function(r,level="transcript",uniquelyMappableLengths=NULL){
  data("annotation")
  level <- match.arg(level, c("transcript","gene"))
  if(is.null(uniquelyMappableLengths))	uniquelyMappableLengths <- level=="transcript"
  if(level=="transcript" & uniquelyMappableLengths){
    message("Using uniquely-mappable-lengths (instead of whole transcript length) to calculate FPKM")
    if("Length" %in% names(r))	r$Length <- NULL
    r <- merge(unique_lengths, r, by.x="V1",by.y="row.names")
    row.names(r) <- r$V1
    r$V1 <- NULL
    names(r)[1] <- "Length"
    fsize <- 0
  }else{
    if(tolower(names(r)[1])!="length"){
      if(level=="tx"){
	r <- merge(txInfo,r,by="row.names")
      }else{
	r <- merge(geneLengths,r,by="row.names")
      }
      row.names(r) <- r$Row.names
      r$Row.names <- NULL
    }
    fsize <- 220
  }
  for(i in 2:ncol(r))  r[,i] <- counts2fpkm(r[,i],r$Length,fsize)
  r$Length <- NULL
  return(r)
}

#' Writes a table of metrics comparing two vectors
#'
#' @param x A numeric vector of observed values.
#' @param y A numeric vector of expected values (same length as 'x').
#' @param ANALYSIS_NAME The name of the analysis pipeline
#' @param fname The filename for saving the metrics
#' @param clean Logical, whether to clean the data (default TRUE), see \code{\link{isCleanData}}.
#' @param norm Logical, whether to linearly normalize the two vectors before comparison (default FALSE).
#' @param incZeros Logical, whether to keep 0 values (default TRUE).
#'
#' @return Nothing.
#'
#' @export
simStats <- function(x,y,ANALYSIS_NAME, fname,clean=T,norm=F,incZeros=F){
  if(clean){
    w <- isCleanData(x,y)
    x <- x[w]
    y <- y[w]
  }
  if(incZeros){
    x[which(x==0)] <- min(x[which(x>0)])/2
    y[which(y==0)] <- min(y[which(y>0)])/2
  }
  if(.checkPkg(rminer)){
    library(rminer)
    svg(paste(fname,"REC.svg",sep="."),width=4,height=4)
      mgraph2(y,x,graph="REC",PTS=100,main=paste(ANALYSIS_NAME,fname,sep=" - "))
    dev.off()
    png(paste(fname,"REC.png",sep="."),width=400,height=400)
      mgraph2(y,x,graph="REC",PTS=100,main=paste(ANALYSIS_NAME,fname,sep=" - "))
    dev.off()
  }
  if(norm){
    mod <- plGetModel(y,x)
    y2 <- y/as.numeric(coef(mod))
  }else{
    y2 <- y
  }
  m <- list();
#  rec <- RECcurve(y-x)  # replaces rminer code
#  write.table(as.data.frame(t(c(ANALYSIS_NAME, fname, as.numeric(rec[,1])))),paste(fname,"REC",sep="."),row.names=F,col.names=F,sep="\t",quote=F)  
  m[["COR"]] <- cor(x,y)
  m[["MdAE"]] <- MdAE(x,y2)
  m[["SP"]] <- cor(x,y,method="spearman")
  m[["euclidean"]] <- sqrt(sum((x-y2)^2))
  med <- median(y,na.rm=T)
  m[["COR.belowM"]] <- cor(x[which(y<=med)],y[which(y<=med)])
  m[["COR.aboveM"]] <- cor(x[which(y>med)],y[which(y>med)])
  m[["SP.belowM"]] <- cor(x[which(y<=med)],y[which(y<=med)],method="spearman")
  m[["SP.aboveM"]] <- cor(x[which(y>med)],y[which(y>med)],method="spearman")
  m[["MdAE.belowM"]] <- MdAE(x[which(y<=med)],y[which(y<=med)])
  m[["MdAE.aboveM"]] <- MdAE(x[which(y>med)],y[which(y>med)])
  m[["euclidean.belowM"]] <- sqrt(sum((x[which(y<=med)]-y2[which(y<=med)])^2))
  m[["euclidean.aboveM"]] <- sqrt(sum((x[which(y>med)]-y2[which(y>med)])^2))

  m2 <- as.numeric(m)
  names(m2) <- names(m)
  write.table(as.data.frame(t(c(ANALYSIS_NAME,m2))),paste(fname,"metrics",sep="."),col.names=T,row.names=F,sep="\t",quote=F)
}

#' A wrapper that performs the whole benchmark analysis
#'
#' Performs the whole series of benchmark analysis (see \code{\link{txLevel}}, \code{\link{compareWithNanostring}}, \code{\link{analyzeSpikein}}, \code{\link{compareWithPCR}}, and \code{\link{compareSimulated}}).
#' The function expects quantification files with the right column headers ("AJ80" and so on, or "s1","s2" and so on for the simulated data) to be in the rpath folder, and to bear some kind of recognizable name.
#' To avoid confusion, we suggest using the filename 'transcripts.quant' for transcript-level quantification, 'genes.quant' for gene-level quantification (optional), and 'simulated.quant' for transcript-level quantification of the simulated dataset (optional).
#' If you want to specify files manually, use the individual underlying functions.
#'
#' @param rpath The path were the quantification files are stored, and where the output files will be saved.
#' @param ANALYSIS_NAME The name of the analysis pipeline
#' @param qt A string indicating the unit of the expression matrix (either "FPKM", "TPM" or "COUNTS").
#'
#' @return Nothing but saves a bunch of files in 'rpath'.
#'
#' @examples
#' # first we create a directory and put the example quantification file in it:
#' data(exampledata)
#' dir.create("example")
#' write.table(exampleTranscriptLevel,"w12.transcripts.quant",col.names=T,row.names=T,sep="\t",quote=F)
#' write.table(exampleGeneLevel,"w12.genes.quant",col.names=T,row.names=T,sep="\t",quote=F)
#' # run the wrapper, specifying that folder:
#' benchmarkWrapper("example", "tophat.featureCount", qt="COUNTS")
#'
#' @export
benchmarkWrapper <- function(rpath, ANALYSIS_NAME, qt){
  qt <- match.arg(toupper(qt), c("FPKM","COUNTS","TPM"))
  owd <- getwd()
  setwd(rpath)
  data("benchmark_report")
  if(file.exists("w12.transcripts.quant")){
    message("### Analyzing the 12-samples dataset...")
    message("### - spike-in quantification...")
    analyzeSpikein(ANALYSIS_NAME, "w12.transcripts.quant", qt=qt)
    message("### - transcript-level quantification...")
    compareWithNanostring(ANALYSIS_NAME, "w12.transcripts.quant", qt=qt)
    message("### - gene-level quantification...")
    if(file.exists("w12.genes.quant")){
      compareWithNanostring(ANALYSIS_NAME, "w12.genes.quant", qt=qt)
      message("### - comparing with RT-qPCR...")
      compareWithPCR(ANALYSIS_NAME, "w12.genes.quant", qt=qt)
    }else{
      message("(gene-level file not found - creating it from the transcript-level...)")
      compareWithNanostring(ANALYSIS_NAME, convertTx2Genes("w12.transcripts.quant"), qt=qt)
      message("### - comparing with RT-qPCR...")
      compareWithPCR(ANALYSIS_NAME, convertTx2Genes("w12.transcripts.quant"), qt=qt)
    }
    message("")
  }
  if(file.exists("w6.transcripts.quant")){
    message("### Analyzing the 6-samples dataset...")
    message("### - spike-in quantification...")
    analyzeSpikein(ANALYSIS_NAME, "w6.transcripts.quant", qt=qt)
    message("### - transcript-level quantification...")
    compareWithNanostring(ANALYSIS_NAME, "w6.transcripts.quant", qt=qt)
    message("### - gene-level quantification...")
    if(file.exists("w6.genes.quant")){
      compareWithNanostring(ANALYSIS_NAME, "w6.genes.quant", qt=qt)
    }else{
      message("(gene-level file not found - creating it from the transcript-level...)")
      compareWithNanostring(ANALYSIS_NAME, convertTx2Genes("w6.transcripts.quant"), qt=qt)
    }
    message("")
  }
  if(!file.exists("w12.transcripts.quant") & !file.exists("w6.transcripts.quant")){
    message("Attempting to find files...")
    message("### Analyzing spike-in quantification...")
    analyzeSpikein(ANALYSIS_NAME, qt=qt)
    message("### Comparing gene-level quantification with Nanostring...")
    compareWithNanostring(ANALYSIS_NAME, .getGeneQt(), qt=qt)
    message("### Comparing transcript-level quantification with Nanostring...")
    compareWithNanostring(ANALYSIS_NAME, .getTxQt(), qt=qt)
    message("### Comparing gene-level quantification with qPCR...")
    compareWithPCR(ANALYSIS_NAME, qt=qt)
    message("")
  }
  if(file.exists("simulated.quant")){
    message("### Comparing quantification of the simulated dataset with ground truth...")
    compareSimulated(ANALYSIS_NAME)
    sim_report <- benchmark_has_sim
    message("")
  }else{
    sim_report <- benchmark_no_sim
  }
  tabheaders <- ""
  tabcontents <- ""
  if(file.exists(.plfilename("w12",ANALYSIS_NAME,"transcript.stats"))){
    tabheaders <- paste(tabheaders, gsub("$$DATASET$$","w12",benchmark_tabheaders,fixed=T))
    tabcontents <- paste(tabcontents, gsub("$$PREFIX$$",.plfilename("w12",ANALYSIS_NAME,""),sub("$$QPCR$$",benchmark_qpcr,benchmark_main_tabs,fixed=T),fixed=T))
  }
  if(file.exists(.plfilename("w6",ANALYSIS_NAME,"transcript.stats"))){
    tabheaders <- paste(tabheaders, gsub("$$DATASET$$","w6",benchmark_tabheaders,fixed=T))
    tabcontents <- paste(tabcontents, gsub("$$PREFIX$$",.plfilename("w6",ANALYSIS_NAME),"",sub("$$QPCR$$","",benchmark_main_tabs,fixed=T),fixed=T))
  }
  write(benchmark_js,"lib.js")
  o <- gsub("$$ANALYSIS_NAME$$",ANALYSIS_NAME,benchmark_page,fixed=T)
  o <- sub("$$TAB_HEADERS$$",tabheaders,o,fixed=T)
  o <- sub("$$SIMULATED_CONTENT$$",sim_report,o,fixed=T)
  o <- sub("$$TABS$$",tabcontents,o,fixed=T)
  write(o,"index.html")
  message("Analysis done!")
  message(paste("See results in ",getwd(),"/index.html", sep=""))
  browseURL("index.html")
  setwd(owd)
}




.checkPkg <- function(package){
  return(tryCatch(require(package, character.only=TRUE),error = function(e) FALSE))
}

# returns the type of the dataset
.checkDataset <- function(x, analysis="",checksamples=T,checklevel=T){
  di <- list(analysis=analysis)
  if(checklevel){
    if(length(grep("^NM_",row.names(x))) > 0 | length(grep("^DQ516748",row.names(x)))>0){
      di[["level"]] <- "transcript"
    }else{
      if(length(grep("EIF4H",row.names(x))) == 1 | length(grep("^ERCC-",row.names(x)))>1){
	di[["level"]] <- "gene"
      }else{
	stop("Row names appear to be neither refseq ids nor gene symbols!")
      }
    }
  }else{
    di[["level"]] <- NULL
  }
  if(checksamples){
    di[["sample_colors"]] <- "black"
    if(length(grep("^AJ",colnames(x)))>1){
      di[["dataset"]] <- "w12"
      di[["samples"]] <- c("AJ80","AJ81","AJ82","AJ83","AJ84","AJ86","AJ87","AJ89","AJ90","AJ91","AJ92","AJ93")
      di[["sample_colors"]] <- c('#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02','#A6761D','#666666','black','orange','lightgoldenrod',"violet")
    }else{
      if(length(grep("shCTR$",colnames(x)))==3){
	di[["dataset"]] <- "w6"
	di[["samples"]] <- c("w306o.shCTR","w306o.sh2","c3391S.shCTR","c3391S.sh2","CFG.shCTR","CFG.sh2")
	di[["sample_colors"]] <- c( "#000000FF", "#FF0000FF", "#00CD00FF", "#0000FFFF", "#BEBEBEFF", "#FF00FFFF")
      }else{
        if(all(paste("s",gsub("s","",colnames(x),fixed=T),sep="") %in% paste("s",1:8,sep=""))){
            colnames(x) <- paste("s",gsub("s","",colnames(x),fixed=T),sep="")
            di[["samples"]] <- paste("s",1:8,sep="")
            di[["sample_colors"]] <- c( "#000000FF", "#FF0000FF", "#00CD00FF", "#0000FFFF", "#00FFFFFF", "#FF00FFFF", "#FFFF00FF", "#BEBEBEFF")
            di[["dataset"]] <- "sim"
	}else{
          if(all(colnames(x) %in% paste(rep(LETTERS[1:4],each=5),rep(1:5,4),sep="_"))){
            di[["samples"]] <- paste(rep(LETTERS[1:2],each=5),rep(1:5,4),sep="_")
            di[["dataset"]] <- "seqc"
          }else{
            stop("Could not identify the dataset. Make sure that the column names are as they should be.")
          }
	}
      }
    }
    if(!all(di[["samples"]] %in% colnames(x)))	stop(paste("There are some missing samples! Expected samples:\n", paste(di[["samples"]], collapse=", ")))
  }else{
    di[["samples"]] <- colnames(x)
  }
  di[["ptitle"]] <- paste(di$dataset,":",di$analysis," (",di$level,")",sep="")
  di[["fname"]] <- .plfilename(di$dataset,di$analysis,di$level)
  return(di)
}

.plfilename <- function(...){
  gsub("[^[:alnum:]_.]", "", gsub(" ","_",paste(...,collapse=".",sep=".")))
}
