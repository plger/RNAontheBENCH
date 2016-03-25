# #####
# functions to perform differential expression analysis and perform benchmarking plots on the results on the basis of the spike-ins or nanostring data
# For sleuth, should be run from within the folder containing the kallisto/salmon/sailfish results subfolders (with bootstrap); then use "sleuthWrapper" (below)
# For other, load your count data matrix (or other) and use "deSpikein" (below)
# #####


#' Runs a spike-in differential expression analysis.
#'
#' Runs a differential expression analysis and compares it to real differences between spike-in mixes. 
#' Sleuth is handled in a different function (see \code{\link{sleuthWrapper}})  
#' A warning is given if the specified/expected mix distribution does not match the observed one.
#'
#' @param dat The counts matrix or data.frame, with gene symbols or transcript Refseq IDs as row.names, and sample names as column headers.
#' @param method The differential expression method to use. Either 'edgeR', 'DESeq', 'DESeq2', 'EBSeq', 'voom', 't' (t-test), or 'logt' (t-test on log-transformed values). A string indicating the unit of the expression matrix (either "FPKM", "TPM" or "COUNTS").
#' @param norm The normalization method to use (either "linear", or any of the methods supported by edgeR's \code{\link[edgeR]{calcNormFactors}}. Defaults to 'TMM'.
#' @param quantification A string indicating the name of the analysis/pipeline form which the quantification comes. Will be used in filenames, plot titles, etc.
#' @param homogenize.mixes logical, whether the two spike-in mixes should be homogenized for the purpose of calculating normalization factors.
#' @param saveResults Logical, whether to save the results of the DEA in the current working directory (default FALSE).
#' @param savePlot Logical, whether to save the plot in the current working directory (default FALSE).
#' @param mix1 A character vector indicating the column names of `dat` that have been spiked with mix 1. If you are using the SEQC data or the dataset at the basis of this package, leave this to NULL.
#'
#' @return A data.frame with the results of the differential expression analysis, as well as a plot.
#'
#' @examples
#' data(exampledata)
#' res <- deSpikein(exampleGeneLevel, method="edgeR", norm="TMM", 
#'   quantification="Tophat-featureCounts")
#'
#' @export
deSpikein <- function(dat, method="edgeR", norm="TMM", quantification="", homogenize.mixes=T, saveResults=F, savePlot=F, mix1=NULL){
  method=match.arg(method, c("t","logt","edgeR","voom","DESeq","DESeq2","EBSeq"))
  if(method %in% c("DESeq","DESeq2","EBSeq")){
    if(!.checkPkg(method))	stop(paste("The package ",method," should first be installed to call this with method='",method,"'",sep=""))
  }
  
  if(!is.null(mix1)){
    if(!all(mix1 %in% colnames(dat)))   stop("The samples with mix1 were provided manually, but cannot be found in the counts matrix!")
    message("Using specified mix1 columns.")
  }else{
    dinf <- .checkDataset(dat)
    if(dinf$dataset == "sim")	stop("This appears to be simulated data, in which spike-ins cannot be analyzed...")
    # indicates which samples have the spike-in mix 1 (others have spike-in mix 2)
    if(dinf$dataset == "w12"){
        mix1 <- c("AJ86","AJ90","AJ91")
    }else{
        if(dinf$dataset == "w6"){
            mix1 <- c("CFG.sh2","c3391S.sh2","w306o.shCTR","w306o.sh2")
        }else{
            mix1 <- paste("A",1:5,sep="")
            dat <- dat[,paste(rep(c("A","B"),each=5),1:5,sep="_")]
        }
    }
  }
  
  data("ercc")
  if(!(sum(as.character(ercc$TX.ID) %in% row.names(dat))+sum(as.character(ercc$ERCC.ID) %in% row.names(dat))>1))    stop("The data given appears not to contain spike-ins, or not having proper row.names.")
  if( sum(as.character(ercc$TX.ID) %in% row.names(dat)) > sum(as.character(ercc$ERCC.ID) %in% row.names(dat)) ){
    # transcript names
    row.names(ercc) <- ercc$TX.ID
  }else{
    # gene names
    row.names(ercc) <- ercc$ERCC.ID
  }
  
  # ensure that all spike-ins are there, set to 0 the missing ones
  sp <- merge(dat,ercc,by="row.names",all.y=T)[,1:(ncol(dat)+1)]
  row.names(sp) <- sp$Row.names
  sp$Row.names <- NULL
  sp <- as.matrix(sp)
  sp[which(is.na(sp) | is.nan(sp) | is.infinite(sp) | sp<0)] <- 0

  # put everything in the same order
  genes <- row.names(sp)[which(row.names(sp) %in% row.names(ercc))]
  sp <- as.data.frame(sp[genes,])
  ercc  <- ercc[genes,]

  d <- as.data.frame(runTests(sp, as.factor(names(sp) %in% mix1), de=method, norm=norm, homogenize.mixes=homogenize.mixes))

  if(saveResults) write.table(d,paste("diff",quantification,method,"tsv",sep="."),row.names=T,col.names=T,sep="\t",quote=F)

  d <- d[genes,]  
  
  name2 <- gsub("%","pc",gsub("(",".",gsub(")",".",quantification,fixed=T),fixed=T),fixed=T)
  if(savePlot)	png(paste("diff",name2,norm,method,"png",sep="."),width=600,height=800)
  deSpikein.compare(d,quantification,norm,method)
  if(savePlot)	dev.off()

  return(d)
}


#' deNanostring
#'
#' Runs a differential expression analysis and compares it to a log-t-test performed on the nanostring data.
#'
#' @param rnaseq Either the (gene-level) count matrix, or a character indicating the location of the file.
#' @param method The differential expression method to use. Either 'edgeR', 'DESeq', 'DESeq2', 'EBSeq', 'voom', 't' (t-test), or 'logt' (t-test on log-transformed values).
#' @param norm The normalization method to use (either "linear", or any of the methods supported by edgeR's \code{\link[edgeR]{calcNormFactors}}. Defaults to 'TMM'.
#' @param quantification A string indicating the quantification that was used to produce the expression matrix (use for labeling)
#' @param threshold The p-value threshold applied on the Nanostring data to identify whether a gene is differentially-expressed or not.
#'
#' @return A deNanostring.compare plot
#'
#' @examples
#' data(exampledata)
#' deNanostring(exampleGeneLevel, method="edgeR", norm="TMM", 
#'   quantification="Tophat-featureCounts")
#'
#' @export
deNanostring <- function(rnaseq=NULL, method="edgeR", norm="TMM", quantification="", threshold=0.01){
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
    cn <- as.logical(c(3,1,1,1,3,1,3,1,1)-1)
    names(cn) <- c("AJ80", "AJ81", "AJ82", "AJ83", "AJ84", "AJ87", "AJ89", "AJ92", "AJ93")
    if(!all(names(cn) %in% colnames(ra)))   stop(paste("There are missing samples, or the column headers are not of the right format. Expected columns:",paste(names(cn),collapse=", ")))
    ra <- ra[,names(cn)]
    res <- runTests(ra, cn, norm=norm, de=method, homogenize.mixes=F)
    deNanostring.compare(res, method, norm, quantification, threshold=threshold)
}

#' deNanostring.compare
#'
#' Compares the p-values and estimated foldchanges from a DEA analysis to a log-t-test performed on the nanostring data.
#'
#' @param results A dataframe as produced by the \code{\link{deNanostring}} function, with gene symbols as row names, and containing the columns "p" (p-value) and "log2FC".
#' @param method The differential expression method used (for reporting only)
#' @param norm The normalization method used (for reporting only)
#' @param quantification The quantification used (for reporting only) 
#' @param threshold The p-value threshold applied on the Nanostring data to identify whether a gene is differentially-expressed or not.
#'
#' @return A layout with 3 plots
#'
#' @export
deNanostring.compare <- function(results, method, norm, quantification="", threshold=0.01){
    if(!all(c("log2FC","p") %in% names(results))) stop("The results dataframe does not appear to be of the correct format. Make sure that the 'log2FC' and 'p' (for p-value) columns are present.")
    data("w12_nanostring")
    if(!all(as.character(nanostring$gene) %in% row.names(results))){
        if(!all(as.character(row.names(nanostring)) %in% row.names(results))){
            stop("The results dataframe does not appear to contain all genes that were measured by nanostring!")
        }else{
            # transcript-level
            na <- nanostring[,3:ncol(nanostring)]
        }
    }else{
        # gene-level
        na <- aggregate(nanostring[,3:ncol(nanostring)],by=list(gene=nanostring$gene),FUN=sum)
        row.names(na) <- na$gene
        na$gene <- NULL
    }
    cn <- c(3,1,1,1,3,1,3,1,1)
    na <- na[,c("AJ80", "AJ81", "AJ82", "AJ83", "AJ84", "AJ87", "AJ89", "AJ92", "AJ93")]
    na$p <- apply(na,1,FUN=function(x){ x <- as.numeric(x); t <- try(t.test(log(x[which(cn==3)]),log(x[which(cn==1)])),silent=T); if(class(t)=="try-error") return(NA); t$p.value})
    na$p[which(is.na(na$p))] <- 1
    na$fc <- apply(na,1,FUN=function(x){ x <- as.numeric(x); mean(x[which(cn==3)])/mean(x[which(cn==1)])})
    results <- results[row.names(na),]
    results$p[which(is.na(results$p))] <- 1
    results$log2FC[which(is.na(results$log2FC))] <- 0
    results$log2FC[which(is.infinite(results$log2FC))] <- 3*sign(results$log2FC[which(is.infinite(results$log2FC))])
    layout(matrix(1:4,nrow=2,byrow=T))
    frame()
    legend("topleft",bty="n",legend=c("Differential expression analysis,","comparison with Nanostring","(7dup vs WBS lines)","",
    paste("Quantification:",quantification),
    paste("Normalization:",norm),
    paste("DEA Method:",method),"",
    paste("Accuracy at FDR<",threshold,": ",round(sum((results$p<threshold)==(na$p<threshold))*100/nrow(na)),"%",sep=""),
    paste("Spearman cor of p-values:",round(cor(results$p, na$p, method="spearman"),3)),
    paste("Pearson cor of log2(FC):",round(cor(results$log2FC, log2(na$fc)),3))
    ))
    plot(-log10(results$p),-log10(na$p),pch=15,col=maketrans("black"),ylab="Nanostring -log10(p-value)",xlab=paste("RNAseq",method,"-log10(p-value)"),main="Comparison of p-values")
    abline(v=-log10(threshold),lty="dashed")
    abline(h=-log10(threshold),lty="dashed")
    legend("bottomright", bty="n", legend=c(sum(results$p < threshold & na$p > threshold)))
    legend("topleft", bty="n", legend=c(sum(results$p > threshold & na$p < threshold)))
    legend("topright", bty="n", legend=c(sum(results$p < threshold & na$p < threshold)))
    ROC(results$p,na$p<threshold,"ROC curve")
    plot(results$log2FC, log2(na$fc), pch=15, col=maketrans("black"),ylab="Nanostring log2(FC)",xlab=paste("RNAseq",method,"log2(FC)"),main="Comparison of foldchanges")
    abline(a=0,b=1,lty="dashed")
    legend("bottomright",bty="n",legend=c(paste("Pearson:",round(cor(results$log2FC, log2(na$fc)),3)),paste("MdAE:",round(MdAE(2^results$log2FC, na$fc),4))))
}

.nfwrapper <- function(x, method="TMM"){
  x <- homomixes(x)
  nf <- calcNormFactors(x, method=method)
  1/(nf/colSums(x)*mean(colSums(x)))
}

#' Runs a spike-in differential expression analysis using sleuth.
#'
#' Runs a differential expression analysis using sleuth and compares it to real differences between spike-in mixes. Currently only for the 12-samples dataset.
#'
#' @param name A string indicating the name of the analysis/pipeline form which the quantification comes. Will be used in filenames, plot titles, etc.
#' @param folder A character vector containing the path to the Salmon/Sailfish/Kallisto quantification folders. If missing, will look for names matching "AJ" in the current working directory.
#' @param norm The normalization method to use (either "linear", or any of the methods supported by edgeR's \code{\link[edgeR]{calcNormFactors}}. Defaults to 'TMM'.
#' @param savePlot Logical, whether to save the plot in the current working directory (default TRUE).
#'
#' @return A data.frame with the results of the differential expression analysis, as well as a plot.
#'
#' @export
sleuthWrapper <- function(name, folders=NULL, norm="TMM", savePlot=F){
  if(!.checkPkg("sleuth"))	stop("The sleuth package must be installed to run this functions")
  library(sleuth)
  if(is.null(folders)){
    folders <- list.files(pattern="AJ")
    if(length(folders)==0){
        folders <- paste(rep(LETTERS[1:2],each=5),rep(1:5,2),sep="")
        if(!all(file.exists(folders))){
            folders <- paste(rep(LETTERS[1:2],each=5),rep(1:5,2),sep="_")
            if(!all(file.exists(folders))){
                stop("Could not find folders; please specify the names!")
            }
        }
    }
  }
  if(!file.exists(paste(folders[[1]],"abundance.h5",sep="/"))){
    message(".h5 files not found... trying to create them (assuming Salmon/Sailfish data)")
    package_sf_as_kal(folders)
  }
  if(!file.exists(paste(folders[[1]],"abundance.h5",sep="/")))	stop("Could not find .h5 files, and could not create them either...")
  design <- data.frame(path=folders, mix=2, stringsAsFactors=F)
  samples <- c("AJ80","AJ81","AJ82","AJ83","AJ84","AJ86","AJ87","AJ89","AJ90","AJ91","AJ92","AJ93")
  if(all(samples %in% folders)){
    row.names(design) = design$sample <- sapply(design$path,FUN=function(x){ strsplit(x,".",fixed=T)[[1]][[1]]})
    design$mix[which(design$sample %in% c("AJ86","AJ90","AJ91"))] <- 1
  }else{
    samples <- paste(rep(LETTERS[1:4],each=5),rep(1:5,4),sep="")
    if(all(gsub("_","",folders,fixed=T) %in% samples)){
        row.names(design) = design$sample <- gsub("_","",folders,fixed=T)
        design <- design[samples[1:10],]
        design$mix[1:5] <- 1
    }else{
        stop("Could not identify the dataset.")
    }
  }
  data(ercc)  
  
  if(is.null(norm)){
    so <- sleuth_prep(design, ~as.factor(mix))
  }else{
    if(norm=="linear"){
      so <- sleuth_prep(design, ~as.factor(mix), norm_fun_counts=function(x){ 1/getLinearNormalizers(homomixes(x)) }, norm_fun_tpm=function(x){ 1/getLinearNormalizers(homomixes(x))})
    }else{
      library(edgeR)
      so <- sleuth_prep(design, ~as.factor(mix), norm_fun_counts=function(x){ .nfwrapper(x, method=norm)}, norm_fun_tpm=function(x){ .nfwrapper(x, method=norm)})
    }
  }
  if(FALSE){
    print("Size factors (counts):")
    print(so$est_counts_sf)
    print("Size factors (TPM):")
    print(so$tpm_sf)
  }
  
  so <- sleuth_fit(so)
  so <- sleuth_wt(so, which_beta = 'as.factor(mix)2')
  res <- sleuth_results(so,'as.factor(mix)2',show_all=T)
  m <- merge(res,ercc[,c("ERCC.ID","TX.ID")],by.x="target_id",by.y="TX.ID",all.y=T)

  row.names(ercc) <- ercc$ERCC.ID
  d <- data.frame(row.names=m$ERCC.ID, log2FC=log2(exp(-1*m$b)), p=m$pval)
  d <- d[as.character(ercc$ERCC.ID),]

  write.table(d,.plfilename("diff",name,norm,"sleuth","tsv"),row.names=T,col.names=T,sep="\t",quote=F)

  if(savePlot)	png(.plfilename("diff",name,norm,"sleuth","png"),width=600,height=800)
  deSpikein.compare(d,name,norm,"sleuth")
  if(savePlot)	dev.off()
}

# t-test wrapper; returns a p-value if possible
.try.t <- function(data, groups, logt=F, var.equal=T){
  if(logt)	data = log2(data+0.001)
  t <- try(t.test(as.numeric(data[which(groups==unique(groups)[[1]])]),as.numeric(data[which(groups==unique(groups)[[2]])]),var.equal=var.equal), silent=T)
  if(class(t)=="htest")	return(as.numeric(t$p.value))
  return(NA)
}



#' Runs a differential expression analysis.
#'
#' Runs a differential expression analysis using the selected method. This is used by the \code{\link{differentialExpression}} function. 
#' This is simply a wrapper to make sure that the different analysis methods use the same normalization method and return results in the same format.
#'
#' @param data The expression matrix or data.frame, with gene symbols or transcript Refseq IDs as row.names.
#' @param groups A logical vector (or coercible to logical) of length ncol(data), indicating to which group each sample (i.e. column in 'data') belongs. There can be only two groups.
#' @param norm The normalization method to use (either "linear", or any of the methods supported by edgeR's \code{\link[edgeR]{calcNormFactors}}. Defaults to 'TMM'.
#' @param de The differential expression method to use. Either 'edgeR', 'DESeq', 'DESeq2', 'EBSeq', 'voom', 't' (t-test), or 'logt' (t-test on log-transformed values). A string indicating the unit of the expression matrix (either "FPKM", "TPM" or "COUNTS").
#' @param homogenize.mixes logical, whether the two spike-in mixes should be homogenized for the purpose of calculating normalization factors.
#'
#' @return A data.frame with the results of the differential expression analysis.
#'
#' @export
runTests <- function(data, groups, norm="TMM", de="edgeR", homogenize.mixes=T){
  de=match.arg(de, c("t","logt","edgeR","voom","DESeq","DESeq2","EBSeq"))
  if(de %in% c("DESeq","DESeq2","EBSeq")){
    if(!.checkPkg(de))	stop(paste("The package ",de," should first be installed to call this with method='",de,"'",sep=""))
  }
  groups <- as.logical(groups)
  library(edgeR)
  if(homogenize.mixes){
    z <- homomixes(data)
  }else{
    z <- data
  }
  nf <- switch(norm,
      "linear"=getLinearNormalizers(z),
      calcNormFactors(z,method=norm))
 if(de=="t" | de=="logt"){
    data <- t(t(data)/nf)
    d <- data.frame(row.names=row.names(data), 
	  log2FC=as.numeric(apply(data,1,FUN=function(x){ getLogFC(x,groups)})),
	  p=as.numeric(apply(data,1,FUN=function(x){ return(as.numeric(.try.t(x,groups,de=="logt",T))) })))
    return(d)
 }
 if(de=="EBSeq"){
  res <- EBTest(Data=as.matrix(data),Conditions=groups,maxround=5,sizeFactors=nf)
  d <- data.frame(row.names=row.names(data), 
	      log2FC=as.numeric(apply(t(t(data)/sf),1,FUN=function(x){ getLogFC(x,groups)})),
	      p=res$PPMatWith0[,1])
  return(d)
 }
 if(de=="edgeR" || de=="voom"){
    y <- DGEList(counts=data,group=as.factor(groups), norm.factors=nf)
    if(de=="edgeR"){
      y <- estimateCommonDisp(y)
      y <- estimateTagwiseDisp(y)
      y <- exactTest(y)
      y <- y$table
      return(data.frame(row.names=row.names(y), log2FC=y$logFC, p=y$PValue))
    }else{
      mm <- model.matrix(~groups)
      v <- voom(y, mm)
      fit <- lmFit(v, mm)
      fit <- eBayes(fit)
      y<-topTable(fit,coef="groupsTRUE",n=nrow(data))
      return(data.frame(row.names=row.names(y), log2FC=y$logFC, p=y$P.Value))
    }
  }
  if(de=="DESeq"){
    library(DESeq)
    cds <- newCountDataSet(round(as.matrix(data)),groups)
    sizeFactors(cds) <- 1/(nf/colSums(data)*mean(colSums(data)))
    cds = estimateDispersions( cds )
    res = as.data.frame(nbinomTest( cds, "TRUE","FALSE" ))
    return(data.frame(row.names=res$id, log2FC=-1*res$log2FoldChange, p=res$pval))
  }
  colData <- data.frame(row.names=names(data), groups=groups)
  if(de=="DESeq2"){
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(round(as.matrix(data)), colData=colData, design = ~ groups)
    sizeFactors(dds) <- 1/(nf/colSums(data)*mean(colSums(data)))
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    res <- as.data.frame(results(dds))
    return(data.frame(row.names=row.names(res), log2FC=res$log2FoldChange, p=res$pvalue))    
  }
}

#' Homogenizes the two spike-in mixes.
#'
#' Divides spike-ins quantifications by their expected foldchange, thereby making all samples in principle equal. This is used for the purpose of calculating normalization factors.
#'
#' @param sp Expression matrix or data.frame, with spike-in IDs as row.names.
#' @param mix1 The column headers of the samples that were spiked with mix1 (optional). This is detected automatically for the datasets included in this study.
#'
#' @return A data.frame including only the spike-ins, and with the two mixes 'homogenized'.
#'
#' @export
homomixes <- function(sp, mix1=NULL){
  data("ercc")
  if(is.null(mix1)){
    if(length(grep("AJ",colnames(sp)))>1){
        mix1 <- c("AJ86","AJ90","AJ91")
    }else{
        if(length(grep("sh2",colnames(sp)))>1){
            mix1 <- c("CFG.sh2","c3391S.sh2","w306o.shCTR","w306o.sh2")
        }else{
            mix1 <- paste("A",1:5,sep=ifelse("A_1" %in% colnames(sp),"_",""))
        }
    }
  }
  if(length(grep("^DQ459",row.names(sp)))>1){
    row.names(ercc) <- ercc$TX.ID
  }else{
    row.names(ercc) <- ercc$ERCC.ID
  }
  sp <- sp[which(row.names(sp) %in% row.names(ercc)), ]  
  ercc <- ercc[row.names(sp),]
  for(i in which(colnames(sp) %in% mix1))	sp[,i] <- as.numeric(sp[,i])/as.numeric(ercc$expected.fold.change.ratio)
  sp
}

#' deSpikein.compare
#'
#' A wrapper to create spike-in DEA benchmarking plots
#'
#' @param d A data.frame with at least the columns "p" (for p-value) and "log2FC", as produced in the \code{\link{differentialExpression}} function.
#' @param quantification Name of the pipeline that produced the underlying quantification, used in plot titles.
#' @param norm Normalization method, used in plot titles.
#' @param method DEA method, used in plot titles.
#'
#' @return Nothing.
#'
#' @export
deSpikein.compare <- function(d,quantification,norm="",method=""){
  data("ercc")
  if("DQ459412" %in% row.names(d)){
    # transcript names
    row.names(ercc) <- ercc$TX.ID
  }else{
    # gene names
    row.names(ercc) <- ercc$ERCC.ID
  }

  if(!all(row.names(ercc) %in% row.names(d)))   message("Some spike-ins are missing. They will be set to p-value 1 and log2FC 0.")
  d <- d[row.names(ercc),]
  row.names(d) <- row.names(ercc)

  d$log2FC[which(is.na(d$log2FC))] <- 0
  d$p[which(is.na(d$p))] <- 1
  if(!("fdr" %in% colnames(d))) d$fdr <- p.adjust(d$p,method="fdr")
  d$p2 <- log10(d$p)
  pceil=-10
  d$p2[which(d$p2 < pceil)] <- pceil
  d$pch <- 20
  d$pch[which(d$p2 <= pceil)] <- 60

  layout(matrix(1:6,nrow=3))
  par(cex=0.8)
  ROC(d$p,ercc$expected.fold.change.ratio!=1)
  myc(d$p2,ercc$expected.fold.change.ratio!=1,title="C: Sensitivity and specificity by p-value")
  fcdev(d,ercc$log2.Mix.1.Mix.2.,title="E: Foldchange deviations")
  hf <- hist(d$fdr,breaks=20,plot=F)
  hp <- hist(d$p,breaks=20,plot=F)
  hist(d$p,breaks=20,xlim=c(0,1),ylim=c(0,max(c(hf$counts,hp$counts))),xlab="",col="grey",main=paste("B: ",quantification,"-",method,"(",norm,")",sep=""))
  points(1:20/20-0.025,c(hf$counts,rep(0,20-length(hf$counts))),lwd=2,col="blue",type="b",pch=20)
  legend("top",fill=c("grey","blue"),legend=c("p-values","q-values"),bty="n")
  myc(log10(d$fdr),ercc$expected.fold.change.ratio!=1,title="D: Sensitivity and specificity by FDR",xlab="log10(FDR)")
  pval2fc(d,ercc$log2.Mix.1.Mix.2.,title="F: P-values by real foldchanges")
}

#' Plots p-values according to real foldchange.
#'
#' @param d A data.frame as produced in the \code{\link{differentialExpression}} function.
#' @param real.log2fc A vector of length nrow(d) containing the real foldchanges.
#' @param jitter.arrows Logical, whether the position of arrows outside the plotting range should be jittered to better see their number.
#' @param title Plot title.
#'
#' @return Nothing.
#'
#' @export
pval2fc <- function(d,real.log2fc,jitter.arrows=T,title=""){
  if(jitter.arrows) real.log2fc[which(d$pch==60)] <- jitter(real.log2fc[which(d$pch==60)])
  plot(d$p2, real.log2fc,ylab="log2(expected foldchange)",xlab="log10(p-value)",pch=d$pch,cex=1.5,col="#00000064",main=title)
  abline(v=log10(0.05),lty="dashed",lwd=2,col="red")
  abline(v=log10(0.01),lty="dashed",lwd=2,col="red")
}

#' Plots deviation from real foldchange by p-value
#'
#' @param d A data.frame as produced in the \code{\link{differentialExpression}} function.
#' @param real.log2fc A vector of length nrow(d) containing the real foldchanges.
#' @param title Plot title.
#'
#' @return Nothing.
#'
#' @export
fcdev <- function(d,real.log2fc,title=""){
  plot(d$p2, d$log2FC-real.log2fc, xlab="log10(p-value)", ylab="Error in log2(foldchange)", pch=d$pch, main=title)
  abline(h=0)
  abline(v=log10(0.05),lty="dashed",col="red",lwd=2)
  abline(v=log10(0.01),lty="dashed",col="red",lwd=2)
}

#' Plots sensitivity and specificity by p-value
#'
#' @param logp The log10(p-value) of each feature.
#' @param sig A logical vector indicating whether each element in logp is actually differentially-expressed (TRUE) or not.
#' @param xlab Plot xlab, defaults to "log10(p-value)".
#'
#' @return Nothing.
#'
#' @export
myc <- function(logp,sig,title="",xlab="log10(p-value)"){
  if(title=="")	title <- "Sensitivity and specificity by p-value"
  x <- c(0,-0.01,-0.05,-0.1,-0.15,-0.25,-0.5,-0.75,-1,-1.25,-1.5,-1.75,(-4:-10)/2)
  sensitivity <- sapply(x,FUN=function(z){ sum(sig & logp <= z, na.rm=T)/sum(sig) })
  specificity <- sapply(x,FUN=function(z){ (sum(!(logp <= z) & !sig , na.rm=T)+sum(is.na(logp) & !sig))/sum(!sig) })
  plot(x, sensitivity,type="l",lwd=2,ylab="",xlab=xlab,col="blue",ylim=c(0,1),main=title)
  lines(x,specificity,lwd=2,col="red")
  legend("topleft",bty="n",border=c(NA,NA,"black","black",NA,NA,NA),fill=c(NA,NA,"blue","red",NA,NA,NA),
	legend=c("","",	paste("Sensitivity (AUC: ",round(auc(x,sensitivity),3),")",sep=""),
			paste("Specificity (AUC: ",round(auc(x,specificity),3),")",sep="")
			))
}

#' compareSpikeinDEcalls
#'
#' Benchmarks a list of differential expression calls between the two spike-in mixes.
#'
#' @param tests A list of tests, each element of which should be a data.frame with a "p" column indicating the p-value under that test, and have ERCC IDs as row.names.
#' @param thres Numeric value between 0 and 1, indicating the p-value treshold to use. Default 0.01.
#' @param colors A vector of colors to be used (R arbitrary colors if NULL), or "greys" to use different tones of grey.
#'
#' @return A layout of two plots, with ROC curves on the left and a barplot of accuracy measurements on the right.
#'
#' @export
compareSpikeinDEcalls <- function(tests, thres=0.01, colors=NULL){
    if(is.null(colors)) colors <- 1:length(tests)
    if(length(colors)==1)   colors <- gray.colors(length(tests))
    data(ercc)
    sig <- ercc$expected.fold.change.ratio!=1
    names(sig) <- ercc$ERCC.ID
    tests <- lapply(tests, FUN=function(x){ x[names(sig),] })
    layout(matrix(1:2,nrow=1))
    library(pheatmap)
    mROC(tests,sig,thres=thres,colors=colors,main="A: ROC curves")
    d <- data.frame(row.names=names(tests))
    d$accuracy <- sapply(tests,FUN=function(x){ sum( (x$p<thres) == sig)/length(sig) })
    d$sensitivity <- sapply(tests,FUN=function(x){ sum( sig & x$p<thres)/sum(sig) })
    d$specificity <- sapply(tests,FUN=function(x){ sum(!(x$p <= thres | sig))/sum(!sig) })
    d$PPV <- sapply(tests,FUN=function(x){ sum(x$p<thres & sig)/(sum(x$p<thres)) })
    barplot(as.matrix(d),beside=T,col=colors,main="B: accuracy at p<0.01")
    legend("topleft",bty="n",fill=colors, legend=names(tests))
}


#' mROC
#'
#' Plots superposed ROC curves.
#'
#' @param tests A list of tests, each element of which should be a data.frame with a "p" column indicating the p-value under that test.
#' @param sig A logical vector indicating whether each row of elements in 'tests' is actually differentially-expressed (TRUE) or not. If this is a named vector, the names will be used to order the rows in 'tests', otherwise they are assumed to be ordered in the same way.
#' @param show.AUC Logical; whether to show the area under the curve in the figure's legend. Default TRUE.
#' @param thres Numeric value between 0 and 1, indicating the p-value treshold to be plotted on each curve. Set to NULL to plot no threshold.
#' @param lwd Line width, passed to the plotting functions. Default 2.
#' @param rounding The rounding (number of digits) of log-transformed p-values at which to plot (default 3). Increasing this number will increase the resolution of the curve, but also increase execution time.
#' @param colors A vector of colors for the different tests. If NULL (default), R basic colors are used.
#' @param ... Any further argument passed to the initial plot function.
#'
#' @export
mROC <- function(tests, sig, show.AUC=T, thres=0.01, lwd=2, rounding=3, na.rm=T, colors=NULL, ...){
  shapes <- c(1,22,23,2,6,7,10)
  if(is.null(colors)) colors <- 1:length(tests)
  ll <- list()
  if(!is.null(names(sig)))  tests <- lapply(tests, FUN=function(x){ x[names(sig),] })
  ff <- data.frame(row.names=names(tests),x=rep(0,length(tests)),y=rep(0,length(tests)))
  for(i in 1:length(tests)){
    de <- names(tests)[i]
    p <- as.numeric(tests[[i]]$p)
    p[which(is.na(p))] <- 1    
    o <- order(p, decreasing=T)
    p <- p[o]
    sig2 <- sig[o]
    pp <- unique(2^round(log2(p+1),rounding)-1)
    d <- data.frame(y=sapply(pp,FUN=function(x){ sum(sig2 & p <= x, na.rm=na.rm)/sum(sig2, na.rm=na.rm) }),
		  x=sapply(pp,FUN=function(x){ sum(!(p <= x | sig2), na.rm=na.rm)/sum(!sig2, na.rm=na.rm) }))
    d$x <- 1-d$x
    if(!is.null(thres)){
        ff[de,] <- c(1-sum(!(p <= thres | sig2), na.rm=na.rm)/sum(!sig2, na.rm=na.rm), sum(sig2 & p <= thres, na.rm=na.rm)/sum(sig2, na.rm=na.rm))
    }
    ll[[de]] <- d[which(!duplicated(d)),]
  }
  plot(ll[[1]]$x, ll[[1]]$y, xlab="1-Specificity",ylab="Sensitivity",xlim=c(0,1),ylim=c(0,1),col=colors[1],type="l",lwd=lwd, ...)
  if(length(tests)>1){
    for(i in 2:length(tests)){
        lines(ll[[i]]$x, ll[[i]]$y, pch=shapes[i],lty=i,col=colors[i],lwd=lwd)
    }
    if(show.AUC){
        leg <- sapply(1:length(tests),FUN=function(x){ paste(names(tests)[x]," (AUC:",round(auc(ll[[x]]$x, ll[[x]]$y),3),")",sep="") })
    }else{
        leg <- names(tests)
    }
    legend("bottomright",bty="n",col=colors,pch=shapes[1:length(tests)],lwd=lwd,lty=1:length(tests),legend=leg)
  }else{
    legend("bottomright",bty="n",legend=paste("AUC:",round(auc(ll[[1]]$x, ll[[1]]$y),3)))
  }
  if(!is.null(thres))   points(ff$x,ff$y,pch=shapes[1:nrow(ff)],col=colors,lwd=lwd)
}

#' ROC
#'
#' Plots the ROC curve.
#'
#' @param p The p-value.
#' @param sig A logical vector indicating whether each element in 'p' is actually differentially-expressed (TRUE) or not.
#' @param rounding The rounding (number of digits) of log-transformed p-values at which to plot (default 3). Increasing this number will increase the resolution of the curve, but also increase execution time.
#' @param main Plot title.
#'
#' @return Nothing.
#'
#' @export
ROC <- function(p,sig,main=NULL,rounding=5){
  if(is.null(main))	main <- "A: ROC curve"
  o <- order(p, decreasing=T)
  p <- as.numeric(p[o])
  sig <- sig[o]
  pp <- unique(2^round(log2(p+1),rounding)-1)
  d <- data.frame(sensitivity=sapply(pp,FUN=function(x){ sum(sig & p <= x)/sum(sig) }),
		  specificity=sapply(pp,FUN=function(x){ sum(!(p <= x | sig))/sum(!sig) }))
  d <- d[which(!duplicated(d)),]
  plot(1-d$specificity, d$sensitivity,xlab="1-Specificity",ylab="Sensitivity",xlim=c(0,1),ylim=c(0,1),type="b",lwd=2,main=main)
  legend("bottomright",bty="n",legend=c(paste("AUC:",round(auc(1-d$specificity, d$sensitivity),3))))
}

#' auc
#'
#' Computes the area under a curve.
#'
#' @param x x values.
#' @param y corresponding y values.
#' @param dens Number of points.
#'
#' @return The area under the curve.
#'
#' @export
auc <- function (x, y, dens = 100){
    if(dens > length(x)) dens <- length(x)
    w <- which(!is.na(x) & !is.na(y))
    y <- y[w]
    x <- x[w]
    o <- order(x)
    x <- x[o]
    y <- y[o]
    idx = 2:length(x)
    x <- as.vector(apply(cbind(x[idx - 1], x[idx]), 1, function(x) seq(x[1], x[2], length.out = dens)))
    y <- as.vector(apply(cbind(y[idx - 1], y[idx]), 1, function(x) seq(x[1], x[2], length.out = dens)))
    idx = 2:length(x)
    as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2
}
