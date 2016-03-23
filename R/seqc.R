#' seqc.diff
#'
#' Runs a comparison of differential expression calls on the SEQC data, comparing calls across conditions to calls within conditions. 
#' The assumption here (and made by the SEQC consortium) is that differential expression across replicates are false positives.
#' See \code{\link{seqc.diff.example}} for an example usage.
#'
#' @param e The count matrix or data.frame. Anything can be used as row names (e.g. gene symbols, transcript ids...). If NULL, the 'site' parameter is used to fetch the quantification from the seqc package (if installed).
#' @param norm The normalization method to use (see \code{\link{donorm}}). Defaults to 'TMM'.
#' @param de The differential expression method to use. Either 'edgeR', 'DESeq', 'DESeq2', 'EBSeq', 'voom', 't' (t-test), or 'logt' (t-test on log-transformed values).
#' @param site If 'e' is NULL, data from the specified sequencing site will be fetched from the seqc package. Default BGI.
#' @param between.groups A vector of integers from 1 to 5 indicating the replicates to be used for comparison across conditions. Default 1:5.
#' @param inner.groups A list of two vectors, each containing the replicates to be used for comparison inside each condition. Default list(c(1,2,3),c(4,5)).
#' @param do.plot Logical, whether to produce Positives plots (see posplot). Default TRUE.
#'
#' @return A list of differential expression calls.
#'
#' @export
seqc.diff <- function(e=NULL, norm="TMM", de="edgeR", site="BGI", between.groups=1:5, inner.groups=list(c(1,2,3),c(4,5)), do.plot=T){
    if(is.null(e)){
        e <- seqc.prepareData(site)
    }else{
        if(class(e)=="character"){
            e <- read.delim(e,header=T,row.names=1)
        }
        colnames(e) <- gsub("_","",colnames(e))
        if(!all(colnames(e) %in% paste(rep(LETTERS[1:4],each=5),rep(1:5,4),sep=""))){
            stop("Could not recognize samples' names. Column headers should be in the format 'A1', 'A2', etc.")
        }
        if(!all(paste(LETTERS[1:4],unique(c(unlist(inner.groups),between.groups)),sep="") %in% colnames(e))){
            stop("Some samples are missing. Check your data, or change the between.groups and inner.groups arguments.")
        }
    }
    gr <- c(rep(T,length(inner.groups[[1]])),rep(F,length(inner.groups[[2]])))
    message("# Running background differential expression analyses")
    fp <- list()
    for(g in c("A","B","C","D")){
        message(gsub("c","",paste(paste(g,inner.groups,sep=""),collapse=" vs ")))
        fp[[g]] <- suppressMessages(runTests(e[,paste(g,unlist(inner.groups),sep="")], gr, norm, de, F))
    }
    message("# Running differential expression analysis A vs B")
    fp[["AvB"]] <- suppressMessages(runTests(e[,paste(rep(c("A","B"),each=length(between.groups)),rep(between.groups,2),sep="")], rep(c(T,F),each=length(between.groups)), norm, de, F))
    message("# Running differential expression analysis C vs D")
    fp[["CvD"]] <- suppressMessages(runTests(e[,paste(rep(c("C","D"),each=length(between.groups)),rep(between.groups,2),sep="")], rep(c(T,F),each=length(between.groups)), norm, de, F))
    if(do.plot){
        message("# Producing plots")
        auc1 <- posplot(fp[["AvB"]]$p,fp[["A"]]$p,fp[["B"]]$p,col="blue")
        auc2 <- posplot(fp[["CvD"]]$p,fp[["C"]]$p,fp[["D"]]$p,col="red",add=T)
        legend("bottomright",bty="n",fill=c("blue","red"),legend=c(paste(c("A vs B (AUC:","C vs D (AUC:"),c(round(auc1,3),round(auc2,3)),")",sep="")))
    }
    return(fp)
}


#' seqc.diff.sleuth
#'
#' Sleuth wrapper for the analysis of SEQC data.
#'
#' @param folders Folders containing the quantification and bootstraps. The folders should be named "A_1" and so on.
#' @param between.groups A vector of integers from 1 to 5 indicating the replicates to be used for comparison across conditions. Default 1:5.
#' @param inner.groups A list of two vectors, each containing the replicates to be used for comparison inside each condition. Default list(c(1,2,3),c(4,5)).
#' @param do.plot Logical, whether to produce Positives plots (see \code{\link{posplot}}). Default TRUE.
#'
#' @return A list of differential expression calls.
#'
#' @export
seqc.diff.sleuth <- function(folders, between.groups=1:5, inner.groups=list(c(1,2,3),c(4,5)), do.plot=T){
    library(sleuth)
    design <- data.frame(path=paste(rep(LETTERS[1:4],each=5),rep(1:5,4),sep="_"), mix=rep(c(1,2,NA,NA),each=5), cond=rep(LETTERS[1:4],each=5), rep=rep(1:5,4), stringsAsFactors=F)
    design$innerGroup <- relevel(as.factor(design$rep %in% inner.groups[[1]]),"FALSE")
    row.names(design) = design$sample = design$path
    fp <- list()
    message("# Running background differential expression analyses")
    for(g in c("A","B","C","D")){
        message(gsub("c","",paste(paste(g,inner.groups,sep=""),collapse=" vs ")))
        d1 <- design[which(design$cond==g),]
        so <- sleuth_prep(d1, ~innerGroup)
        so <- sleuth_wt(sleuth_fit(so), which_beta = 'innerGroupTRUE')
        res <- sleuth_results(so,'innerGroupTRUE')
        fp[[g]] <- data.frame(row.names=row.names(res), log2FC=-1*log2(exp(res$b)), p=res$pval, fdr=p.adjust(res$pval,method="fdr"))
    }
    for(g in list(c("A","B"),c("C","D"))){
        message("")
        message(paste("# Running differential expression analysis",g[[1]],"vs",g[[2]]))
        d1 <- design[which(design$cond %in% g),]
        d1$cond <- relevel(as.factor(d1$cond),g[[1]])
        coe <- paste('cond',g[[2]],sep="")
        so <- sleuth_prep(d1, ~cond)
        so <- sleuth_wt(sleuth_fit(so), which_beta = coe)
        res <- sleuth_results(so,coe)
        fp[[paste(g[[1]],"v",g[[2]],sep="")]] <- data.frame(row.names=row.names(res), log2FC=-1*log2(exp(res$b)), p=res$pval, fdr=p.adjust(res$pval,method="fdr"))
    }
    if(do.plot){
        message("")
        message("# Producing plots")
        auc1 <- posplot(fp[["AvB"]]$p,fp[["A"]]$p,fp[["B"]]$p,col="blue")
        auc2 <- posplot(fp[["CvD"]]$p,fp[["C"]]$p,fp[["D"]]$p,col="red",add=T)
        legend("bottomright",bty="n",fill=c("blue","red"),legend=c(paste(c("A vs B (AUC:","C vs D (AUC:"),c(round(auc1,3),round(auc2,3)),")",sep="")))
    }
    return(fp)    
}


#' Positives plot
#'
#' Produces a 'positives plot', i.e. a plot analogous to the ROC curve, showing the positives across conditions against the positives within condition, across a significance thresholds.
#'
#' @param p A vector of p-values across conditions.
#' @param fp1 A vector of p-values within condition A.
#' @param fp2 An optional vector of p-values within condition B.
#' @param subsamp The maximum number of datapoints to plot. Default 200.
#' @param pround Rounding of the log10 p-values (default 2 digits). Increasing this will increase the precision of the curve, at a speed in cost.
#' @param add Whether to add the data series on top of the current graph. Default FALSE (produces a new plot).
#' @param xlab Passed to the plot function.
#' @param ylab Passed to the plot function.
#' @param main Passed to the plot function.
#' @param col Passed to the plot function, default black.
#' @param xlim Passed to the plot function.
#' @param ylim Passed to the plot function.
#' @param lwd Passed to the plot function, default 3.
#' @param lty Passed to the plot function, default 1.
#' @param auc.plotted Whether to compute the area under the plotted curve (i.e. within plotting limits) instead of that of the full curve (default FALSE).
#'
#' @return Produces a plot and returns the area under the curve.
#'
#' @export
posplot <- function(p,fp1,fp2=NULL,subsamp=200,pround=2,add=F,xlab="DEGs between replicates",ylab="DEGs between conditions",main="Positives plot",col="black",xlim=NULL,ylim=NULL,lwd=3,lty=1, auc.plotted=F){
    if(is.null(fp2)) fp2 <- fp1
    if(is.null(ylim))   ylim <- c(0,length(p))
    if(is.null(xlim))   xlim <- c(0,length(p))
    p <- as.numeric(p)
    p[which(is.na(p))] <- 1
    fp1[which(is.na(fp1))] <- 1
    fp2[which(is.na(fp2))] <- 1
    y <- c(0,10^unique(round(log10(p[order(p)]),pround)),1)
    d <- data.frame("PTP"=sapply(y,FUN=function(x){ sum(p < x) }),
		  "FP"=sapply(y,FUN=function(x){ sum(fp1 < x | fp2 < x) }))
    d <- d[which(!duplicated(d)),]
    d[nrow(d)+1,] <- c(length(p),length(p))
    if(!auc.plotted)    a <- auc(d[,2],d[,1])/length(p)^2
    w <- which(d[,2] < xlim[[2]]*1.05)
    if(auc.plotted)    a <- auc(d[w,2],d[w,1])/(xlim[[2]]*1.05*max(d[w,1]))
    if(max(w) < nrow(d)) w <- c(w,max(w)+1)
    d <- d[w,]
    if(!is.null(subsamp)){
        if(subsamp>2 & subsamp < nrow(d)) d <- d[c(1,1:(subsamp-2)*floor(nrow(d)/(subsamp-2)),nrow(d)),]
    }
    if(add){
        lines(d[,2],d[,1],lwd=3,lty=lty,col=col)
    }else{
        plot(d[,2],d[,1],type="l",lwd=lwd,lty=lty,main=main,xlab=xlab,ylab=ylab,col=col,xlim=xlim,ylim=ylim)
    }
    return(a)
}



#' seqc.diff.example
#'
#' Compares differential expression analysis (DEA) methods using the \code{\link{diff.seqc}} function, and produces ROC-like plots.
#'
#' @param e The count matrix or data.frame. Anything can be used as row names (e.g. gene symbols, transcript ids...). If NULL, the 'site' parameter is used to fetch the quantification from the seqc package (if installed).
#' @param site If 'e' is NULL, data from the specified sequencing site will be fetched from the seqc package. Default BGI.
#' @param tests A vector of the DEA methods to be compared.
#' @param between.groups A vector of integers from 1 to 5 indicating the replicates to be used for comparison across conditions. Default 1:5.
#' @param inner.groups A list of two vectors, each containing the replicates to be used for comparison inside each condition. Default list(c(1,2,3),c(4,5)).
#' @param do.plot Logical, whether to produce Positives plots (see \code{\link{posplot}}). Default TRUE.
#' @param returnData Logical, whether to return the resulting data.
#'
#' @return Produces a plot, and returns a list if returnData=T.
#'
#' @export
seqc.diff.example <- function(e=NULL, site="BGI", tests=c("edgeR","voom","DESeq2"), between.groups=1:5, inner.groups=list(c(1,2,3),c(4,5)), do.plot=T, returnData=F){
    if(is.null(e)){
        e <- seqc.prepareData(site)
    }else{
        if(class(e)=="character"){
            e <- read.delim(e,header=T,row.names=1)
        }
        colnames(e) <- gsub("_","",colnames(e))
    }
    ps <- list()
    for(de in tests){
        message(paste("\n # ",de))
        ps[[de]] <- seqc.diff(e=e,de=de,do.plot=F, between.groups=between.groups, inner.groups=inner.groups)
    }
    
    if(do.plot) seqc.diff.plot(ps)

    if(returnData) return(ps)
}

#' seqc.prepareData
#'
#' Fetches and prepare SEQC data using the \code{\link{seqc}} package.
#'
#' @param site Sequencing site, default BGI.
#' @param removeERCC Logical; whether to remove the spike-ins (default TRUE)
#'
#' @return A data.frame of counts, with gene symbols as row.names.
#'
#' @export
seqc.prepareData <- function(site="BGI", removeERCC=T){
    site <- match.arg(site, c("AGR","BGI","CNL","COH","MAY","NVS"))
    if(!.checkPkg("seqc"))	stop("The package 'seqc' should first be installed before using this function.")
    library(seqc)
    message(paste("# Fetching and preparing SEQC data from the",site,"site..."))
    e <- get(paste("ILM_refseq_gene",site,sep="_"))
    if(removeERCC) e <- e[which(!e$IsERCC),]
    # we aggregate the different lanes for each sample
    sn <- sapply(names(e)[5:ncol(e)],FUN=function(x){ paste(strsplit(x,"_",fixed=T)[[1]][1:2],collapse="")})
    e2 <- aggregate(t(e[,5:ncol(e)]),by=list(sample=sn),FUN=sum)
    row.names(e2) <- e2$sample
    e2$sample <- NULL
    e <- as.data.frame(t(e2))
    e[,which(colnames(e) %in% paste(rep(LETTERS[1:4],each=5),rep(1:5,5),sep=""))]
}

#' seqc.diff.plot
#'
#' Produces a combination of 'positives plot' (see \code{\link{posplot}}) for a list of differential expression calls on the SEQC data.
#' The x/ylim parameters control the 'zoomed plots'; the default settings are good for gene-level using all replicates.
#' For gene-level with only 3 samples/group, use:
#'  xlim=c(0,60), ylimAB=c(5000,15500), ylimCD=c(2000,12000)
#' For transcript-level with all samples, use:
#'  xlim=c(0,250), ylimAB=c(10000,24000), ylimCD=c(4000,17000)
#' For transcript-level with only 3 samples/group, use:
#'  xlim=c(0,200), ylimAB=c(0,23000), ylimCD=c(0,14000)
#'
#' @param ps A list with softwares as names, and lists as elements, as produced by the 'diff.seqc.example' function. Each sublist contains the results of comparisons between  A vs B, C vs D, and within each group.
#' @param xlim x coordinates for the zoomed plot.
#' @param ylimAB y coordinates for the AvsB zoomed plot.
#' @param ylimCD y coordinates for the CvsD zoomed plot.
#'
#' @return Produces 4 plots and returns a list of accuracy values for each test
#'
#' @export
seqc.diff.plot <- function(ps, xlim=c(0,60), ylimAB=c(15000,19000), ylimCD=c(9000,15000)){
    tests <- names(ps)
    layout(matrix(1:4,nrow=2,byrow=T))
    
    message("Producing plots for A vs B comparison")
    auc1 <- list()    
    auc1[[names(ps)[[1]]]] <- posplot(ps[[1]][["AvB"]][["p"]],ps[[1]][["A"]][["p"]],ps[[1]][["B"]][["p"]],col=1,main="A vs B")
    if(length(tests)>1){
        for(i in 2:length(tests)){
            auc1[[names(ps)[[i]]]] <- posplot(ps[[i]][["AvB"]][["p"]],ps[[i]][["A"]][["p"]],ps[[i]][["B"]][["p"]],lty=i,col=i,add=T)
        }
        for(i in 1:length(tests)){
            points(sum(ps[[i]][["A"]][["p"]]<0.01 | ps[[i]][["B"]][["p"]]<0.01, na.rm=T),sum(ps[[i]][["AvB"]][["p"]]<0.01, na.rm=T),col=i,pch=4,cex=2,lwd=2)
        }
    }
    legend("bottomright",bty="n",col=1:length(tests),lwd=3,lty=1:length(tests),legend=sapply(tests,FUN=function(x){ paste(x," (AUC:",round(auc1[[x]],3),")",sep="")}))
    polygon(c(xlim[[1]]-xlim[[2]],2*xlim[[2]],2*xlim[[2]],xlim[[1]]-xlim[[2]],xlim[[1]]-xlim[[2]]),c(ylimAB[[1]],ylimAB[[1]],ylimAB[[2]],ylimAB[[2]],ylimAB[[1]]),lty="dashed")
    auc1[[names(ps)[[i]]]] <- posplot(ps[[1]][["AvB"]][["p"]],ps[[1]][["A"]][["p"]],ps[[1]][["B"]][["p"]],col=1,main="A vs B",pround=3,xlim=xlim,ylim=ylimAB,auc.plotted=T)
    if(length(tests)>1){
        for(i in 2:length(tests)){
            auc1[[names(ps)[[i]]]] <- posplot(ps[[i]][["AvB"]][["p"]],ps[[i]][["A"]][["p"]],ps[[i]][["B"]][["p"]],lty=i,col=i,add=T,pround=3,xlim=xlim,ylim=ylimAB,auc.plotted=T)
        }
        for(i in 1:length(tests)){
            points(sum(ps[[i]][["A"]][["p"]]<0.01 | ps[[i]][["B"]][["p"]]<0.01, na.rm=T),sum(ps[[i]][["AvB"]][["p"]]<0.01, na.rm=T),col=i,pch=4,cex=2,lwd=2)
        }
    }
    legend("bottomright",bty="n",col=1:length(tests),lwd=3,lty=1:length(tests),legend=sapply(tests,FUN=function(x){ paste(x," (AUC:",round(auc1[[x]],3),")",sep="")}))

    message("Producing plots for C vs D comparison")
    auc2 <- list()
    auc2[[names(ps)[[1]]]] <- posplot(ps[[1]][["CvD"]][["p"]],ps[[1]][["C"]][["p"]],ps[[1]][["D"]][["p"]],col=1,main="C vs D")
    if(length(tests)>1){
        for(i in 2:length(tests)){
            auc2[[names(ps)[[i]]]] <- posplot(ps[[i]][["CvD"]][["p"]],ps[[i]][["C"]][["p"]],ps[[i]][["D"]][["p"]],lty=i,col=i,add=T)
        }
        for(i in 1:length(tests)){
            points(sum(ps[[i]][["C"]][["p"]]<0.01 | ps[[i]][["D"]][["p"]]<0.01, na.rm=T),sum(ps[[i]][["CvD"]][["p"]]<0.01, na.rm=T),col=i,pch=4,cex=2,lwd=2)
        }
    }
    legend("bottomright",bty="n",col=1:length(tests),lwd=3,lty=1:length(tests),legend=sapply(tests,FUN=function(x){ paste(x," (AUC:",round(auc2[[x]],3),")",sep="")}))
    polygon(c(xlim[[1]]-xlim[[2]],2*xlim[[2]],2*xlim[[2]],xlim[[1]]-xlim[[2]],xlim[[1]]-xlim[[2]]),c(ylimCD[[1]],ylimCD[[1]],ylimCD[[2]],ylimCD[[2]],ylimCD[[1]]),lty="dashed")
    auc2[[names(ps)[[1]]]] <- posplot(ps[[1]][["CvD"]][["p"]],ps[[1]][["C"]][["p"]],ps[[1]][["D"]][["p"]],col=1,main="C vs D",pround=3,xlim=xlim,ylim=ylimCD,auc.plotted=T)
    if(length(tests)>1){
        for(i in 2:length(tests)){
            auc2[[names(ps)[[i]]]] <- posplot(ps[[i]][["CvD"]][["p"]],ps[[i]][["C"]][["p"]],ps[[i]][["D"]][["p"]],lty=i,col=i,add=T,pround=3,xlim=xlim,ylim=ylimCD,auc.plotted=T)
        }
        for(i in 1:length(tests)){
            points(sum(ps[[i]][["C"]][["p"]]<0.01 | ps[[i]][["D"]][["p"]]<0.01, na.rm=T),sum(ps[[i]][["CvD"]][["p"]]<0.01, na.rm=T),col=i,pch=4,cex=2,lwd=2)
        }
    }
    legend("bottomright",bty="n",col=1:length(tests),lwd=3,lty=1:length(tests),legend=sapply(tests,FUN=function(x){ paste(x," (AUC:",round(auc2[[x]],3),")",sep="")}))
    v <- rep(0,length(tests))
    s <- data.frame(row.names=tests, AB.FN=v, AB.TP=v, AB.FDR=v, CD.FN=v, CD.TP=v, CD.FDR=v)
    for(de in tests){
        s[de,1] <- sum(ps[[de]][["A"]]$p<0.01, na.rm=T)+sum(ps[[de]][["B"]]$p<0.01, na.rm=T)
        s[de,2] <- sum(ps[[de]][["AvB"]]$p<0.01, na.rm=T)
        s[de,3] <- s[de,1]/s[de,2]
        s[de,4] <- sum(ps[[de]][["C"]]$p<0.01, na.rm=T)+sum(ps[[de]][["D"]]$p<0.01, na.rm=T)
        s[de,5] <- sum(ps[[de]][["CvD"]]$p<0.01, na.rm=T)
        s[de,6] <- s[de,4]/s[de,5]
    }
}