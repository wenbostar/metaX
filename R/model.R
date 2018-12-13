

##' @title Calculate the VIP for PLS-DA
##' @description Calculate the VIP for PLS-DA
##' @rdname calcVIP
##' @docType methods
##' @param x An object of output from \code{plsr}
##' @param ncomp The number of component used in PLS-DA
##' @param ... Additional parameters 
##' @return An vector of VIP value
##' @author Bo Wen \email{wenbostar@@gmail.com}
##' @examples
##' library(pls)
##' x <- matrix(rnorm(1000),nrow = 10,ncol = 100)
##' y <- rep(0:1,5)
##' res <- plsr(y~x)
##' calcVIP(res,2)
setGeneric("calcVIP",function(x,ncomp,...) standardGeneric("calcVIP"))
##' @describeIn calcVIP
setMethod("calcVIP", signature(), function(x,ncomp,...){
    b <- c(x$Yloadings)[1:ncomp]
    TT <- x$scores[,1:ncomp,drop=FALSE]
    SS <- b^2 * colSums(TT^2)
    W <- x$loading.weights[,1:ncomp,drop=FALSE]
    Wnorm2 <- colSums(W^2)
    SSW <- sweep(W^2,2,SS/Wnorm2, "*")
    vips <- sqrt(nrow(SSW) * apply(SSW,1,cumsum) / cumsum(SS))
    if(ncomp >1){
        vip.mat <- as.matrix(t(vips))
    }else{
        vip.mat <- as.matrix(vips)
    }
    return(vip.mat[,ncomp])
    
    
}
)


##' @title runPLSDA
##' @description Validation of the PLS-DA model by using permutation test 
##' statistics
##' @param para An object of metaXpara
##' @param plsdaPara An object of plsDAPara
##' @param auc A logical indicates whether calculate the AUC
##' @param sample Sample class
##' @param valueID The name of column used
##' @param label The label used for plot
##' @param ... additional arguments
##' @return A list object
##' @export
runPLSDA=function(para,plsdaPara,auc=TRUE,sample=NULL,valueID="valueNorm",
                  label="order",...){
    
    ptm <- proc.time()
    
    n <- plsdaPara@nperm
    np <- plsdaPara@ncomp
    t <- plsdaPara@t
    method <- plsdaPara@method
    scale <- plsdaPara@scale
    center <- plsdaPara@center
    validation <- plsdaPara@validation
    k <- plsdaPara@kfold
    cpu <- plsdaPara@cpu
    
    
    outdir <- para@outdir
    prefix <- para@prefix
    result <- list()
    
    ###########################################################################
    ## data preprocess
    ## step 1: Data transformation
    if(t!=0){
        message("log...")
        para <- transformation(para,method=t,valueID=valueID)
        result$log=TRUE
    }else{
        result$log=FALSE
    }
    
    ## step 2: Data scaling
    if(!is.null(scale)){
        message("scale ",scale,"...")
        para <- metaX::preProcess(para = para,scale = scale,center = center,
                                  valueID = valueID)
    }
    
    # sampleList  <- read.delim(para@sampleListFile)
    if(is.null(para@sampleList) || is.na(para@sampleList) ||
       nrow(para@sampleList) ==0){
        sampleList  <- read.delim(para@sampleListFile,stringsAsFactors = FALSE)    
    }else{
        sampleList  <- para@sampleList
    }
    
    peaksData <- para@peaksData
    peaksData <- peaksData[peaksData$class %in% sample,]
    peaksData$class <- as.character(peaksData$class)
    plsData <- dcast(peaksData,sample+class~ID,
                     value.var = valueID)
    plsData$class[is.na(plsData$class)] <- "QC"
    pls.X <- as.matrix(plsData[,-c(1:2)])
    #pls.Y <- as.factor(as.character(plsData$class))
    pls.Y <- plsData$class ## retain the raw data class for bootPLSDA
    
    message("Remove NA value ...")
    numNA <- apply(pls.X,2,function(x){any(is.na(x))})
    message(sum(numNA),"\tfeatures are removed!")
    pls.X <- pls.X[,!numNA]
    
    x <- pls.X
    y <- pls.Y

    y <- as.factor(as.character(y))
    
    result$class <- list()
    
    if(any(is.numeric(y))){
        message("run PLSDA...") 
    }else{
        result$class$rawLabel <- y
        result$class$numericLabel <- as.numeric(as.factor(y))-1
        y <- result$class$numericLabel
    }
    
    ## data preprocess end
    ###########################################################################
    
    ## must put "..." before x,y, so that the parameter from ddply( or sapply)
    ## transfer to fun is omit.
    runInlinePLSDA = function(...,xx,y,ncomp,validation,k=7,
                              method = "oscorespls"){
        sid<-sample(length(y),length(y))
        ## cite "https://github.com/lpantano/isomiRs/blob/master/R/PLSDA.R"
        res <- myPLSDA(xx,y[sid],method=method,ncomp=ncomp,
                       validation=validation,k=k)
        res$cor <- cor(y[sid],y)
        return(unlist(res))    
    }

    p <- myPLSDA(x,y,method=method,ncomp=np,validation=validation,k=k,save=TRUE)
    p$cor <- 1
    result$model <- p$model
    
    p$model <- NULL
    
    result$plsda <- list()
    result$plsda$res <- p
    result$VIP <- calcVIP(result$model,ncomp = np)
    result$x <- x
    result$y <- y
    
    result$ncomp <- np
    result$kfold <- k
    result$validation <- validation
    result$method <- method
    
    dat <- 1:n
    if(cpu==0){
        cpu <- detectCores()
    }
    cl <- makeCluster(getOption("cl.cores", cpu))
    clusterExport(cl, c("runInlinePLSDA"),envir=environment())
    clusterExport(cl, c("myPLSDA"),envir=environment())
    clusterExport(cl, c("plsr"),envir=environment())
    clusterExport(cl, c("R2"),envir=environment())
    clusterExport(cl, c("mvrValstats"),envir=environment())
    
    result$plsda$perm <- parSapply(cl,dat,FUN = runInlinePLSDA,xx=x,y=y,
                                   ncomp=np,
                                   validation=validation,
                                   k=k,method = method)

    stopCluster(cl)
    
    ## pvalue
    result$pvalue <- list()
    result$pvalue$R2 <- sum(result$plsda$perm['R2',] > result$plsda$res$R2)/n
    result$pvalue$Q2 <- sum(result$plsda$perm['Q2',] > result$plsda$res$Q2)/n
    
    message("P-value R2Y: ",result$pvalue$R2)
    message("P-value Q2Y: ",result$pvalue$Q2)
    
    message("R2Y: ",result$plsda$res$R2)
    message("Q2Y: ",result$plsda$res$Q2)
    
    
    ###########################################################################
    
    plotLoading(result$model,fig=paste(outdir,"/",prefix,"-",paste(sample,collapse = "_"),
                                       "-PLSDA-loading.png",sep=""))
    
    fig <- paste(outdir,"/",prefix,"-",paste(sample,collapse = "_"),
                 "-PLSDA.pdf",sep="")
    pngfig <- paste(outdir,"/",prefix,"-",paste(sample,collapse = "_"),
                 "-PLSDA-score.png",sep="")
    
    pdf(fig,width = 5,height = 5)
    
    ## plsda score plot
    if(np >=3){
        
        plotData <- data.frame(x=result$model$scores[,1],
                               y=result$model$scores[,2],
                               z=result$model$scores[,3],
                               sample=plsData$sample,
                               class=result$class$rawLabel)
    }else{
        plotData <- data.frame(x=result$model$scores[,1],
                               y=result$model$scores[,2],
                               sample=plsData$sample,
                               class=result$class$rawLabel)
    }
    plotData$class <- as.character(plotData$class)
    plotData$class[is.na(plotData$class)] <- "QC" 
    sampleList$class <- NULL
    plotData <- merge(plotData,sampleList,by="sample",sort=FALSE)
    ggobj <-ggplot(data = plotData,aes(x=x,y=y,colour=class))+
        geom_hline(yintercept=0,colour="gray")+
        geom_vline(xintercept=0,colour="gray")+
        geom_point()+
        theme_bw()+
        xlab(paste("PC1"," (",sprintf("%.2f%%",explvar(result$model)[1]),") ",
                   sep=""))+
        ylab(paste("PC2"," (",sprintf("%.2f%%",explvar(result$model)[2]),") ",
                   sep=""))+
        #theme_bw()+
        theme(#legend.justification=c(1,1), 
              #legend.position=c(1,1),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())+
              #panel.background=element_rect(fill="#E3E3EE"))+
        #theme(legend.direction = 'horizontal', legend.position = 'top')+
        #stat_ellipse(geom = "polygon", type="euclid",alpha = 0.4, 
        #             aes(fill = class))+
        stat_ellipse(geom = "path",level=0.95)+
        guides(col=guide_legend(nrow=1))+
        theme(legend.position="top")
    
    ggobjpng <- ggobj
    
    if(label == "order"){
        ggobj <- ggobj + geom_text(aes(label=order),size=4,hjust=-0.2)
    }else if(label == "sample"){
        ggobj <- ggobj + geom_text(aes(label=sample),size=4,hjust=-0.2)
    }

    
    print(ggobj)
    
    
    ## 3D PCA plot
    #par(mgp=c(1.6,1,0))
    if(np >=3){
        col <- as.numeric(as.factor(plotData$class))
        s3d <- scatterplot3d(plotData$x,plotData$y,plotData$z,type="h",
                             angle = 24,
                    col.grid="lightblue",lty.hplot=2,pch="",color="gray",
                    xlab = paste("PC1"," (",
                        sprintf("%.2f%%",explvar(result$model)[1]),") ",sep=""),
                    ylab = paste("PC2"," (",
                        sprintf("%.2f%%",explvar(result$model)[2]),") ",sep=""),
                    zlab = paste("PC3"," (",
                        sprintf("%.2f%%",explvar(result$model)[3]),") ",sep="")
        )#color = as.numeric(as.factor(plotData$class)))
        s3d$points(plotData$x,plotData$y,plotData$z, pch = 1,col = col)
        s3d.coords <- s3d$xyz.convert(plotData$x,plotData$y,plotData$z)
        text(s3d.coords$x, s3d.coords$y, labels = plotData$order,
             pos = 4,cex=0.5,col = col)
        
        
        
        classLabel <- levels(as.factor(plotData$class))
        legend(s3d$xyz.convert(max(plotData$x)*0.7, max(plotData$y), 
                               min(plotData$z)), 
               col=as.numeric(as.factor(classLabel)), yjust=0,pch=1,
               legend = classLabel, cex = 0.8)
    }
    
    #dev.off()
    
    ## validation plot
    
    plotdat <- as.matrix(cbind(unlist(result$plsda$res),result$plsda$perm))
    plotdat <- as.data.frame(t(plotdat))
    x1 <- plotdat$cor[1]
    y1 <- plotdat$R2[1]
    y2 <- plotdat$Q2[1]
    plotdat$cor <- abs(plotdat$cor)    
    plotdat <- plotdat[order(plotdat$cor),]
    par(mar=c(3,3,2,1),mgp=c(1.6,0.6,0),cex.lab=1.2,cex.main=0.9)
    plot(plotdat$cor,plotdat$R2,ylim=c(min(plotdat$R2,plotdat$Q2),1),pch=16,
         xlab="Cor",ylab="Value",col="blue")
    points(plotdat$cor,plotdat$Q2,pch=15,col="red")
    
    lm.r <- lm(I(R2-y1)~I(cor-x1)+0,data=plotdat)
    lm.q <- lm(I(Q2-y2)~I(cor-x1)+0,data=plotdat)
    #lines(plotdat$cor,predict(lm.r,data=plotdat$cor),col="blue",lty=2)
    #lines(plotdat$cor,predict(lm.q,data=plotdat$cor),col="red",lty=2)
    int.R <- predict(lm.r,newdata=list(cor=0))+y1
    int.Q <- predict(lm.q,newdata=list(cor=0))+y2
    
    abline(int.R,coef(lm.r),lty=2,col="blue")
    abline(int.Q,coef(lm.q),lty=2,col="red")
    #abline(lm.q,lty=2,col="red")
    legend("bottomright",pch=c(16,15),legend = c("R2","Q2"),col=c("blue","red"))
    
    title(main = paste("Intercepts:","R2=(0.0,",sprintf("%.4f",int.R),
                       "), Q2=(0.0,",sprintf("%.4f",int.Q),")"))
    dev.off()
    
    ## calc AUROC for PLS-DA
    if(auc==TRUE){
        #result$model <- bootPLSDA(pls.X,pls.Y,ncomp = np,sample=sample,... )
    }

    png(filename = pngfig,width = 4,height = 4.2,res = 120,units = "in")    
    print(ggobjpng)
    dev.off()
    
    result$time <- proc.time() - ptm
    
    return(result)
}


##' @title Perform PLS-DA analysis
##' @description Perform PLS-DA analysis
##' @rdname myPLSDA
##' @docType methods
##' @param x A matrix of observations
##' @param y a vector or matrix of responses
##' @param save A logical indicates whether save the \code{plsr} result
##' @param select A logical indicates whether select the best component
##' @param ncomp The number of component used for PLS-DA
##' @param validation See \code{\link{plsr}}
##' @param method See \code{\link{plsr}}
##' @param k k-fold
##' @param ... Additional parameters 
##' @return The PLS-DA result
##' @author Bo Wen \email{wenbostar@@gmail.com}
##' @examples
##' x <- matrix(rnorm(1000),nrow = 10,ncol = 100)
##' y <- rep(0:1,5)
##' res <- myPLSDA(x,y,save=TRUE,ncomp=2,validation="CV",k=7,
##'     method="oscorespls")
setGeneric("myPLSDA",function(x,y,save,select,...) 
    standardGeneric("myPLSDA"))
##' @describeIn myPLSDA
setMethod("myPLSDA", signature(x="ANY",y="ANY",save="logical"),
          function(x,y,ncomp=10,validation="CV",k=7,
                   method = "oscorespls",save=TRUE){
    res <- list()
    plsr.obj <- plsr(y~x,method = method,ncomp=ncomp,validation=validation,
                     segments=k)
    r2res <- pls::R2(plsr.obj,estimate = "all")
    r2q2 <- r2res$val[,1,][,-1]
    res$R2 <- r2q2[1,ncomp]
    res$Q2 <- r2q2[2,ncomp]    
    res$model <- plsr.obj
    return(res)    
})

##' @describeIn myPLSDA
setMethod("myPLSDA", signature(),function(x,y,ncomp=10,validation="CV",k=7,
                                          method = "oscorespls"){
    res <- list()
    plsr.obj <- plsr(y~x,method = method,ncomp=ncomp,validation=validation,
                     segments=k)
    r2res <- pls::R2(plsr.obj,estimate = "all")
    r2q2 <- r2res$val[,1,][,-1]
    res$R2 <- r2q2[1,ncomp]
    res$Q2 <- r2q2[2,ncomp]    
    return(res)    
})

##' @describeIn myPLSDA
setMethod("myPLSDA", signature(x="ANY",y="ANY",select="logical",save="missing"),
          function(x,y,ncomp=10,validation="CV",k=7,method = "oscorespls",
                   select=TRUE){
    res <- list()
    plsr.obj <- plsr(y~x,method = method,ncomp=ncomp,validation=validation,
                     segments=k)
    r2res <- pls::R2(plsr.obj,estimate = "all")
    r2q2 <- r2res$val[,1,][,-1]
    res$R2 <- r2q2[1,]
    res$Q2 <- r2q2[2,]
    res$model <- plsr.obj
    return(res)    
})

## only for bootstrap
## not used now
mybootPLSDA=function(dataset,ind,ncomp=2,validation="CV",k=7,
                     method = "oscorespls"){
    res <- list()
    dataset <- dataset[ind,]
    plsr.obj <- plsr(dataset[,1]~dataset[,-1],method = method,ncomp=ncomp,
                     validation=validation,segments=k)
    r2res <- pls::R2(plsr.obj,estimate = "all")
    r2q2 <- r2res$val[,1,][,-1]
    res$R2 <- r2q2[1,ncomp]
    res$Q2 <- r2q2[2,ncomp]
    #res$model <- plsr.obj
    aa <- c(res$R2,res$Q2)
    names(aa)<-c("R2","Q2")
    #resvip <- calcVIP(x = plsr.obj,ncomp = ncomp)
    return(aa)   
}

##' @title Select the best component for PLS-DA
##' @description Select the best component for PLS-DA
##' @param para A metaXpara object
##' @param np The number of max component
##' @param sample The sample class used
##' @param t Method used to transform the data
##' @param method See \code{\link{plsr}}
##' @param scale Method used to scale the data
##' @param center Centering
##' @param valueID The name of column contained the data
##' @param validation See \code{\link{plsr}}
##' @param k k-fold
##' @param ... Additional parameter
##' @return A list
##' @author Bo Wen \email{wenbostar@@gmail.com}
##' @export
##' @examples
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
##' sampleListFile(para) <- sfile
##' para <- reSetPeaksData(para)
##' para <- missingValueImpute(para)
##' selectBestComponent(para,np=10,sample=c("S","C"),scale="uv",valueID="value")
selectBestComponent=function(para,np=10,sample=NULL,t=1,method = "oscorespls",
                             scale=NULL,center=TRUE,valueID="valueNorm",
                             validation="CV",k=7,...){
    result <- list()
    
    ###########################################################################
    ## data preprocess
    ## step 1: Data transformation
    if(t!=0){
        message("log...")
        para <- transformation(para,method=t,valueID=valueID)
        result$log=TRUE
    }else{
        result$log=FALSE
    }
    
    ## step 2: Data scaling
    if(!is.null(scale)){
        message("scale ",scale)
        para <- preProcess(para = para,scale = scale,center = center,
                           valueID = valueID)
    }
    
    # sampleList  <- read.delim(para@sampleListFile)
    if(is.null(para@sampleList) || is.na(para@sampleList) ||
       nrow(para@sampleList) ==0){
        sampleList  <- read.delim(para@sampleListFile,stringsAsFactors = FALSE)    
    }else{
        sampleList  <- para@sampleList
    }
    
    peaksData <- para@peaksData
    peaksData$class <- as.character(peaksData$class)
    
    peaksData <- peaksData[peaksData$class %in% sample,]
    
    plsData <- dcast(peaksData,sample+class~ID,
                     value.var = valueID)
    plsData$class[is.na(plsData$class)] <- "QC"
    pls.X <- as.matrix(plsData[,-c(1:2)])
    pls.Y <- as.factor(as.character(plsData$class))
    
    message("Remove NA value ...")
    numNA <- apply(pls.X,2,function(x){any(is.na(x))})
    message(sum(numNA),"\tfeatures are removed!")
    pls.X <- pls.X[,!numNA]
    
    x <- pls.X
    y <- pls.Y
    
    if(any(is.numeric(y))){
        message("run PLSDA...") 
    }else{
        result$class$rawLabel <- y
        result$class$numericLabel <- as.numeric(as.factor(y))-1
        y <- result$class$numericLabel
    }
    
    ## data preprocess end
    ###########################################################################
    
    res <- myPLSDA(x,y,method=method,ncomp=np,validation=validation,k=k,
                   select=TRUE)    
    dat <- data.frame(R2=res$R2,Q2=res$Q2,Component=1:np,"R2-Q2"=res$R2-res$Q2,
                      check.names = FALSE)
    #dat$Component <- 1:np
    #names(dat)[1:2] <- c("R2","Q2")
    dat <- melt(dat,value.name = "Value",variable.name = "Metrics",
                id.vars = "Component")
    
    result$res <- res
    result$plotdata <- dat
    
    ###########################################################################
    ## plot
    fig <- paste(para@outdir,"/",para@prefix,"-",
                 paste(sample,collapse = "_"),"-selectLV.pdf",sep="")
    pdf(fig,width = 5,height = 5)
    ggobj <- ggplot(dat,aes(x=Component,y=Value,colour=Metrics))+
        geom_point()+
        geom_line()+
        geom_hline(yintercept=c(0.3,0.5),linetype=2)+
        scale_x_continuous(breaks=1:10)
    print(ggobj)
    dev.off()
    return(result)
}


class2ind=function(x, drop2nd = FALSE) {
    if(!is.factor(x)) stop("'x' should be a factor")
    y <- model.matrix(~ x - 1)
    colnames(y) <- gsub("^x", "", colnames(y))
    attributes(y)$assign <- NULL
    attributes(y)$contrasts <- NULL
    if(length(levels(x)) == 2 & drop2nd) {
        y <- y[,1]
    }
    y
}

##' @title Get a data.frame which contained the peaksData in metaXpara
##' @description Get a data.frame which contained the peaksData in metaXpara
##' @rdname getPeaksTable
##' @docType methods
##' @param para An object of data
##' @param sample Sample class used
##' @param valueID The column name used
##' @return A data.frame
##' @author Bo Wen \email{wenbostar@@gmail.com}
##' @examples
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
##' sampleListFile(para) <- sfile
##' para <- reSetPeaksData(para)
##' res <- getPeaksTable(para)
setGeneric("getPeaksTable",function(para,sample=NULL,valueID="value") 
    standardGeneric("getPeaksTable"))
##' @describeIn getPeaksTable
setMethod("getPeaksTable", signature(para = "metaXpara"), 
          function(para,sample=NULL,valueID="value"){

    message("value = ",valueID)
    
    #sampleList  <- read.delim(para@sampleListFile)
    
    peaksData <- para@peaksData
    if(!is.null(sample)){
        peaksData <- peaksData[peaksData$class %in% sample,]
    }
    peaksData$class <- as.character(peaksData$class)
    pData <- dcast(peaksData,sample+class+batch+order~ID,value.var = valueID)
    pData$class[is.na(pData$class)] <- "QC"
    
    return(pData)
    
}
)


##' @title run OPLS-DA
##' @description Perform the OPLS-DA analysis
##' @param para An object of metaXpara
##' @param oplsdaPara An object of oplsDAPara
##' @param sample Sample class
##' @param valueID The name of column used
##' @param label The label used for plot
##' @param ... additional arguments
##' @return A list object
##' @export
runOPLSDA=function(para,oplsdaPara,auc=TRUE,sample=NULL,valueID="valueNorm",
                  label="order",...){
    
    ptm <- proc.time()
    
    n <- oplsdaPara@nperm
    np <- oplsdaPara@ncomp
    northo <- ifelse(oplsdaPara@northo==0,NA,oplsdaPara@northo)
    t <- oplsdaPara@t
    scale <- oplsdaPara@scale
    center <- oplsdaPara@center
    kfold <- oplsdaPara@kfold
    
    
    outdir <- para@outdir
    prefix <- para@prefix
    result <- list()
    
    ###########################################################################
    ## data preprocess
    ## step 1: Data transformation
    if(t!=0){
        message("log...")
        para <- transformation(para,method=t,valueID=valueID)
        result$log=TRUE
    }else{
        result$log=FALSE
    }
    
    ## step 2: Data scaling
    if(!is.null(scale)){
        message("scale ",scale,"...")
        para <- metaX::preProcess(para = para,scale = scale,center = center,
                                  valueID = valueID)
    }
    
    # sampleList  <- read.delim(para@sampleListFile)
    if(is.null(para@sampleList) || is.na(para@sampleList) ||
       nrow(para@sampleList) ==0){
        sampleList  <- read.delim(para@sampleListFile,stringsAsFactors = FALSE)    
    }else{
        sampleList  <- para@sampleList
    }
    
    peaksData <- para@peaksData
    peaksData <- peaksData[peaksData$class %in% sample,]
    peaksData$class <- as.character(peaksData$class)
    plsData <- dcast(peaksData,sample+class~ID,
                     value.var = valueID)
    plsData$class[is.na(plsData$class)] <- "QC"
    pls.X <- as.matrix(plsData[,-c(1:2)])
    #pls.Y <- as.factor(as.character(plsData$class))
    pls.Y <- plsData$class ## retain the raw data class for bootPLSDA
    
    message("Remove NA value ...")
    numNA <- apply(pls.X,2,function(x){any(is.na(x))})
    message(sum(numNA),"\tfeatures are removed!")
    pls.X <- pls.X[,!numNA]
    
    x <- pls.X
    y <- pls.Y
    
    y <- as.factor(as.character(y))
    
    result$class <- list()
    
    if(any(is.numeric(y))){
        message("run OPLS-DA ...") 
    }else{
        result$class$rawLabel <- y
        #result$class$numericLabel <- as.numeric(as.factor(y))-1
        result$class$numericLabel <- y
        y <- result$class$numericLabel
    }
    
    ## data preprocess end
    ###########################################################################
    
    ## must put "..." before x,y, so that the parameter from ddply( or sapply)
    ## transfer to fun is omit.
    p <- opls(x,y,predI=np,orthoI=northo,crossvalI=kfold,log10L=FALSE,permI=n,
              scaleC="none",...)
    result$model <- p
    
    #result$plsda <- list()
    #result$plsda$res <- p
    #result$VIP <- calcVIP(result$model,ncomp = np)
    result$x <- x
    result$y <- y
    
    result$ncomp <- np
    result$kfold <- kfold
    result$noc <- northo
    
    result$R2Y <- result$model@summaryDF[,"R2Y(cum)"]
    result$Q2 <- result$model@summaryDF[,"Q2(cum)"]
    
    result$perm <- data.frame(cor=result$model@suppLs[["permMN"]][, "sim"],
                              R2Y=result$model@suppLs[["permMN"]][, "R2Y(cum)"],
                              Q2Y=result$model@suppLs[["permMN"]][, "Q2(cum)"])
    
    ## pvalue
    result$pvalue <- list()
    result$pvalue$R2 <- sum(result$perm$R2Y[result$perm$cor<1] > result$R2Y)/n
    result$pvalue$Q2 <- sum(result$perm$Q2Y[result$perm$cor<1] > result$Q2)/n
    
    message("P-value R2Y: ",result$pvalue$R2)
    message("P-value Q2Y: ",result$pvalue$Q2)
    
    message("R2Y: ",result$R2Y)
    message("Q2Y: ",result$Q2)
    
    
    ###########################################################################
    
    #plotLoading(result$model,fig=paste(outdir,"/",prefix,"-",paste(sample,collapse = "_"),
    #                                   "-OPLSDA-loading.png",sep=""))
    
    fig <- paste(outdir,"/",prefix,"-",paste(sample,collapse = "_"),
                 "-OPLSDA.pdf",sep="")
    pdf(fig,width = 3.4,height = 3.4)
    
    ## plsda score plot
    
    plotData <- data.frame(x=result$model@scoreMN[,1],
                           y=result$model@orthoScoreMN[,1],
                           sample=plsData$sample,
                           class=result$class$rawLabel)
    
    plotData$class <- as.character(plotData$class)
    plotData$class[is.na(plotData$class)] <- "QC" 
    sampleList$class <- NULL
    plotData <- merge(plotData,sampleList,by="sample",sort=FALSE)
    ggobj <-ggplot(data = plotData,aes(x=x,y=y,colour=class))+
        geom_hline(yintercept=0,colour="gray")+
        geom_vline(xintercept=0,colour="gray")+
        geom_point()+
        #theme(legend.position="top")+
        #guides(col=guide_legend(nrow=1))+
        theme_bw()+
        xlab(paste("t1"," (",sprintf("%.2f%%",100*result$model@modelDF$R2X[1]),") ",
                              sep=""))+
        ylab("to1")+
        
        #xlab(paste("PC1"," (",sprintf("%.2f%%",explvar(result$model)[1]),") ",
        #           sep=""))+
        #ylab(paste("PC2"," (",sprintf("%.2f%%",explvar(result$model)[2]),") ",
        #           sep=""))+
        #theme_bw()+
        theme(#legend.justification=c(1,1), 
              #legend.position=c(1,1),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())+
              #panel.background=element_rect(fill="#E3E3EE"))+
        #theme(legend.direction = 'horizontal', legend.position = 'top')+
        #stat_ellipse(geom = "polygon", type="euclid",alpha = 0.4, 
        #             aes(fill = class))+
        stat_ellipse(geom = "path",level=0.95)+
        guides(col=guide_legend(nrow=1))+
        theme(legend.position="top")
    
    if(label == "order"){
        ggobj <- ggobj + geom_text(aes(label=order),size=4,hjust=-0.2)
    }else if(label == "sample"){
        ggobj <- ggobj + geom_text(aes(label=sample),size=4,hjust=-0.2)
    }
    
    print(ggobj)
    
    ## validation plot
    
    x1 <- result$model@suppLs[["permMN"]][, "sim"]
    y1 <- result$model@suppLs[["permMN"]][, "R2Y(cum)"]    
    y2 <- result$model@suppLs[["permMN"]][, "Q2(cum)"]
    plotdat <- data.frame(cor=x1,R2=y1,Q2=y2)
    plotdat$cor <- abs(plotdat$cor)    
    plotdat <- plotdat[order(plotdat$cor),]
    par(mar=c(3,3,2,1),mgp=c(1.6,0.6,0),cex.lab=1.2,cex.main=0.9)
    plot(plotdat$cor,plotdat$R2,ylim=c(min(plotdat$R2,plotdat$Q2),1),pch=16,
         xlab="Cor",ylab="Value",col="blue",cex=0.8)
    points(plotdat$cor,plotdat$Q2,pch=15,col="red",cex=0.8)
    
    lm.r <- lm(I(R2-result$R2Y)~I(cor-1)+0,data=plotdat)
    lm.q <- lm(I(Q2-result$Q2)~I(cor-1)+0,data=plotdat)
    #lines(plotdat$cor,predict(lm.r,data=plotdat$cor),col="blue",lty=2)
    #lines(plotdat$cor,predict(lm.q,data=plotdat$cor),col="red",lty=2)
    int.R <- predict(lm.r,newdata=list(cor=0))+y1
    int.Q <- predict(lm.q,newdata=list(cor=0))+y2
    
    abline(int.R,coef(lm.r),lty=2,col="blue")
    abline(int.Q,coef(lm.q),lty=2,col="red")
    #abline(lm.q,lty=2,col="red")
    legend("bottomright",pch=c(16,15),legend = c("R2","Q2"),col=c("blue","red"))
    
    #title(main = paste("Intercepts:","R2=(0.0,",sprintf("%.4f",int.R),
    #                   "), Q2=(0.0,",sprintf("%.4f",int.Q),")"))
    dev.off()
    
    ## calc AUROC for PLS-DA
    if(auc==TRUE){
        #result$model <- bootPLSDA(pls.X,pls.Y,ncomp = np,sample=sample,... )
    }
    
    result$time <- proc.time() - ptm
    
    return(result)
}
