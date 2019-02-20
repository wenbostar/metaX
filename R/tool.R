
calcRatio=function(x=NULL,class=NULL,weights=NULL,
                   method=c("mean","median","weight"),
                   case=NA,control=NA, ...){
    
    r1 <- x[class==case]
    r2 <- x[class==control]
    
    r1 <- r1[!is.nan(r1)]
    r1 <- r1[!is.infinite(r1)]
    r2 <- r2[!is.nan(r2)]
    r2 <- r2[!is.infinite(r2)]
    
    if(method=="mean"){
        ratio <- mean(r1,na.rm = TRUE)/mean(r2,na.rm = TRUE)
    }else if(method=="median"){
        ratio <- median(r1,na.rm = TRUE)/median(r2,na.rm = TRUE)
    }else{
        stop("method must be in 'mean', 'median' and 'weight'!")
    }
    
    
    
}




##Hotellings T2 statistic and ellipse calculation function##
HotE = function(x, y, len = 200,alfa=0.95){
    N <- length(x)
    A <- 2
    mypi <- seq(0, 2 * pi, length = len)
    r1 <- sqrt(var(x) * qf(alfa,2, N-2) * (2*(N^2-1)/(N*(N - 2))))
    r2 <- sqrt(var(y) * qf(alfa,2, N-2) * (2*(N^2-1)/(N *(N-2))))
    cbind(r1 * cos(mypi) + mean(x), r2*sin(mypi)+mean(y))
}

##' @title permutePLSDA
##' @description Validation of the PLS-DA model by using permutation test 
##' statistics
##' @param x a matrix of observations.
##' @param y a vector or matrix of responses.
##' @param n number of permutations to compute the PLD-DA p-value based on R2 
##' magnitude. Default n=100
##' @param np the number of components to be used in the modelling.
##' @param outdir output dir
##' @param prefix the prefix of output figure file 
##' @param tol tolerance value based on maximum change of cumulative R-squared 
##' coefficient for each additional PLS component. Default tol=0.001
##' @param cpu 0
##' @param ... additional arguments
##' @return pvalue
##' @export
permutePLSDA=function(x,y,n=100,np=2,outdir = "./", prefix="metaX",tol=0.001,
                      cpu=0,...){
    
    p <- plsDA(variables = x,group = y,autosel = FALSE,comps = np)
    r2 <- NA
    if(p$R2[dim(p$R2)[1],3] <= tol){
        a.perm <- which(p$R2[,3] < tol)
        if(a.perm[1] > 1){ 
            r2 <- p$R2[a.perm[1]-1,4] 
        }else if(a.perm[1] == 1){
            r2 <- p$R2[1,4] 
        }
    }else{
        r2 <- p$R2[dim(p$R2)[1],4]
    }
    
    
    
    fig <- paste(outdir,"/",prefix,"-permutation_PLSDA.pdf",sep="")
    pdf(fig,width = 4,height = 4)
    
    
    dat <- 1:n
    
    if(cpu==0){
        cpu <- detectCores()
    }
    cl <- makeCluster(getOption("cl.cores", cpu))
    clusterExport(cl, c("run_PLSDA"),envir=environment())

    res<-parSapply(cl,dat,FUN = run_PLSDA,xx=x,y=y,np=np,tol=tol)
    stopCluster(cl)


    ## sometimes, maybe there are rows that the values are nan.    
    trueTotal <- sum(!is.na(res))
    message("True permutation number is ",trueTotal)
    pvalue <- sum(r2 <= res,na.rm = TRUE)/trueTotal
    message("p-value = ",pvalue)
    plotDat <- data.frame(R2=res)

    
    ggobj <- ggplot(data = plotDat, aes(x=R2))+
        geom_histogram(colour="white")+
        geom_vline(xintercept = r2,colour="red",size=1)+
        ylab("Count")+
        scale_y_continuous(expand = c(0,0))+
        theme(legend.justification=c(1,1), 
              legend.position=c(1,1),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background=element_rect(fill="#E3E3EE"))
    print(ggobj)
    dev.off()
    return(pvalue)
}

run_PLSDA = function(...,xx,y,np,tol=0.001){
    ## must put "..." before x,y, so that the parameter from ddply( or sapply)
    ## transfer to fun is omit.
    sid<-sample(length(y),length(y))
    p <- plsDA(variables = xx,group = y[sid],autosel = FALSE,comps = np)
    R2value <- NA
    if(p$R2[dim(p$R2)[1],3] <= tol){
        a.perm <- which(p$R2[,3] < tol)
        if(a.perm[1] > 1){ 
            R2value <- p$R2[a.perm[1]-1,4] 
        }else if(a.perm[1] == 1){
            R2value <- p$R2[1,4] 
        }
    }else{
        R2value <- p$R2[dim(p$R2)[1],4]
    }
    R2value
}


countMissingValue = function(x,ratio,omit.negative=TRUE){
    if(omit.negative){
        rp <- (sum(x<=0 | is.na(x))/length(x))>=ratio
    }else{
        rp <- ( sum(is.na(x)) / length(x) )>=ratio
    }
    return(rp)
}

##' @title metaXpipe
##' @description metaXpipe
##' @rdname metaXpipe
##' @docType methods
##' @param para A metaXpara object.
##' @param plsdaPara A plsDAPara object.
##' @param cvFilter Filter peaks which cv > cvFilter in QC samples.
##' @param remveOutlier Remove outlier samples.
##' @param outTol The threshold to remove outlier samples.
##' @param doQA Boolean, setting the argument to TRUE will perform plot 
##' quality figures. 
##' @param doROC A logical indicates whether to calculate the ROC
##' @param qcsc QC-based batch correction, 0=none,1=QC-RLSC(Quality 
##' control-robust loess signal correction),2=SVR(SVR normalization).
##' @param nor.method Normalization method.
##' @param pclean Boolean, setting the argument to TRUE to perform 
##' data cleaning 
##' @param t Data transformation method. See \code{\link{transformation}}.
##' @param scale Data scaling method.
##' @param idres A file containing the metabolite identification result
##' @param nor.order The order of normalization, only valid when \code{qcsc} is 
##' TRUE. 1: before QC-based batch correction, 2: after QC-based batch correction.
##' @param out.rmqc Boolean, setting the argument to TRUE to remove the QC 
##' samples for the csv file.
##' @param saveRds Boolean, setting the argument to TRUE to save some objects to
##' disk for debug. Only useful for developer. Default is TRUE.
##' @param cpu The number of cpu used, default is all available cpus.
##' @param missValueRatioQC The cutoff of the ratio of miss value for the 
##' features in QC, default is 0.5.
##' @param missValueRatioSample The cutoff of the ratio of miss value for the 
##' features in sample, default is 0.8.
##' @param pcaLabel The label used for PCA score plot, "none","order" and 
##' "sample" are supported. Default is "none"
##' @param ... Other argument
##' @return A metaXpara object.
##' @exportMethod
##' @author Bo Wen \email{wenbostar@@gmail.com}
##' @examples
##' \dontrun{
##' ## example 1: no QC sample
##' library(faahKO)
##' xset <- group(faahko)
##' xset <- retcor(xset)
##' xset <- group(xset)
##' xset <- fillPeaks(xset)
##' peaksData <- as.data.frame(groupval(xset,"medret",value="into"))
##' peaksData$name <- row.names(peaksData)
##' para <- new("metaXpara")
##' rawPeaks(para) <- peaksData
##' ratioPairs(para) <- "KO:WT"
##' outdir(para) <- "test"
##' sampleListFile(para) <- system.file("extdata/faahKO_sampleList.txt", 
##'     package = "metaX")
##' plsdaPara <- new("plsDAPara")
##' p <- metaXpipe(para,plsdaPara=plsdaPara)
##' 
##' ## example 2: has QC samples
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
##' sampleListFile(para) <- sfile
##' ratioPairs(para) <- "S:C"
##' plsdaPara <- new("plsDAPara")
##' p <- metaXpipe(para,plsdaPara=plsdaPara)
##' }
setGeneric("metaXpipe",function(para,plsdaPara,cvFilter = 0.3,remveOutlier=TRUE,
                              outTol=1.2,doQA=TRUE,doROC=TRUE,qcsc=0,nor.method="pqn",
                              pclean=TRUE,t=1,scale="uv",idres=NULL,center=TRUE,
                              nor.order=1,out.rmqc=FALSE,saveRds=TRUE,cpu=0,
                              missValueRatioQC=0.5,missValueRatioSample=0.8,
                              pcaLabel="none",classCol=NULL,
                              fig_width_boxplot = 9.5,...) 
    standardGeneric("metaXpipe"))
##' @describeIn metaXpipe
setMethod("metaXpipe", signature(para = "metaXpara"),
          function(para,plsdaPara,cvFilter = 0.3,remveOutlier=TRUE,outTol=1.2,
                   doQA=TRUE,doROC=TRUE,qcsc=0,nor.method="pqn",
                   pclean=TRUE,t=1,scale="uv",idres=NULL,center=TRUE,nor.order=1,
                   out.rmqc=FALSE,saveRds=TRUE,cpu=0,
                   missValueRatioQC=0.5,missValueRatioSample=0.8,
                   pcaLabel="none",classCol=NULL,
                   fig_width_boxplot = 9.5,...){
    
    if(is.null(para@sampleList) || is.na(para@sampleList) ||
        nrow(para@sampleList) ==0){
        sampleList  <- read.delim(para@sampleListFile,stringsAsFactors = FALSE)    
    }else{
        sampleList  <- para@sampleList
    }
              
    checkSampleList(sampleList)
    
    check_input_data_format(para@rawPeaks)
              
    if(is.null(para@ratioPairs)){
        stop("Please set the value of ratioPairs!")
    }
    
    ## note
    plsdaPara@scale <- scale
    plsdaPara@t <- t
    
              
    ## directory structure, data (all figures and other files), report.html
    raw_outdir <- para@outdir
    ## all figures and other files
    para@outdir <- paste(para@outdir,"/data",sep="") 
    makeDirectory(para)
    
    
    ## create report 
    report<-newCustomReport("metaX Report")    
    report<-setReportTitle(report,
                           "Metabolomic Data Analysis with metaX")
    
    ##The first section
    ## ref : http://www.mdpi.com/2218-1989/3/3/552/pdf
    s1 <-newSection("Introduction")
    s1 <-addTo(s1,newParagraph(
        "Metabolomics is a growing and powerful technology capable of detecting 
        hundreds to thousands of metabolites in tissues and biofluids. 
        metabolomics alterations represent changes in the phenotype and 
        molecular physiology."))
    
    s2<-newSection("Methods and Data")
    ##Introduce the experiment design and data analysis
    s2_sub1 <- newSubSection(asStrong("Summary of Data Set"))
    
    # expDesign <- read.delim(para@sampleListFile)
    
    if(is.null(para@sampleList) || is.na(para@sampleList) ||
       nrow(para@sampleList) ==0){
        expDesign  <- read.delim(para@sampleListFile,stringsAsFactors = FALSE)    
    }else{
        expDesign  <- para@sampleList
    }
    
    expDesign$class <- as.character(expDesign$class)
    expDesign$class[is.na(expDesign$class)] <- "QC"
    
    
    expClass <- ddply(expDesign,.(class),summarise,n=length(sample))
    expClassTB <- newTable(expClass,"Sample group information.")
    
    expBatch <- ddply(expDesign,.(batch),summarise,n=length(sample))
    expBatchTB <- newTable(expBatch,"Experiment batch information.")
    
    expTB <- newTable(expDesign,"Sample information.")
    s2_sub1 <- addTo(s2_sub1,expClassTB,expBatchTB,expTB)
    
    
    s2 <- addTo(s2,s2_sub1)
    
    s3<-newSection("Results")
    ##The first subsection
    s3_sub1 <- newSubSection(asStrong("Data quality analysis"))
    s3_sub1 <- addTo(s3_sub1,
                     newParagraph("This part contains the basic data quality."))
              
    ## 1. peak picking
    ## whether para has contained the rawPeaks
    if(ncol(para@rawPeaks)<=0){
        para <- peakFinder(para)
        if(saveRds){
            save(para,file = paste(para@outdir,"/",para@prefix,"-peakFinder.rda",
                                   sep=""))
        }
    }
    message(date(),"\treset the peaksData...")
    para <- reSetPeaksData(para = para)
    makeMetaboAnalystInput(para,rmQC=out.rmqc,valueID="value",prefix="raw")
    
    ## 2. pre-processing, filter noise peaks
    if(hasQC(para)){
        message(date(),"\tfilter peaks in QC sample:>=",missValueRatioQC)
        pre_para <- para
        para <- filterQCPeaks(para = para,ratio = missValueRatioQC)
        rpeak <- length(unique(pre_para@peaksData$ID)) - length(unique(para@peaksData$ID))
        
        s3_sub1 <- addTo(s3_sub1,
                         newParagraph(paste("Remove peaks in QC samples with missing 
            value greater than ",100*missValueRatioQC, " percent:",sep=""),rpeak,"."))
        rm(pre_para)
        
    }
    message(date(),"\tfilter peaks in non-QC sample:>=",missValueRatioSample)
    pre_para <- para
    para <- filterPeaks(para = para,ratio = missValueRatioSample)
    rpeak <- length(unique(pre_para@peaksData$ID)) - length(unique(para@peaksData$ID))
    
    s3_sub1 <- addTo(s3_sub1,
                     newParagraph(paste("Remove peaks in non-QC samples with missing 
            value greater than ",100*missValueRatioSample," percent:",sep=""),rpeak,"."))
    rm(pre_para)
    
    message("The CV distribution before normalization: ")
    printCV(para,valueID = "value")
    
    ## 3. Quality control
    if(doQA==TRUE){
        message("plot peak number distribution...")
        fig <- plotPeakNumber(para)
        fig1 <- paste("data/",basename(fig$fig),sep="")
        fig2 <- paste("data/",basename(fig$highfig),sep="")
        s3_sub1 <- addTo(s3_sub1,
                         newFigure(fig1,"Peak number distribution.",
                                   fileHighRes = fig2))
        
        
        message("plot CV distribution...")
        fig <- metaX::plotCV(para)
        fig1 <- paste("data/",basename(fig$fig),sep="")
        fig2 <- paste("data/",basename(fig$highfig),sep="")
        s3_sub1 <- addTo(s3_sub1,
                         newFigure(fig1,"Peak CV distribution.",
                                   fileHighRes = fig2))
        s3_sub1 <- addTo(s3_sub1,
                         newTable(fig$stat,"CV stat."))
        
        
        message("plot missing value distribution...")
        fig <- plotMissValue(para,height = 3.7,width = 9.15,...)
        fig1 <- paste("data/",basename(fig$fig),sep="")
        fig2 <- paste("data/",basename(fig$highfig),sep="")
        s3_sub1 <- addTo(s3_sub1,
                         newFigure(fig1,"Missing value distribution.",
                                   fileHighRes = fig2))
        
        
        message("plot peak intensity distribution...")
        fig <- plotIntDistr(para,width = fig_width_boxplot)
        fig1 <- paste("data/",basename(fig$fig),sep="")
        fig2 <- paste("data/",basename(fig$highfig),sep="")
        s3_sub1 <- addTo(s3_sub1,
                         newFigure(fig1,"Peak intensity distribution.",
                                   fileHighRes = fig2))
        
        
                    
        
        if(hasQC(para)){
            message("plot correlation heatmap...")
            #saveRDS(para,file = "para.rda")
            fig <- plotCorHeatmap(para = para,valueID = "value",anno = TRUE,
                                  samples = NA, ## QC sample
                                  height = 6,width = 6,cluster = FALSE,
                                  classCol=classCol)
            fig1 <- paste("data/",basename(fig$fig),sep="")
            fig2 <- paste("data/",basename(fig$highfig),sep="")
            s3_sub1 <- addTo(s3_sub1,
                             newFigure(fig1,"Correlation heatmap.",
                                       fileHighRes = fig2))
            
        }
        message("plot TIC distribution...")
        fig <- plotPeakSumDist(para,valueID = "value",height = 3.7,width = 5.5)
        fig1 <- paste("data/",basename(fig$fig),sep="")
        fig2 <- paste("data/",basename(fig$highfig),sep="")
        s3_sub1 <- addTo(s3_sub1,
                         newFigure(fig1,"TIC distribution.",
                                   fileHighRes = fig2))
        
        ticTable <- para@peaksData %>% group_by(sample,batch,class) %>% 
            dplyr::summarise(tic=sum(value)) %>% 
            group_by(batch) %>% 
            dplyr::summarise(n00=quantile(tic)[1],
                      n25=quantile(tic)[2],
                      n50=quantile(tic)[3],
                      n75=quantile(tic)[4],
                      n100=quantile(tic)[5])
        print(ticTable)
        
        message("plot average intensity distribution...")
        fig <- plotPeakMeanDist(para,valueID = "value")
        fig1 <- paste("data/",basename(fig$fig),sep="")
        fig2 <- paste("data/",basename(fig$highfig),sep="")
        s3_sub1 <- addTo(s3_sub1,
                         newFigure(fig1,"Average intensity distribution.",
                                   fileHighRes = fig2))
        
        meanTable <- para@peaksData %>% group_by(sample,batch,class) %>% 
            dplyr::summarise(meanIntensity=mean(value,na.rm = TRUE)) %>% 
            group_by(batch) %>% 
            dplyr::summarise(n00=quantile(meanIntensity)[1],
                             n25=quantile(meanIntensity)[2],
                             n50=quantile(meanIntensity)[3],
                             n75=quantile(meanIntensity)[4],
                             n100=quantile(meanIntensity)[5])
        print(meanTable)

    }
    
    message(date(),"\tmissing value inputation...")
    para <- missingValueImpute(para,cpu=cpu,...)
    s3_sub1 <- addTo(s3_sub1,newParagraph("The missing value were imputed by ",
                                          para@missValueImputeMethod,"."))
    ## 4. remove outlier samples
    ## outlier sample remove
    if(remveOutlier){
        message(date(),"\tauto remove outlier sample...")
        pr <- metaX::normalize(para = para, method = nor.method,
                                 valueID="value",...)
        ## maybe have missing value
        if(nor.method=="combat"){
            message(date(),"\tDo missing value imputation after combat normalization...")
            pr <- missingValueImpute(x = pr,valueID="value",cpu=cpu)
        }
        pr <- transformation(pr,method = t,valueID = "value")
        pr <- metaX::preProcess(pr,scale = scale,center = center,
                                valueID = "value")
        ## 
        if(saveRds){
            saveRDS(para,file = paste(para@outdir,"/",para@prefix,"-para_before_removeOutlier.rds",
                                      sep=""))
        }
        removeSampleNames <- autoRemoveOutlier(pr,outTol=outTol,scale="none",
                                               center=FALSE,valueID = "value")
        para <- removeSample(para,rsamples = removeSampleNames)
        #para@peaksData <- para@peaksData[!para@peaksData$sample %in% 
        #                                     removeSampleNames,]
        
        if(length(removeSampleNames)>=1){
            s3_sub1 <- addTo(s3_sub1,newParagraph("Remove outlier sample."))
            reTB <- dplyr::filter(expDesign,sample %in% removeSampleNames)
            s3_sub1 <- addTo(s3_sub1,newTable(reTB,"The sample removed."))
        }else{
            s3_sub1 <- addTo(s3_sub1,newParagraph("Don't find outlier sample."))
        }
    }
    

    
    if(hasQC(para) && qcsc!=0){
        
        if(nor.order==1){
            if(qcsc==1){
                message(date(),"\tnormalization before QC-RLSC...")
                para <- metaX::normalize(para = para, method = nor.method,
                                         valueID="value",...)
                
                ## maybe have missing value
                if(nor.method=="combat"){
                    message(date(),"\tDo missing value imputation after combat normalization...")
                    para <- missingValueImpute(x = para,valueID="value",cpu=cpu)
                }
                
                message(date(),"\tQC-RLSC...")
                res <- doQCRLSC(para = para)
                para <- res$metaXpara
                message("The CV distribution after normalization: ")
                printCV(para,valueID = "valueNorm")
                prefix_tmp <- para@prefix
                para@prefix <- paste(para@prefix,"-norm",sep="")
                metaX::plotCV(para,valueID="valueNorm")
                para@prefix <- prefix_tmp
                para <- filterQCPeaksByCV(para,cvFilter = cvFilter,
                                          valueID = "valueNorm")
                fig <- plotQCRLSC(para = para)
            }else if(qcsc==2){
                message(date(),"\tnormalization before SVR-based batch correction...")
                para <- metaX::normalize(para = para, method = nor.method,
                                         valueID="value",...)
                ## maybe have missing value
                if(nor.method=="combat"){
                    message(date(),"\tDo missing value imputation after combat normalization...")
                    para <- missingValueImpute(x = para,valueID="value",cpu=cpu)
                }
                message(date(),"\tSVR-based batch correction...")
                res <- svrNormalize(para = para,ntop=5,cpu=cpu)
                para <- res$metaXpara
                message("The CV distribution after normalization: ")
                printCV(para,valueID = "valueNorm")
                prefix_tmp <- para@prefix
                para@prefix <- paste(para@prefix,"-norm",sep="")
                metaX::plotCV(para,valueID="valueNorm")
                para@prefix <- prefix_tmp
                para <- filterQCPeaksByCV(para,cvFilter = cvFilter,
                                          valueID = "valueNorm")
                
                fig <- plotQCRLSC(para = para)
            }else if(qcsc==3){
                message(date(),"\tnormalization before ComBat batch correction...")
                para <- metaX::normalize(para = para, method = nor.method,
                                         valueID="value",...)
                ## maybe have missing value
                #if(nor.method=="combat"){
                #    message(date(),"\tDo missing value imputation after combat normalization...")
                #    para <- missingValueImpute(x = para,valueID="value",cpu=cpu)
                #}
                message(date(),"\tComBat batch correction...")
                para <- batchCorrect(para,valueID = "value",use_class = FALSE,
                                     cpu=cpu,
                                     impute_method=para@missValueImputeMethod)
                para@peaksData$valueNorm <- para@peaksData$value
                message("The CV distribution after normalization: ")
                printCV(para,valueID = "valueNorm")
                prefix_tmp <- para@prefix
                para@prefix <- paste(para@prefix,"-norm",sep="")
                metaX::plotCV(para,valueID="valueNorm")
                para@prefix <- prefix_tmp
                para <- filterQCPeaksByCV(para,cvFilter = cvFilter,
                                          valueID = "valueNorm")
                
            }else{
                stop("Not valid qcsc value!")
            }
        }else{
            if(qcsc==1){
                message(date(),"\tQC-RLSC...")
                res <- doQCRLSC(para = para)
                para <- res$metaXpara
                fig <- plotQCRLSC(para = para)
                message(date(),"\tnormalization after QC-RLSC...")
                para <- metaX::normalize(para = para, method = nor.method,
                                         valueID="valueNorm",...)
                
                if(nor.method=="combat"){
                    message(date(),"\tDo missing value imputation after combat normalization...")
                    para <- missingValueImpute(x = para,valueID="valueNorm",cpu=cpu)
                }
                message("The CV distribution after normalization: ")
                printCV(para,valueID = "valueNorm")
                prefix_tmp <- para@prefix
                para@prefix <- paste(para@prefix,"-norm",sep="")
                metaX::plotCV(para,valueID="valueNorm")
                para@prefix <- prefix_tmp
                para <- filterQCPeaksByCV(para,cvFilter = cvFilter,
                                          valueID = "valueNorm")
            }else if(qcsc==2){
                message(date(),"\tSVR-based batch correction...")
                res <- svrNormalize(para = para,ntop=5,cpu=cpu)
                para <- res$metaXpara
                fig <- plotQCRLSC(para = para)
                message(date(),"\tnormalization after SVR-based batch correction...")
                para <- metaX::normalize(para = para, method = nor.method,
                                         valueID="valueNorm",...)
                if(nor.method=="combat"){
                    message(date(),"\tDo missing value imputation after combat normalization...")
                    para <- missingValueImpute(x = para,valueID="valueNorm",cpu=cpu)
                }
                message("The CV distribution after normalization: ")
                printCV(para,valueID = "valueNorm")
                prefix_tmp <- para@prefix
                para@prefix <- paste(para@prefix,"-norm",sep="")
                metaX::plotCV(para,valueID="valueNorm")
                para@prefix <- prefix_tmp
                para <- filterQCPeaksByCV(para,cvFilter = cvFilter,
                                          valueID = "valueNorm")
            }else if(qcsc==3){
                
                message(date(),"\tComBat batch correction...")
                para <- batchCorrect(para,valueID = "value",use_class = FALSE,
                                     cpu=cpu,
                                     impute_method=para@missValueImputeMethod)
                para@peaksData$valueNorm <- para@peaksData$value
                
                message(date(),"\tnormalization after ComBat batch correction...")
                para <- metaX::normalize(para = para, method = nor.method,
                                         valueID="valueNorm",...)
                
                
                message("The CV distribution after normalization: ")
                printCV(para,valueID = "valueNorm")
                prefix_tmp <- para@prefix
                para@prefix <- paste(para@prefix,"-norm",sep="")
                metaX::plotCV(para,valueID="valueNorm")
                para@prefix <- prefix_tmp
                para <- filterQCPeaksByCV(para,cvFilter = cvFilter,
                                          valueID = "valueNorm")
                
            }else{
                stop("Not valid qcsc value!")
            }
            
        }
        
        if(qcsc!=3){
            s3_sub1 <- addTo(s3_sub1,newParagraph("The data was normalized by QC-based batch correction."))
            s3_sub1 <- addTo(s3_sub1,newTable(res$cvBatch,"CV summary after QC-based batch correction for each batch."))
            s3_sub1 <- addTo(s3_sub1,newTable(res$cvAll,"CV summary after QC-based batch correction for all samples."))
            fig1 <- paste("data/",basename(fig$fig),sep="")
            fig2 <- paste("data/",basename(fig$highfig),sep="")
            s3_sub1 <- addTo(s3_sub1,
                             newFigure(fig1,"QC-based batch correction figure.",
                                       fileHighRes = fig2))
        }
        
    }else{
        para@peaksData$valueNorm <- para@peaksData$value
        para <- metaX::normalize(para = para, method = nor.method,
                                 valueID="valueNorm",...)
        if(nor.method=="combat"){
            message(date(),"\tDo missing value imputation after combat normalization...")
            para <- missingValueImpute(x = para,valueID="valueNorm",cpu=cpu)
        }
        
        message("The CV distribution after normalization: ")
        printCV(para,valueID = "valueNorm")
        prefix_tmp <- para@prefix
        para@prefix <- paste(para@prefix,"-norm",sep="")
        metaX::plotCV(para,valueID="valueNorm")
        para@prefix <- prefix_tmp
        para <- filterQCPeaksByCV(para,cvFilter = cvFilter,
                                  valueID = "valueNorm")
    }
    
    
    if(pclean==TRUE){
        if(hasQC(para)){
            message("data clean...")
            para <- dataClean(para,valueID = "valueNorm",...)
        }
    }
    
    if(doQA==TRUE){
        if(hasQC(para)){
            s3_sub1 <- addTo(s3_sub1,
                             newParagraph("Data quality check after 
                                          normalization and pre-processing."))
            
            
            message("plot correlation heatmap...")
            prefix <- para@prefix
            para@prefix <- paste(prefix,"-nor",sep="")
            #save(para,file="para.rda")
            fig <- plotCorHeatmap(para = para,valueID = "valueNorm",anno = TRUE,
                                  samples = NA, # QC samples
                                  height = 6,width = 6,cluster = FALSE,
                                  classCol=classCol)
            fig1 <- paste("data/",basename(fig$fig),sep="")
            fig2 <- paste("data/",basename(fig$highfig),sep="")
            s3_sub1 <- addTo(s3_sub1,
                             newFigure(fig1,"Correlation heatmap.",
                                       fileHighRes = fig2))
            
            
            message("plot peak intensity distribution...")
            pp <- para
            pp@peaksData$value <- pp@peaksData$valueNorm
            fig <- plotIntDistr(pp,width = fig_width_boxplot)
            fig1 <- paste("data/",basename(fig$fig),sep="")
            fig2 <- paste("data/",basename(fig$highfig),sep="")
            s3_sub1 <- addTo(s3_sub1,
                             newFigure(fig1,"Peak intensity distribution.",
                                       fileHighRes = fig2))
            
            message("plot TIC distribution...")
            fig <- plotPeakSumDist(para,valueID = "valueNorm",height = 3.7,width = 5.5)
            fig1 <- paste("data/",basename(fig$fig),sep="")
            fig2 <- paste("data/",basename(fig$highfig),sep="")
            s3_sub1 <- addTo(s3_sub1,
                             newFigure(fig1,"TIC distribution.",
                                       fileHighRes = fig2))
            
            message("plot average intensity distribution...")
            fig <- plotPeakMeanDist(para,valueID = "valueNorm")
            fig1 <- paste("data/",basename(fig$fig),sep="")
            fig2 <- paste("data/",basename(fig$highfig),sep="")
            s3_sub1 <- addTo(s3_sub1,
                             newFigure(fig1,"Average intensity distribution.",
                                       fileHighRes = fig2))
            
            
            
            fig <- plotHeatMap(pp,valueID="valueNorm",log=TRUE,rmQC=FALSE,
                               #scale="row",
                               #clustering_distance_rows="euclidean",
                               #clustering_distance_cols="euclidean",
                               #clustering_method="ward.D2",
                               classCol=classCol,
                               show_colnames=FALSE)
            fig1 <- paste("data/",basename(fig$fig),sep="")
            fig2 <- paste("data/",basename(fig$highfig),sep="")
            s3_sub1 <- addTo(s3_sub1,
                             newFigure(fig1,"Heatmap.",
                                       fileHighRes = fig2))
            
            
            rm(pp)
            para@prefix <- prefix
        }else{
            message("plot correlation heatmap...")
            prefix <- para@prefix
            para@prefix <- paste(prefix,"-nor-batchcheck",sep="")
            #save(para,file="debug_p.rda")
            fig_tmp <- plotCorHeatmap(para = para,valueID = "valueNorm",anno = TRUE,
                                      samples = NULL,
                                      height = 6,width = 6,sortBy="batch",cluster = FALSE,
                                      classCol=classCol)
            para@prefix <- prefix 
        }
    }
    
    
    makeMetaboAnalystInput(para,rmQC=out.rmqc,valueID="valueNorm",
                           prefix="norm")
    #save(para,t,file="test.rda")
    save(para,t,scale,center,pcaLabel,classCol,file=paste(para@outdir,"/","pca.rda",sep=""))
    ppca <- transformation(para,method = t,valueID = "valueNorm")
    ppca <- metaX::preProcess(ppca,scale = scale,center = center,
                              valueID = "valueNorm")
    fig <- metaX::plotPCA(ppca,valueID = "valueNorm",scale = "none",batch = TRUE,
                          rmQC = FALSE,label = pcaLabel,classColor = classCol,...)
    prefix_bak <- ppca@prefix
    ppca@prefix <- paste(ppca@prefix,"-nobatch",sep="")
    fig_nobatch <- metaX::plotPCA(ppca,valueID = "valueNorm",scale = "none",batch = FALSE,
                          rmQC = FALSE,label = pcaLabel,classColor = classCol,...)
    
    
    
    ## plot heatmap
    plotHeatMap(ppca,valueID="valueNorm",log=FALSE,rmQC=FALSE,
                #scale="row",
                #clustering_distance_rows="euclidean",
                #clustering_distance_cols="euclidean",
                #clustering_method="ward.D2",
                classCol=classCol,
                show_colnames=FALSE)
    
    ppca@prefix <- prefix_bak
    
    
    fig1 <- paste("data/",basename(fig$fig),sep="")
    fig2 <- paste("data/",basename(fig$highfig),sep="")
    s3_sub1 <- addTo(s3_sub1,
                     newFigure(fig1,"PCA.",
                               fileHighRes = fig2))
    s3 <- addTo(s3,s3_sub1)
    
    ## remove QC
    ppca@prefix <- paste(ppca@prefix,"-noqc",sep="")
    fig <- metaX::plotPCA(ppca,valueID = "valueNorm",scale = "none",batch = TRUE,
                          rmQC = TRUE,label = pcaLabel,classColor = classCol,...)
    ## no QC, no batch
    ppca@prefix <- paste(ppca@prefix,"-nobatch",sep="")
    fig <- metaX::plotPCA(ppca,valueID = "valueNorm",scale = "none",batch = FALSE,
                          rmQC = TRUE,label = pcaLabel,classColor = classCol,...)
    #plotPLSDA(para)
    
    if(!is.null(para@ratioPairs) && para@ratioPairs != ""){
        ## only perform this analysis when a user provides comparison group
        ## information. 
        para <- peakStat(para = para,plsdaPara = plsdaPara,doROC=doROC,pcaLabel=pcaLabel,classColor = classCol)
        
        ## add identification result
        if(!is.null(idres)){
            para <- addIdentInfo(para=para,file=idres)
        }
        
        s3_sub2 <- newSubSection(asStrong("Metabolite quantification and identification."))
        s3_sub2 <- addTo(s3_sub2,
                         newParagraph("This part contains the quantification and identification."))
        
        qf <- paste("data/",para@prefix,"-quant.txt",sep="")
        show_quant <- para@quant[1:20,]
        names(show_quant) <- gsub(pattern = "_",replacement = "\n",x=names(show_quant))
        s3_sub2 <- addTo(s3_sub2,
                         newTable(show_quant,"Metabolite quantification result.",file=qf))
        
        s3 <- addTo(s3,s3_sub2)
    }
    
    
    if(saveRds){
        saveRDS(para,file = paste(para@outdir,"/",para@prefix,"-result.rds",
                                  sep=""))
    }
    
    report<-addTo(report,s1,s2,s3)
    message("Write report to file:",paste(raw_outdir,"/report",sep=""))
    writeReport(report, filename = paste(raw_outdir,"/report",sep=""))
    
    
    message("Print information about the current R session:")
    writeLines(capture.output(sessionInfo()), 
               paste(para@outdir,"/",para@prefix,"-sessionInfo.txt",sep=""))
    
    return(para)         
             
}
)


.tTest=function(x,y,log2=TRUE){
    if(log2){
        sv <- !is.na(x) & x>0 & is.finite(x)
        x <- log2(x[sv])
        y <- y[sv]
    }
    z <- t.test(x~y)
    return(z$p.value)
}

.tTestPair=function(x,y,log2=TRUE,pair){
    #save(x,y,pair,file="ttest.rda")
    if(log2){
        sv <- !is.na(x) & x>0 & is.finite(x)
        x <- log2(x[sv])
        y <- y[sv]
        pair <- pair[sv]
    }
    dat <- data.frame(x=x,y=y,pair=pair,stringsAsFactors = FALSE)
    dat2 <- spread(data = dat,value = x,key = y) %>% mutate(pair=NULL)
    
    z <- t.test(dat2[,1],dat2[,2],paired=TRUE)
    return(z$p.value)
}

.wilcoxTestPair=function(x,y,pair){

    dat <- data.frame(x=x,y=y,pair=pair,stringsAsFactors = FALSE)
    dat2 <- spread(data = dat,value = x,key = y) %>% mutate(pair=NULL)
    
    z <- wilcox.test(dat2[,1],dat2[,2],paired=TRUE)
    return(z$p.value)
}


isValid=function(x){
    return(x>0 & !is.na(x) & is.finite(x))
    
}

tuneSpline = function(x,y,span.vals=seq(0.1,1,by=0.05)){
    mae <- numeric(length(span.vals))
    
    crossEva <- function(span,x,y) {
        
        fun.fit <- function(x,y,span) {smooth.spline(x = x,y =y ,spar = span)}
        fun.predict <- function(fit,x0) {predict(fit,x0)$y}
        y.cv <- bootstrap::crossval(x,y,fun.fit,fun.predict,span=span,
                                    ngroup = length(x))$cv.fit
        fltr <- !is.na(y.cv)
        return(mean(abs(y[fltr]-y.cv[fltr])))
        
    }
    
    mae <- sapply(span.vals,crossEva,x=x,y=y)
    span <- span.vals[which.min(mae)]
    return(span)
}


.glogfit=function (x, a = 1, inverse = FALSE) {
    if (inverse) {
        out <- 0.25 * exp(-x) * (4 * exp(2 * x) - (a * a))
    }else{
        out <- log((x + sqrt(x^2 + a^2))/2)
    }
    return(out)
}

myLoessFit = function(x,y,newX,span.vals=seq(0.1,1,by=0.05),log=TRUE,a=1){
    if(log==TRUE){
        y <- .glogfit(y,a=a)
    }
    #sp.obj <- smooth.spline(x,y,spar = tuneSpline(x,y,span.vals = span.vals))
    sp.obj <- smooth.spline(x,y,cv = TRUE)
    
    valuePredict=predict(sp.obj,newX)
    if(log==TRUE){
        valuePredict$y <- .glogfit(valuePredict$y,a = a,inverse = TRUE)
    }
    return(valuePredict$y)
}



##' @title Remove samples from the metaXpara object
##' @description Remove samples from the metaXpara object
##' @rdname removeSample
##' @docType methods
##' @param para A metaXpara object.
##' @param rsamples The samples needed to be removed
##' @param ... Other argument
##' @return A metaXpara object.
##' @exportMethod
##' @author Bo Wen \email{wenbostar@@gmail.com}
##' @examples
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
##' sampleListFile(para) <- sfile
##' para <- reSetPeaksData(para)
##' new_para <- removeSample(para,rsamples=c("batch01_QC01"))
setGeneric("removeSample",function(para,rsamples,...) 
    standardGeneric("removeSample"))
##' @describeIn removeSample
setMethod("removeSample", signature(para = "metaXpara"),
          function(para,rsamples,...){
    if(is.null(para@peaksData)){
        para@rawPeaks <- para@rawPeaks[,! c(names(para@rawPeaks) %in% rsamples)]
    }else{
        para@peaksData <- dplyr::filter(para@peaksData,!sample %in% rsamples)
    }
    return(para)
})


##' @title Add identification result into metaXpara object
##' @description Add identification result into metaXpara object
##' @rdname addIdentInfo
##' @docType methods
##' @param para A metaXpara object.
##' @param file The file name which contains the identification result
##' @param ... Other argument
##' @return A metaXpara object.
##' @exportMethod
##' @author Bo Wen \email{wenbostar@@gmail.com}
setGeneric("addIdentInfo",function(para,file,...) 
    standardGeneric("addIdentInfo"))
##' @describeIn addIdentInfo
setMethod("addIdentInfo", signature(para = "metaXpara",file="character"),
          function(para,file,...){
    
    if(str_detect(tolower(file),pattern = ".csv$")){
        message("The identification result file is csv format!")
        res <- read.csv(file,stringsAsFactors=FALSE)
    }else{
        message("The identification result file is txt format!")
        res <- read.delim(file,stringsAsFactors=FALSE)
    }
    ## whether there is a column named "ID"
    checkID <- which(names(res)=="ID")
    if(length(checkID)==0){
        ## no found, ## the first column must be peak ID
        names(res)[1] <- "ID"
    }else if(length(checkID)>=2){
        stop("There are more than 2 column named ID,please check your data!\n")
    }else{
        ## don't need to do anything.
    }
    
    ## please note that maybe there are more than one rows corresponding to 
    ## one peak
    
    ## only output merged result to a file
    if(!is.null(para@quant)){
        ## merge
        ## firstly, check whether there is duplicate column in the two table
        cID <- sum(names(res) %in% names(para@quant))
        if(cID==1){
            idres <- dplyr::left_join(para@quant,res,by="ID")
            message("Peak number: ",length(unique(para@quant$ID)))
            message("Identified peak: ",sum(unique(para@quant$ID) %in% res$ID))
            rfile <- paste(para@outdir,"/",para@prefix,
                           "-quant-identification.txt",sep="")
            message("Merge identification and quantification result and save to", 
                " file ",rfile," ...")
            write.table(idres,file=rfile,col.names = TRUE,row.names = FALSE,
                        quote=FALSE,sep="\t")
            
        }else{
            stop("There are more than two column contained in quantification", 
                 "result! Please check your data!\n")
        }
        
    }
    
    para@idres <- res
    return(para)
})


##' @title Plot boxplot for each feature
##' @description Plot boxplot for each feature
##' @rdname plotPeakBox
##' @docType methods
##' @param para A metaXpara object
##' @param samples Sample class used
##' @param log Whether log transform or not
##' @param ... Additional parameters
##' @return The output figure name.
##' @exportMethod
##' @author Bo Wen \email{wenbostar@@gmail.com}
##' @examples
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)[1:20,]
##' sampleListFile(para) <- sfile
##' para <- reSetPeaksData(para)
##' para <- missingValueImpute(para)
##' addValueNorm(para) <- para
##' plotPeakBox(para,samples=c("S","C"))
setGeneric("plotPeakBox",function(para,samples,log=FALSE,...) 
    standardGeneric("plotPeakBox"))
##' @describeIn plotPeakBox
setMethod("plotPeakBox", signature(para = "metaXpara"),
          function(para,samples,log=FALSE,...){
    
    para@peaksData$class <- as.character(para@peaksData$class)
    peaksData <- dplyr::filter(para@peaksData,class %in% samples)
    
    ## whther the data is valid
    if(log==TRUE){
        ## can't have zero or negative value
        if(sum(peaksData$valueNorm<=0)>=1){
            message("value <= 0 ",sum(peaksData$valueNorm<=0))
            stop("cannot handle with <=0 value when log...\n")
        }else{
            peaksData$valueNorm <- log2(peaksData$valueNorm)
        }
    }
    
    fig <- paste(para@outdir,"/",para@prefix,"-peakboxplot-",
                 paste(samples,collapse = "_"),".pdf",sep="")
    pdf(fig,width = 4,height = 4)
    
    figres <- ddply(peaksData,.(ID),function(dat){
        dat$class <- as.factor(dat$class)
        ggobj <- ggplot(data=dat,aes(x=class,y=valueNorm))+
                    geom_boxplot(width = 0.3)+
                    #theme(legend.justification=c(1,1), 
                    #      legend.position=c(1,1))+
                    geom_jitter(position = position_jitter(width = 0.1),
                                aes(colour=class))+
                    ggtitle(dat$ID[1])+
                    ylab("Intensity")
                                
        print(ggobj)
        NA
    })
    dev.off()
    #return(para)
})



##' @title Plot the total peak intensity distribution
##' @description Plot the total peak intensity distribution
##' @rdname plotPeakSumDist
##' @docType methods
##' @param para A metaXpara object.
##' @param valueID The name of the column used
##' @param width Width of the figure
##' @param height Height of the figure
##' @param ... Other argument
##' @return The output figure name.
##' @exportMethod
##' @author Bo Wen \email{wenbostar@@gmail.com}
##' @examples
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
##' sampleListFile(para) <- sfile
##' para <- reSetPeaksData(para)
##' plotPeakSumDist(para)
setGeneric("plotPeakSumDist",function(para,valueID="value",width=6,height=4,...) 
    standardGeneric("plotPeakSumDist"))
##' @describeIn plotPeakSumDist
setMethod("plotPeakSumDist", signature(para = "metaXpara"),
          function(para,valueID="value",width=8,height=5,...){
    
    dat <- ddply(para@peaksData,.(sample,batch,class,order),here(summarise),
                   sum=sum(get(valueID)))
    
    dat$batch <- as.character(dat$batch)
    dat$class <- as.character(dat$class)
    dat$class[is.na(dat$class)] <- "QC"
    
    outliers.coef <- 1.5
    qs <- quantile(dat$sum,c(.25,.75),na.rm=TRUE)
    iqr <- qs[2] - qs[1]
    dat$outlier <- dat$sum < qs[1]-outliers.coef*iqr | 
        dat$sum > qs[2]+outliers.coef*iqr
    
    fig <- paste(para@outdir,"/",para@prefix,"-peakSumDist.png",sep="")
    highfig <- sub(pattern = "png$",replacement = "pdf",x = fig)
    pdf(file = highfig,width = width,height = height)
    
    write.table(dat,file=paste(para@outdir,"/",para@prefix,"-peakSumDist.txt",sep=""),
                quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
    
    ggobj <- ggplot(data=dat,aes(x=order,y=sum,colour=class,shape=batch))+
        geom_point()+
        theme_bw()+
        scale_shape_manual(values=1:n_distinct(dat$batch))+
        geom_text(aes(label=ifelse(outlier,order,"")),hjust=-0.2,size=4)+
        ylab("Total peak intensity")
    
    print(ggobj)
    dev.off()
    png(filename = fig,width = width,height = height,units = "in",res = 150)
    print(ggobj)
    dev.off()
    res <- list(fig=fig,highfig=highfig)
    return(res)
    
})



##' @title Plot the average peak intensity distribution
##' @description Plot the average peak intensity distribution
##' @rdname plotPeakMeanDist
##' @docType methods
##' @param para A metaXpara object.
##' @param valueID The name of the column used
##' @param width Width of the figure
##' @param height Height of the figure
##' @param ... Other argument
##' @return The output figure name.
##' @exportMethod
##' @author Bo Wen \email{wenbostar@@gmail.com}
##' @examples
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
##' sampleListFile(para) <- sfile
##' para <- reSetPeaksData(para)
##' plotPeakMeanDist(para)
setGeneric("plotPeakMeanDist",function(para,valueID="value",width=6,height=4,...) 
    standardGeneric("plotPeakMeanDist"))
##' @describeIn plotPeakSumDist
setMethod("plotPeakMeanDist", signature(para = "metaXpara"),function(para,valueID="value",width=8,height=5,...){
    dat <- ddply(para@peaksData,.(sample,batch,class,order),here(summarise),
                 mean=mean(get(valueID),na.rm = TRUE))
    
    dat$batch <- as.character(dat$batch)
    dat$class <- as.character(dat$class)
    dat$class[is.na(dat$class)] <- "QC"
    
    outliers.coef <- 1.5
    qs <- quantile(dat$mean,c(.25,.75),na.rm=TRUE)
    iqr <- qs[2] - qs[1]
    dat$outlier <- dat$mean < qs[1]-outliers.coef*iqr | 
        dat$mean > qs[2]+outliers.coef*iqr
    
    fig <- paste(para@outdir,"/",para@prefix,"-peakMeanDist.png",sep="")
    highfig <- sub(pattern = "png$",replacement = "pdf",x = fig)
    pdf(file = highfig,width = width,height = height)
    
    write.table(dat,file=paste(para@outdir,"/",para@prefix,"-peakMeanDist.txt",sep=""),
                quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
    
    ggobj <- ggplot(data=dat,aes(x=order,y=mean,colour=class,shape=batch))+
        geom_point()+
        theme_bw()+
        scale_shape_manual(values=1:n_distinct(dat$batch))+
        geom_text(aes(label=ifelse(outlier,order,"")),hjust=-0.2,size=4)+
        ylab("Average peak intensity")
    
    print(ggobj)
    dev.off()
    png(filename = fig,width = width,height = height,units = "in",res = 150)
    print(ggobj)
    dev.off()
    res <- list(fig=fig,highfig=highfig)
    return(res)

})

##' @title dataClean
##' @description dataClean
##' @rdname dataClean
##' @docType methods
##' @param para A metaXpara object.
##' @param valueID The name of the column used
##' @param sd.factor The factor used to filter peak based on SD
##' @param snr The threshold to filter peak
##' @param ... Other argument
##' @return A metaXpara object.
##' @exportMethod
##' @author Bo Wen \email{wenbostar@@gmail.com}
##' @examples
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
##' sampleListFile(para) <- sfile
##' para <- reSetPeaksData(para)
##' p <- dataClean(para)
setGeneric("dataClean",function(para,valueID="value",sd.factor=3,snr=1,...) 
    standardGeneric("dataClean"))
##' @describeIn dataClean
setMethod("dataClean", signature(para = "metaXpara"),
          function(para,valueID="value",sd.factor=3,snr=1,...){
    
    peaksData <- para@peaksData
    ## true sample
    peaksData1 <- dplyr::filter(peaksData,!is.na(class))
    ## QC sample
    peaksData2 <- dplyr::filter(peaksData,is.na(class)) 
    
    dat1 <- ddply(peaksData1,.(ID),here(summarise),
                  mean=mean(get(valueID),na.rm = TRUE),
                  sd=sd(get(valueID),na.rm = TRUE))
    
    dat2 <- ddply(peaksData2,.(ID),here(summarise),
                  meanQC=mean(get(valueID),na.rm = TRUE),
                  sdQC=sd(get(valueID),na.rm = TRUE))
    
    dat <- dplyr::full_join(dat1,dat2,by="ID")
    dat_remove <- dplyr::mutate(dat,remove=abs(meanQC-mean) >sd.factor*sd,
                                SNR=sd/sdQC)
    
    tmpPeaksdData <- peaksData
    tmpPeaksdData$val <- peaksData[,valueID]
    kw_pvalue <- ddply(tmpPeaksdData,.(ID),summarise,
                       pvalue=kruskal.test(val~batch)$p.value)
    dat <- dplyr::full_join(dat_remove,kw_pvalue,by="ID")
    
    rID <- dat_remove$ID[dat_remove$remove==TRUE | dat_remove$SNR < snr]
    
    save2file <- paste(para@outdir,"/",para@prefix,"-dataClean.txt",sep="")
    message("Save removed data to file:",save2file)
    write.table(dat_remove[dat_remove$ID%in%rID,],file = save2file,
                quote = FALSE,sep="\t",
                row.names = FALSE)
    
    message("remove ID: ",length(rID))
    
    ## output the removed feature
    rPeak <- dplyr::filter(peaksData,ID%in%rID)
    rPeak$class <- as.character(rPeak$class)
    rPeak$class[is.na(rPeak$class)] <- "QC"
    rPeak$batch <- as.character(rPeak$batch)
    
    fig <- paste(para@outdir,"/",para@prefix,"-dataClean.pdf",sep="")
    pdf(fig,width = 6,height = 6)
    for(id in unique(rPeak$ID)){
        pdat <- dplyr::filter(rPeak,ID==id)
        pdat$val <- pdat[,valueID]
        tmp <- dplyr::filter(dat_remove,ID==id)
        gtitle <- paste(id,tmp$remove,tmp$SNR,sep="_")
        ggobj <- ggplot(pdat,aes(x=order,y=val,
                                 colour=class,shape=batch))+
            geom_point()+
            scale_shape_manual(values=1:n_distinct(pdat$batch))+
            ylab(valueID)+
            ggtitle(label = gtitle)
        print(ggobj)    
        
    }
    dev.off()
    
    para@peaksData <-dplyr::filter(peaksData,!ID%in%rID)
    return(para)
})



##' @title checkQCPlot
##' @description Plot figure for quantification and identification result
##' @param f1 a file contained the quantification result of metaX
##' @param f2 a file contained the metabolite identification result
##' @param fig the file name of output figure
##' @param group the group name
##' @return none
##' @export
checkQCPlot=function(f1,f2=NULL,fig="test.png",group=NULL){
    
    
    a <- read.delim(f1,stringsAsFactors = FALSE)
    if(is.null(group)){
        print(unique(a$sample))
        stop("Please set the group value!")
    }
    a <- a %>% filter(sample==group)
    
    a$rt <- sapply(a$ID,function(x){as.numeric(strsplit(x,split = "_")[[1]][1])})
    a$mz <- sapply(a$ID,function(x){
        x=gsub(pattern = "m/z",replacement = "",x=x)
        x=gsub(pattern = "n",replacement = "",x=x)
        as.numeric(strsplit(x,split = "_")[[1]][2])})
    
    nbar = 60
    xmin = min(a$rt)
    xmax = max(a$rt)
    xmean = (xmin+xmax)/2.0
    
    png(filename = fig,width = 1000,height = 800,res = 120)
    if(!is.null(f2)){
        par(mfrow=c(5,1),mgp=c(1.6,0.6,0))
    }else{
        par(mfrow=c(4,1),mgp=c(1.6,0.6,0))
    }
    
    
    par(mar=c(0,4,0.2,0.5))
    plot(a$rt,a$VIP,type="h",col="gray",xlim=c(xmin,xmax),xaxt="n",xlab="",ylab="VIP")
    text(xmean,0.7*(par("usr")[4]-par("usr")[3]),labels = "VIP vs RT")
    
    d1 <- a %>% filter(abs(log2(ratio))>=log2(1.2),t.test_p.value <= 0.05) 
    
    if(nrow(d1)>=1){
        points(d1$rt,d1$VIP,col="red",cex=0.6)
    }
    
    if(!is.null(f2)){
        idres <- read.delim(f2,stringsAsFactors = FALSE)
        idres <- idres %>% filter(!is.na(Compound.ID))
        d2 <- d1 %>% filter(ID %in% idres$ID)
        if(nrow(d2)>=1){
            points(d2$rt,d2$VIP,col="blue",cex=0.6)
        }
        
        ## identification
        par(mar=c(0,4,0,0.5))
        d3 <- a %>% filter(ID %in% idres$ID)
        plot(a$rt,a$mz,type="p",col="gray",xlim=c(xmin,xmax),xaxt="n",xlab="",cex=0.6,ylab="m/z")
        text(xmean,0.7*(par("usr")[4]-par("usr")[3]),labels = "MZ vs RT")
        points(d3$rt,d3$mz,cex=0.6,col="blue")
    }
    
    plot(a$rt,a$mz,type="p",col="gray",xlim=c(xmin,xmax),xaxt="n",xlab="",cex=0.6,ylab="m/z")
    text(xmean,0.7*(par("usr")[4]-par("usr")[3]),labels = "MZ vs RT(red=DEB)")
    if(nrow(d1)>=1){
        points(d1$rt,d1$mz,col="red",cex=0.6)
    }
    
    
    ### barplot
    if(nrow(d1)>=1){
        par(mar=c(0,4,0,0.5))
        d1$rt %>% hist(nclass=nbar,xlim=c(xmin,xmax),main="",xaxt="n",xlab="")
        text(xmean,0.7*(par("usr")[4]-par("usr")[3]),labels = "DEB distribution")
        box()
    }
    
    
    par(mar=c(3,4,0,0.5))
    a$rt %>% hist(nclass=nbar,xlim=c(xmin,xmax),main="",xlab="rt")
    text(xmean,0.7*(par("usr")[4]-par("usr")[3]),labels = "All features distribution")
    box()
    dev.off()
    
}


##' @title checkPvaluePlot
##' @description Plot pvalue distribution
##' @param f1 a file contained the quantification result of metaX
##' @param fig the file name of output figure
##' @param group the group name
##' @return none
##' @export
checkPvaluePlot=function(file=NULL,group=NULL,fig="pvalue.png"){
    a <- read.delim(file,stringsAsFactors = FALSE)
    a <- a %>% filter(sample==group)
    png(filename = fig,width = 600,height = 600,res = 120)
    par(mfrow=c(4,1),mgp=c(1.6,0.6,0))
    
    a <- a %>% select(ID,t.test_p.value,wilcox.test_p.value,t.test_p.value_BHcorrect,wilcox.test_p.value_BHcorrect)
    a <- a %>% gather(class,pvalue,-ID) %>% 
        mutate(class=gsub(pattern="_p.value",replacement="",x=class)) %>%
        mutate(class=gsub(pattern="correct",replacement="",x=class))
    gg <- ggplot(data=a,aes(x=pvalue))+
            geom_histogram(colour="blue",fill="white")+
            facet_grid(class~.,scales="free_y")+
            theme_bw()
    print(gg)
    dev.off()
    
}


checkSampleList=function(file=NULL){
    if(!is.data.frame(file)){
        a <- read.delim(file,stringsAsFactors = FALSE)
    }else{
        a <- file
    }
    ## check sample name
    sort_sample <- sort(a$sample)
    ns1 = length(sort_sample)
    ns2 = length(unique(sort_sample))
    stop_flag = FALSE
    
    message("Check sample ...")
    if(ns1!=ns2){
        message("\tFind duplicate sample(s) in your sample list:")
        print(sort_sample[duplicated(sort_sample)])
        stop_flag = TRUE
    }
    
    ## check order duplicate
    sort_order <- sort(a$order)
    no1 = length(sort_order)
    no2 = length(unique(sort_order))
    
    message("Check order ...")
    if(no1!=no2){
        message("\tFind duplicate order(s) in your sample list:")
        print(sort_order[duplicated(sort_order)])
        stop_flag = TRUE
    }
    
   
    
    ## check file name
    message("Check file name of sample list ...")
    file_name = c("sample","batch","class","order")
    ni = sum(file_name %in% names(a))
    if(ni!=4){
        message("\tPlease provide valid file name:")
        print(file_name)
        stop_flag = TRUE
    }
    
    if(stop_flag){
        stop("Please check your sample list. It's not valid format!")
        
    }
    
}

check_input_data_format=function(x){
    
    cat("Check sample names and feature names!\n")
    if(is.data.frame(x)){
        a <- x
    }else{
        a <- read.delim(x,stringsAsFactors = FALSE,check.names = FALSE)
    }
    ## sample name
    nsample <- sum(table(names(a)) >=2)
    if(nsample >= 1){
        stop("There are multiple samples having the same name!!!\n")
    }
    
    nfeature <- sum(table(a$name) >=2)
    if(nfeature >= 1){
        stop("There are multiple features having the same name!!!\n")
    }
}


##' @title Plot chromatogram of total ion count
##' @description Plot chromatogram of total ion count
##' @param file
##' @param corBatch Code colour according to batch number
##' @param corOrder Code colour according to injection order
##' @return none
##' @export
##' @author Bo Wen \email{wenbostar@@gmail.com}
plotTIC=function(file,corBatch=FALSE,corOrder=TRUE){
    
    ## 
    a <- read.delim(file,stringsAsFactors = FALSE)
    
    msdata <- a %>% group_by(sample,batch,order) %>% 
                    do(readMSfile(msfile=.$sample[1]))
    
    msdata$batch <- as.character(msdata$batch)
    msdata$order <- as.character(msdata$order)
    
    if(corBatch){
        ggplot(data = msdata,aes(x = retentionTime,y = totIonCurrent, 
                                 colour = batch)) +
            geom_line() + 
            theme(legend.position="bottom") +
            guides(col = guide_legend(ncol = 10))
    }else if(corOrder){
        ggplot(data = msdata,aes(x = retentionTime,y = totIonCurrent, 
                                 colour = order)) +
            geom_line() + 
            theme(legend.position="bottom")+
            guides(col = guide_legend(ncol = 10))
    }
    
}


readMSfile=function(msfile){
    msdata <- header(openMSfile(msfile)) %>% filter(msLevel==1)
    return(msdata)
}


printCV=function(para,valueID="valueNorm"){
    
    peaksData <- para@peaksData
    if(any(is.na(peaksData$class))){
        peaksData$class <- as.character(peaksData$class)
        peaksData$class[is.na(peaksData$class)]<-"QC"
    }
    
    message("CV distribution: ")
    if(valueID=="value"){
        cvstat <- peaksData %>% 
            group_by(batch,class,ID) %>% 
            summarize(cv=sd(value,na.rm=TRUE)/mean(value,na.rm=TRUE)) %>% 
            ungroup() %>% 
            group_by(batch,class) %>% 
            summarize(median_cv=median(cv,na.rm=TRUE),mean_cv=mean(cv,na.rm=TRUE),n=n()) %>%
            as.data.frame()
        print(cvstat)
        cvstat2 <- peaksData %>% 
            group_by(class,ID) %>% 
            summarize(cv=sd(value,na.rm=TRUE)/mean(value,na.rm=TRUE)) %>% 
            ungroup() %>% 
            group_by(class) %>% 
            summarize(median_cv=median(cv,na.rm=TRUE),mean_cv=mean(cv,na.rm=TRUE),n=n()) %>%
            as.data.frame()
        print(cvstat2)
    }else{
        cvstat <- peaksData %>% 
            group_by(batch,class,ID) %>% 
            summarize(cv=sd(valueNorm,na.rm=TRUE)/mean(valueNorm,na.rm=TRUE)) %>% 
            ungroup() %>% 
            group_by(batch,class) %>% 
            summarize(median_cv=median(cv,na.rm=TRUE),mean_cv=mean(cv,na.rm=TRUE),n=n()) %>%
            as.data.frame()
        print(cvstat)
        cvstat2 <- peaksData %>% 
            group_by(class,ID) %>% 
            summarize(cv=sd(valueNorm,na.rm=TRUE)/mean(valueNorm,na.rm=TRUE)) %>% 
            ungroup() %>% 
            group_by(class) %>% 
            summarize(median_cv=median(cv,na.rm=TRUE),mean_cv=mean(cv,na.rm=TRUE),n=n()) %>%
            as.data.frame()
        print(cvstat2)
    }
}


##' @title Get a table from a metaXpara object.
##' @description Get a table from a metaXpara object.
##' @param para A metaXpara object
##' @param valueID The ID of value to extract
##' @return A data.frame
##' @export
##' @author Bo Wen \email{wenbostar@@gmail.com}
getTable=function(para,valueID="valueNorm",outfile=NULL, zero2NA = FALSE){
    if(zero2NA == TRUE){
        cat("Set zero as NA!\n")
        para@peaksData[,valueID][para@peaksData[,valueID] <= 0] <- NA
    }
    
    res <- para@peaksData %>% select_("ID","sample",valueID) %>% 
        spread_(key = "sample",value = valueID)
    n <- sum(is.na(para@peaksData[,valueID])) + sum(para@peaksData[,valueID]<=0,na.rm = TRUE)
    cat("value with NA or <= 0:",n,"\n")
    
    if(!is.null(outfile)){
        cat("Export data to file:",outfile,"\n")
        write_tsv(res,path = outfile)
    }
    return(res)
}


## put histograms on the diagonal
.panel.hist <- function(x, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    x <- x[!is.infinite(x) & !is.na(x)]
    h <- hist(x, nclass=25,plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
.panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    
    #& x>0 & y>0
    useValue <- !is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y) 
    x <- x[useValue]
    y <- y[useValue]
    r <- abs(cor(x, y,method = "sp"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

.panel.point <- function(x, y){
    #x >0 & y >0 #& is.finite(x) & is.finite(y)
    useValue <- !is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y) 
    x <- x[useValue]
    y <- y[useValue]
    points(x,y,cex=0.5,pch=16,col=rgb(255/255,69/255,0,0.4))
    abline(a = 0,b = 1)
    lines(smooth.spline(x,y),col="blue")
}


plotPairs=function(dat,log2=TRUE,addOne=TRUE,...){
    if(addOne){
        dat = dat + 1
    }
    if(log2){
        dat <- log2(dat)    
    }
    pairs(dat,upper.panel = .panel.cor,lower.panel = .panel.point,
          gap=0,
          #xlim=c(-3,3),ylim=c(-3,3),
          #labels=c("","IsobariQ","ProteinPilot"),
          cex.labels=1.3,font.labels=2,...)    

    
}











