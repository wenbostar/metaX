

##' @title Power Analysis 
##' @description Power Analysis
##' @param para An metaXpara object
##' @param group A vector of sample names
##' @param valueID The column name used
##' @param log A logical indicating whether transform the data with log2
##' @param maxInd max sample number
##' @param fdr The FDR threshold
##' @param showPlot Whether or not to print the figure to screen
##' @return An value
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @export
##' @examples
##' \dontrun{ 
##' library(reshape2)
##' library(dplyr)
##' a <- read.csv("http://www.metaboanalyst.ca/MetaboAnalyst/resources/data/power_example.csv")
##' peaksData <- melt(a,id.vars = c("Diet","Sample"),
##'     value.name = "value",variable.name = "ID")
##' peaksData <- dplyr::rename(peaksData,class=Diet,sample=Sample)
##' para <- new("metaXpara")
##' peaksData(para) <- peaksData
##' para <- missingValueImpute(para)
##' para <- metaX::normalize(para)
##' para <- transformation(para,valueID = "value")
##' para <- preProcess(para,scale = "pareto",valueID="value")
##' powerAnalyst(para,group=c("case","control"),log=FALSE,maxInd=200)
##' }
powerAnalyst=function(para,group,valueID="value",log=TRUE,maxInd=1000,fdr=0.1,showPlot=FALSE){
    
    if(length(group)!=2){
        stop("Only valid for two class group!\n")
    }
    peaksData <- para@peaksData
    peaksData <- dplyr::filter(peaksData,class %in% group)
    if(log==TRUE){
        peaksData[,valueID] <- log2(peaksData[,valueID])
    }
    
    myInlineTtest=function(x,y){
        tmp <- try(t.test(x~y, var.equal = TRUE))
        if(class(tmp) == "try-error") {
            message("x:")
            print(x)
            message("y:")
            print(y)
            return(NA)
        }else{
            return(tmp$statistic)
        }
    }
    
    statres <- ddply(peaksData,.(ID),here(summarise),
                     statistic=myInlineTtest(get(valueID),class))
    statres <- dplyr::filter(statres,!is.na(statistic))
    
    stats <- statres$statistic
    
    nsample <- dplyr::select(peaksData,sample,class) %>% dplyr::distinct()
    nsample <- table(nsample$class)
    n1 <- nsample[1]
    n2 <- nsample[2]
    print(nsample)
    pdD <- pilotData(statistics = stats,
                     samplesize = sqrt(n1+n2),
                     distribution="t",
                     df=n1+n2-2)
    
    
    ##
    Jpred <- c(3, 6, 10, 16, 24, 40, 60, 100, 150, seq(200, 1000, 100))
    inx <- which(Jpred == min(Jpred[Jpred>=maxInd]))
    Jpred <- Jpred[1:inx]
    
    res <- round(length(pdD@statistics)/2)
    ssD <- sampleSize(pdD, method="congrad", 
                      control=list(from=-6, to=6, resolution=res))
    N <- sqrt(Jpred/2)
    
    pi0 <- ssD@pi0
    if(fdr >= pi0){
        fdr <- signif(pi0-pi0/10, 3)
    }
    
    pwrD <- predictpower(ssD, samplesizes=N, alpha=fdr)
    
    expPWR <- predictpower(ssD, samplesizes=sqrt((n1+n2)/2/2), alpha=fdr)
    
    fig <- paste(para@outdir,"/",para@prefix,"-power.pdf",sep="")
    pdf(file = fig,width = 7,height  = 5)
    plot(pdD)
    ## density of effect sizes
    plot(ssD)
    
    ## 
    pdata <- data.frame(x=Jpred,y=pwrD)
    ggobj <- ggplot(pdata,aes(x=x,y=y))+geom_point()+geom_line()+
        xlab("Sample Size (per group)")+
        ylab("Predicted power")+
        geom_hline(yintercept=0.8,colour="gray")+
        geom_vline(xintercept=((n1+n2)/2),colour="red")+
        theme_bw()
    
    print(ggobj)
    
    dev.off()
    message(fig)
    if(showPlot){
        print(ggobj)
    }
    return(expPWR)
}