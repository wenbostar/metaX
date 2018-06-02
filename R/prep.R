

##' @title Pre-Processing
##' @description Pre-Processing
##' @param para An metaX object
##' @param t The method for transformation, 0=none, 1=log, 2=Cube root, 3=glog
##' @param scale The method of scaling: "auto", "range", "pareto", "vast", "level", "power","none"
##' @param center Centering
##' @param valueID The name of column used for transformation
##' @param ... Additional parameter
##' @return An new metaX object
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @export
##' @examples
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
##' sampleListFile(para) <- sfile
##' para <- reSetPeaksData(para)
##' para <- missingValueImpute(para)
##' para <- preProcess(para,valueID = "value",scale="uv")
preProcess=function(para,t=0,scale=c("auto", "range", "pareto", "vast", "level", "power","none"),
                    center=TRUE,valueID="valueNorm"){
    p <- transformation(para=para,method = t,valueID = valueID)
    p <- dataScaling(para=p,method = scale,center = center,valueID = valueID)
    return(p)
}

#Generalized log transformation function
logTransform=function(x){
    min.val <- min(abs(x[x!=0 & !is.na(x)]))/10;
    log2((x + sqrt(x^2 + min.val^2))/2)
}


# glog=function (x, a = 1, inverse = FALSE) {
#     if (inverse) {
#         out <- 0.25 * exp(-x) * (4 * exp(2 * x) - (a * a))
#     }else{
#         out <- log((x + sqrt(x^2 + a^2))/2)
#     }
#     return(out)
# }

## a data matrix ([data.frame object] row: molecules, col: samples or replicates)
dataScaling=function (para, method = c("auto", "range", "pareto", "vast", "level", "power","none"),
                      center=TRUE,valueID="value"){
    peaksData <- para@peaksData
    xyData <- peaksData %>% select_("ID","class","sample",valueID) %>% spread_("ID",valueID) 
    if(center==TRUE){
        rData <- xyData %>% select(-class,-sample) %>% scale(center=TRUE,scale = FALSE) %>% t()    
    }else{
        rData <- xyData %>% select(-class,-sample) %>% t()
    }
    
    if(method == "auto"){
        ## for each row - metabolite
        res <- apply(rData, 1, function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))
    }else if(method == "range"){
        res <- apply(rData, 1, function(x) (x - mean(x, na.rm = TRUE))/(range(x)[2] -  range(x)[1]))
    }else if(method == "pareto"){
        res <- apply(rData, 1, function(x) (x - mean(x, na.rm = TRUE))/sqrt(sd(x,na.rm = TRUE)))
    }else if(method == "vast"){
        res <- apply(rData, 1, function(x) mean(x, na.rm = TRUE) * (x - mean(x, na.rm = TRUE))/(sd(x, na.rm = TRUE)^2))
    }else if(method == "level"){
        res <- apply(rData, 1, function(x) (x - mean(x, na.rm = TRUE))/mean(x,na.rm = TRUE))
    }else if(method == "power"){
        res <- apply(rData, 1, function(x) sqrt(x) - mean(sqrt(x)))
    }else{
        res <- apply(rData, 1, function(x) {x})
    }
    
    scaleData <- cbind(xyData[,c("class","sample")],res)
    # res <- melt(data = scaleData,value.name = valueID,
    #            id.vars = names(xyData)[1:2],variable.name  = "ID")
    # res <- scaleData %>% gather_("ID",valueID,-class,-sample)
    res <- scaleData %>% gather(key="ID",value="intensity",-class,-sample) %>% 
                mutate_(.dots = setNames(list(~intensity),valueID)) %>% 
                mutate(intensity=NULL)
    peaksData[,valueID] <- NULL
    mres <- inner_join(res,peaksData,by=c("sample","class","ID"))
    
    if(nrow(mres)!=nrow(peaksData)){
        stop("scale error!\n")
    }
    para@peaksData <- mres
    return(para)
}


##' @title Data transformation
##' @description Data transformation
##' @param para An metaX object
##' @param method The method for transformation, 0=none, 1=log, 2=Cube root, 3=glog
##' @param valueID The name of column used for transformation
##' @param ... Additional parameter
##' @return An new metaX object
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @export
##' @examples
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
##' sampleListFile(para) <- sfile
##' para <- reSetPeaksData(para)
##' para <- missingValueImpute(para)
##' para <- transformation(para,valueID = "value")
transformation=function(para,method=1,valueID="valueNorm",...){
    if(method==1){
        ## Log transformation (generalized logarithm transformation or glog)
        para@peaksData[,valueID] <- logTransform(para@peaksData[,valueID])
    }else if(method==2){
        x <- para@peaksData[,valueID]
        min.val <- min(abs(x[x!=0 & !is.na(x)]))/10;
        x[is.na(x) | x <=0] <- min.val
        norm.data <- abs(x)^(1/3)
        para@peaksData[,valueID] <- norm.data
    }else if(method==0){
        ## don't transformate
    }else if(method==3){
        para <- doglog(para,valueID = valueID)
    }else{
        message("Please provide valid method!")
    }
    return(para)
}

.glog <- function(y,alpha,lambda){
    z <- log((y-alpha)+sqrt((y-alpha)^2+lambda))
}

.jglog <- function(y,y0,lambda){
    z <- .glog(y,y0,lambda)
    D <- log(sqrt((y-y0)^2+lambda))
    
    gmn <- exp(apply(D,2,mean,na.rm=T))
    zj <- z*gmn
    return(zj)
}

.SSE <- function(lambda,alpha,y){
    N <- dim(y)[2]
    len <- dim(y)[1]
    
    z <- .jglog(y,alpha,lambda)
    s <- 0
    mean_spec <- apply(z,1,mean,na.rm=T)
    
    s <- sum((z-mean_spec)^2,na.rm=T)
    
    #cat(lambda,"\t",s,"\n")
    
    return(s)
}

doglog=function(para,valueID="value",useQC=FALSE){
    
    if(useQC==TRUE){
        qcPeaksData <- para@peaksData %>% filter(is.na(class)) %>% 
            select_("ID",valueID,"sample") %>% 
            spread_("ID",valueID) %>%
            select(-sample)
    }else{
        qcPeaksData <- para@peaksData %>%  
            select_("ID",valueID,"sample") %>% 
            spread_("ID",valueID) %>%
            select(-sample)
    }
    samplePeaksData <- para@peaksData %>% 
        select_("ID",valueID,"sample") %>% 
        spread_("ID",valueID)
    row.names(samplePeaksData) <- samplePeaksData$sample
    
    samplePeaksData <- samplePeaksData %>% select(-sample)
    
    x <- t(qcPeaksData)
    
    y0 <- 0
    N <- ncol(x)
    L <- max(dim(x))
    
    scal_fact <- 1
    pow_fact <- 1
    offset <- min(x,na.rm=T)
    x <- x-offset
    
    step_threshold <- 1e-16
    
    small <- min(x,na.rm=T) #which is 0
    
    if(min(apply(t(x),2,var,na.rm=T))==0){
        varbs <- apply(t(x),2,var,na.rm=T)
        newminVar <- sort(unique(varbs))[2]
    }else{
        newminVar <- min(apply(t(x),2,var,na.rm=T))
    }
    
    low_lim <- -small^2
    upper_lim <- max(pmax(apply(t(x),2,var,na.rm=T),max(apply(t(x),2,var,na.rm=T))/newminVar))
    
    lambda <- optimize(.SSE,interval=c(low_lim,upper_lim),y0,x,tol=step_threshold)
    
    lambda <- as.numeric(lambda[[1]])
    
    lambda_std <- 5.0278*10^(-09)
    
    error_flag=F
    
    if(abs(upper_lim-lambda)<=1e-5){
        cat("Error!Lambda tending to infinity!Using standard\n")
        error_flag=T
    } else if(abs(low_lim-lambda)<=1e-5){
        cat("Error!Lambda tending to -infinity!Using standard\n")
        error_flag=T
    }
    
    x <- samplePeaksData
    x <- t(x)	
    N <- dim(x)[2]
    
    if(error_flag)
    {
        lambda <- lambda_std
        scal_fact <- apply(x,2,sum,na.rm=T)
        scal_fact <- mean(scal_fact)
        scal_fact <- 1/scal_fact
    }
    
    x <- x*scal_fact
    x <- x^pow_fact
    x <- x-min(x,na.rm=T)
    
    Z <- .glog(x,y0,lambda)
    
    X <- as.data.frame(t(Z)) 
    X$sample <- row.names(X)
    X <- X %>% gather("ID","intensity",-sample) 
    pData <- para@peaksData %>% select(-one_of(valueID)) %>% 
        mutate(ID=as.character(ID),sample=as.character(sample))
    
    gdata <- inner_join(pData,X,by=c("ID","sample")) %>% 
        mutate_(.dots = setNames(list(~intensity),valueID)) %>% 
        mutate(intensity=NULL)
    if(nrow(gdata)!=nrow(para@peaksData)){
        stop("Error in glog transformation!")
    }
    para@peaksData <- gdata
    return(para)
    
}

