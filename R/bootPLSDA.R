
##' @title Fit predictive models for PLS-DA
##' @description Fit predictive models for PLS-DA
##' @param x An object where samples are in rows and features are in columns. 
##' This could be a simple matrix, data frame. 
##' @param y A numeric or factor vector containing the outcome for each sample.
##' @param ncomp The maximal number of component for PLS-DA
##' @param sample A vector contains the sample used for the model
##' @param test The data set (data.frame) for testing. If the data contains a 
##' column with the name "class", this column is the sample class. 
##' @param split Whether split the data as train and test set. Default is 0 
##' which indicates not split the data. 
##' @param method The resampling method: boot, boot632, cv, repeatedcv, LOOCV, 
##' LGOCV (for repeated training/test splits), none (only fits one model to the 
##' entire training set), oob (only for random forest, bagged trees, bagged 
##' earth, bagged flexible discriminant analysis, or conditional tree forest 
##' models), "adaptive_cv", "adaptive_boot" or "adaptive_LGOCV"
##' @param repeats For repeated k-fold cross-validation only: the number of 
##' complete sets of folds to compute
##' @param number Either the number of folds or number of resampling iterations
##' @param ... Arguments passed to the classification or regression routine
##' @return A list object
bootPLSDA=function(x,y,ncomp=2,sample=NULL,test=NULL,split=0,
                   method = "repeatedcv",repeats=250,number=7,...){

    ## x: each row is a sample, each column is a peak
    ## y: a vector of response
    # dataset <- cbind(y,x)
    ## dataset: The data as a vector, matrix or data frame. If it is a matrix or 
    ## data frame then each row is considered as one multivariate observation.
    #bootres <- boot(data = dataset,statistic = mybootPLSDA,R = 1000)
    
    ##
    ptm <- proc.time()
    
    rp <- sample
    if(is.numeric(y)){
        y <- paste("C",y,sep="")
        rp <- paste("C",sample,sep="")
    }
    y <- as.character(y)
    
    ind <- y %in% rp
    x <- x[ind,]
    if(is.matrix(y) | is.data.frame(y)){
        y <- y[ind,]    
    }else{
        y <- y[ind]
    }
    y <- factor(x = as.character(y),levels = rp)
    
    
    ## result data to return
    res <- list()
    
    ctrl <- trainControl(method = method,
                         repeats = repeats,
                         number = number,
                         classProbs = TRUE,
                         summaryFunction = twoClassSummary)
    
    
    ## split samples for train and testing
    ## Firstly, if test is not NULL, then use the data in test for testing
    if(!is.null(test)){
        ## the first column of test must be class
        res$test <- test
        
        res$xTrain <- x
        res$yTrain <- y
        
        res$xTest  <- test[,-1]
        res$yTest  <- test[,1]
        
        
    }else if(split!=0){
        res$trainIndex <- createDataPartition(y, p = split,list = FALSE)
        res$xTrain <- x[trainIndex,]
        res$yTrain <- y[trainIndex]
        
        res$xTest  <- x[-trainIndex,]
        res$yTest <- y[-trainIndex]
        
    }else{
        res$xTrain <- x
        res$yTrain <- y
        
    }
    
    plsFit <- train(res$xTrain,res$yTrain,
                    method = "pls",
                    tuneLength = ncomp,
                    trControl = ctrl,
                    #validation = "CV",
                    #segments=7,
                    metric = "ROC",...)
    
    ## need to predict
    if(!is.null(test) || split!=0){
        ## predict
        trueClass <- res$yTest
        if(is.numeric(trueClass)){
            trueClass <- paste("C",trueClass,sep="")
        }
        trueClass <- as.character(trueClass)
        trueClass <- factor(x = as.character(trueClass),levels = rp)
        
        ## confusionMatrix
        predictResConfm <- predict(plsFit, newdata = res$xTest)
        predictResProb  <- predict(plsFit, newdata = res$xTest, type = "prob")
        rocRes <- pROC::roc(trueClass,predictResProb[,1],levels = rev(rp))
        par(mar=c(3,3,1,1),mgp=c(1.6,0.6,0))
        plot(rocRes, type = "S", print.thres = .5)
        
        res$predictResConfm <- predictResConfm
        res$predictResProb <- predictResProb
        res$rocRes <- rocRes
        res$confusionMatrix <- confusionMatrix(predictResConfm, trueClass)
        print(res$confusionMatrix)
        
    }
    
    
    res$model <- plsFit
    res$x <- x
    res$y <- y
    res$time <- proc.time() - ptm
    return(res)
    
}


##' @title Create predictive models
##' @description Create predictive models
##' @param para An object of \code{metaXpara} 
##' @param method Method for model construction
##' @param group Sample class used
##' @param valueID The name of column used
##' @param ... Additional arguments.
##' @return A list object
##' @examples 
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
##' sampleListFile(para) <- sfile
##' para <- reSetPeaksData(para)
##' para <- missingValueImpute(para)
##' para <- transformation(para,method = 1,valueID = "value")
##' para <- metaX::preProcess(para,scale = "uv",center = TRUE,
##'                         valueID = "value")
##' rs <- createModels(para,method="plsda",group=c("S","C"),valueID="value")
createModels=function(para,method="plsda",group=NA,valueID="value",...){
    
    message("value = ",valueID)
    peaksData <- para@peaksData
    pdata <- dcast(peaksData,sample+class~ID,value.var = valueID)
    pdata$sample <- NULL
    
    if(method=="plsda"){
        res <- bootPLSDA(x = pdata[,-1],y=pdata[,1],sample = group,...)
        return(res) 
    }
    
}

