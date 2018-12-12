

# "pearson", "kendall", "spearman"
# This function is used to calculate the correlation between protein and mRNA.
calcCorBetweenProteinAndRNA=function(para1,para2,log=TRUE,select_by=1,top_n=1000,
                                     geneset = NULL,
                                     outdir = "./",prefix = "test",use_class=NULL,
                                     reg=FALSE,
                                     cor_method="spearman",valueID="value"){
    
    if(!is.null(use_class)){
        cat("Only use class ",use_class,"\n")
        if(reg == FALSE){
            para1@peaksData <- para1@peaksData %>% dplyr::filter(class==use_class)
            para2@peaksData <- para2@peaksData %>% dplyr::filter(class==use_class)
        }else{
            para1@peaksData <- para1@peaksData %>% dplyr::filter(grepl(pattern = use_class,x = class))
            para2@peaksData <- para2@peaksData %>% dplyr::filter(grepl(pattern = use_class,x = class))
        }
        cat("Used samples:",length(unique(para1@peaksData$sample)),length(unique(para2@peaksData$sample)),"\n")
    }else{
        cat("Use all samples!\n")
    }
    
    if(!is.null(geneset)){
        cat("Only use specified gene set\n")
        para1@peaksData <- para1@peaksData %>% dplyr::filter(ID %in% geneset)
        para2@peaksData <- para2@peaksData %>% dplyr::filter(ID %in% geneset)  
    }
    
    x1 <- para1@peaksData %>% dplyr::select(ID,sample,!!valueID)
    x2 <- para2@peaksData %>% dplyr::select(ID,sample,!!valueID)
    m <- merge(x1,x2,by=c("ID","sample"))
    
    ## firstly, plot a scatterplot: x axis is sd, y axis is cor
    names(m)[3:4] <- c("x","y")
    m$x[m$x <= 0] <- NA
    m$y[m$y <= 0] <- NA
    if(log==TRUE){
        m$x <- log2(m$x)    
        m$y <- log2(m$y)
    }
    
    ## filter genes with all values are missing
    na_row <- apply(m[,3:4],1,function(x){any(is.na(x))})
    m <- m[!na_row,]
    
    ## Sample wise correlation
    save(m,file = "m.rda")
    sample_wise_cor <- m %>% group_by(sample) %>% 
        dplyr::summarise(cor=cor(x,y,use = "com",method = cor_method)) %>%
        dplyr::mutate(label=prefix)
    
    ## Gene wise correlation
    res <- m %>% group_by(ID) %>% dplyr::summarise(sd_x=sd(x,na.rm = TRUE),
                                            sd_y=sd(y,na.rm=TRUE),
                                            cor=cor(x,y,use = "com",method = cor_method))
    #save(para1,para2,valueID,na_row,m,res,file="test111.rda")
    
    png(paste(outdir,"/",prefix,"-cor-sd.png",sep=""),width = 600,height = 1200,res=120)
    cor_name <- paste("Cor (",cor_method,")",sep="")
    par(mfrow=c(4,2),mar=c(3,3,2,1),mgp=c(1.6,0.6,0))
    plot(res$sd_x,res$cor, xlab="SD (1)",ylab=cor_name,col=rgb(0.5,0.4,0.7,0.4),pch=15,cex=0.6)
    plot(res$sd_y,res$cor, xlab="SD (2)",ylab=cor_name,col=rgb(0.5,0.4,0.7,0.4),pch=15,cex=0.6)
    plot(res$sd_x,res$sd_y, xlab="SD (1)",ylab="SD (2)",col=rgb(0.5,0.4,0.7,0.4),pch=15,cex=0.6)
    hist(res$sd_x,nclass=50,xlab="SD (1)",main="")
    hist(res$sd_y,nclass=50,xlab="SD (2)",main="")
    cor_mean <- sprintf("%.4f",mean(res$cor,na.rm=TRUE))
    cor_md <- sprintf("%.4f",median(res$cor,na.rm=TRUE))
    hist(res$cor,nclass=50,xlab=cor_name,main=paste("Data point: ",nrow(res),"\nmean = ",cor_mean,", median = ",cor_md,sep = ""))
    dev.off()
    
    rr <- data.frame(n=nrow(res),
               n5=sum(res$cor>=0.5,na.rm = TRUE),
               n6=sum(res$cor>=0.6,na.rm = TRUE),
               n7=sum(res$cor>=0.7,na.rm = TRUE),
               n8=sum(res$cor>=0.8,na.rm = TRUE))
    
    fres <- list()
    fres$feature_wise <- rr
    fres$sample_wise <- sample_wise_cor
    return(fres)
}









