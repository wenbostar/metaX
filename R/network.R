


##' @title Differential correlation network analysis
##' @description Differential correlation network analysis
##' @param para A metaXpara object
##' @param valueID The name of the column that used for plot
##' @param group Samples used for plot
##' @param cor.method Method for correlation:"pearson","spearman" or "kendall"
##' @param threshold A threshold of significance levels of 
##' differential correlation
##' @param p.adjust.methods c("local", holm", "hochberg", "hommel", 
##' "bonferroni", "BH", "BY", "fdr", "none")
##' @param plot Whether to plot network figure
##' @param cluster.method The function tries to find dense subgraph. 
##' 1=cluster_fast_greedy,2=cluster_walktrap,3=cluster_edge_betweenness,
##' 4=cluster_optimal,5=cluster_leading_eigen,6=cluster_spinglass,
##' 7=cluster_label_prop,8=cluster_louvain,9=cluster_infomap
##' @param ... Additional parameter
##' @return The name of result file
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
##' resfile <- cor.network(para,group=c("S","C"))
cor.network=function(para,group,valueID="value",cor.method="spearman",
                     threshold=0.1,p.adjust.methods="BH",plot=TRUE,
                     graph_format="gml",mark.groups=TRUE,top.groups=1,
                     cluster.method=1,find.largest.component=TRUE,...){
    if(length(group)!=2){
        print(group)
        stop("We only handle two classes network!\n")
    }
    #sampleList  <- read.delim(para@sampleListFile)
    
    peaksData <- para@peaksData
    peaksData <- dplyr::filter(peaksData,class %in% group)
    xyData <- dcast(peaksData,sample+class~ID,value.var = valueID)
    xyData <- dplyr::mutate(xyData,sample=NULL)
    
    cond1 <- dplyr::filter(xyData,class %in% group[1]) %>% 
        dplyr::mutate(class=NULL) %>% t() 
    cond2 <- dplyr::filter(xyData,class %in% group[2]) %>% 
        dplyr::mutate(class=NULL) %>% t()
    

    resfile <- paste(para@outdir,"/",para@prefix,"-correlation_network.txt",
                     sep="")
    resfig <- paste(para@outdir,"/",para@prefix,"-correlation_network.pdf",
                     sep="")
    if(!is.na(pmatch(p.adjust.methods,table = "local"))){
        pdf(resfig)
    }
    
    comp.2.cc.fdr(output.file=resfile, data1 = cond1,data2 = cond2, 
                  method=cor.method,threshold=threshold,
                  p.adjust.methods = p.adjust.methods)
    
    if(!is.na(pmatch(p.adjust.methods,table = "local"))){
        dev.off()
    }
    
    #res <- read.delim(file = resfile,check.names=FALSE)
    message("Write the result to ",resfile)
    
    if(plot==TRUE){
        resfig <- paste(para@outdir,"/",para@prefix,"-",paste(group,collapse = "-"),
                        "-diff_network.pdf",sep="")
        message("Save the network figure to file: ",resfig)
        pdf(resfig,width = 5,height = 5)
        #png(resfig,width = 164,height = 164,units = "mm",res = 120)
        par(mar=c(0,0,0,0))
        a <- read.delim(resfile,stringsAsFactors = FALSE)
        gg <- graph_from_data_frame(d = a,directed = FALSE,
                                    vertices = unique(c(a$molecule.X,a$molecule.Y)))
        if(is.null(mark.groups)){
            plot(gg,vertex.size=4,vertex.label=NA,...)
        }else{
            
            if(find.largest.component==FALSE){
                if(cluster.method==1){
                    cl_net <- cluster_fast_greedy(g2)
                }else if(cluster.method==2){
                    cl_net <- cluster_walktrap(g2)
                }else if(cluster.method==3){
                    cl_net <- cluster_edge_betweenness(g2)
                }else if(cluster.method==4){
                    cl_net <- cluster_optimal(g2)
                }else if(cluster.method==5){
                    cl_net <- cluster_leading_eigen(g2)
                }else if(cluster.method==6){
                    cl_net <- cluster_spinglass(g2)
                }else if(cluster.method==7){
                    cl_net <- cluster_label_prop(g2)
                }else if(cluster.method==8){
                    cl_net <- cluster_louvain(g2)
                }else if(cluster.method==9){
                    cl_net <- cluster_infomap(g2)
                }
                
                vcol <- c("blue","red","green","orange","yellow","gray","black")
                cc <- membership(cl_net)
                cc[cc>top.groups] <- top.groups+1
                vcol[top.groups+1] <- "white"
                #ID2group <- data.frame(ID=names(cc),color = vcol[cc],NO=as.vector(cc))
                ID2group <- data.frame(ID=names(cc),
                                       color = vcol[cc],
                                       membership=as.vector(membership(cl_net)))
                centmetrics <- data.frame(ID=V(gg)$name,
                                          betweenness = betweenness(gg) ,
                                          closeness = closeness(gg), 
                                          degree = degree(gg,mode = "all"))
                ID2group <- merge(ID2group,centmetrics,by="ID")
                id2gfile <- paste(para@outdir,"/",para@prefix,"-",paste(group,collapse = "-"),
                                  "-diff_network_group.txt",sep="")
                write.table(ID2group,file=id2gfile,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
                plot(gg,vertex.label=NA,mark.groups = cl_net[1:top.groups],
                     vertex.color = vcol[cc],...)
            }else{
                
                ## find the largest component
                cc <- components(gg)
                g2 <- induced_subgraph(gg,which(cc$membership==which.max(cc$csize)))
                V(gg)[!V(gg)$name %in% V(g2)$name]$color <- "gray"
                
                ## detect communities for the largest component
                if(cluster.method==1){
                    cl_net <- cluster_fast_greedy(g2)
                }else if(cluster.method==2){
                    cl_net <- cluster_walktrap(g2)
                }else if(cluster.method==3){
                    cl_net <- cluster_edge_betweenness(g2)
                }else if(cluster.method==4){
                    cl_net <- cluster_optimal(g2)
                }else if(cluster.method==5){
                    cl_net <- cluster_leading_eigen(g2)
                }else if(cluster.method==6){
                    cl_net <- cluster_spinglass(g2)
                }else if(cluster.method==7){
                    cl_net <- cluster_label_prop(g2)
                }else if(cluster.method==8){
                    cl_net <- cluster_louvain(g2)
                }else if(cluster.method==9){
                    cl_net <- cluster_infomap(g2)
                }
                cc <- membership(cl_net)
                
                ## the number of detected communities 
                ncommutities = length(unique(membership(cl_net)))
                
                ## generate colours for plot
                vcol <- rainbow(ncommutities)
                
                ## merge information into a data.frame
                ID2group <- data.frame(ID=names(cc),
                                       color = vcol[cc],
                                       membership=as.vector(membership(cl_net)),stringsAsFactors = FALSE)
                ggname <- data.frame(ID=V(gg)$name,stringsAsFactors = FALSE)
                
                m <- left_join(ggname,ID2group,by="ID")
                m$color[is.na(m$color)] <- "gray"
                m$membership[is.na(m$membership)] <- max(cc)+1
                
                V(gg)$color <- m$color
                V(gg)$membership <-as.factor(m$membership) 
                
                
                centmetrics <- data.frame(ID=V(gg)$name,
                                          betweenness = betweenness(gg) ,
                                          closeness = closeness(gg), 
                                          degree = degree(gg,mode = "all"),
                                          stringsAsFactors = FALSE)
                m <- left_join(m,centmetrics,by="ID")
                id2gfile <- paste(para@outdir,"/",para@prefix,"-",paste(group,collapse = "-"),
                                  "-diff_network_group.txt",sep="")
                write.table(m,file=id2gfile,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
                plot(gg,...)
            }
        }
        dev.off()
        gfile <- paste(para@outdir,"/",para@prefix,"-",paste(group,collapse = "-"),
                       "-diff_network.",graph_format,sep="")
        message("Save the graph object to file: ",gfile)
        write_graph(gg,file=gfile,format = graph_format)
    }
    
    
    return(resfile)
}


##' @title Plot correlation network map
##' @description  Plot correlation network map
##' @param para A metaXpara object
##' @param valueID The name of the column that used for plot
##' @param group Samples used for plot
##' @param cor.thr Threshold of correlation
##' @param degree.thr Threshold of degree of node
##' @param size.factor Node size factor for plot 
##' @param layout layout for plotting
##' @param showPlot Whether or not to print the figure to screen
##' @param graph_format The file format for graph data to save. Currently, 
##' "edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "gml", "dot" and 
##' "leda" are supported.
##' @param ... Additional parameter
##' @return An object of igraph
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
##' gg <- plotNetwork(para,group=c("S","C"),degree.thr = 10,cor.thr = 0.8)
plotNetwork=function(para,group,valueID="value",cor.thr=0.95,degree.thr=10,
                     size.factor=0.5,layout=layout_in_circle,showPlot=FALSE,
                     graph_format="gml",...){
    peaksData <- para@peaksData
    peaksData <- dplyr::filter(peaksData,class %in% group)
    xyData <- dcast(peaksData,sample+class~ID,value.var = valueID)
    xyData <- dplyr::mutate(xyData,sample=NULL,class=NULL) %>% t() 
    
    resfig <- paste(para@outdir,"/",para@prefix,"-",paste(group,collapse = "-"),
                    "-network.png",sep="")
    message("Save the network figure to file: ",resfig)
    #pdf(resfig,width = 10,height = 10)
    png(resfig,width = 164,height = 164,units = "mm",res = 120)
    par(mar=c(0,0,0,0))
    #igraph.options(vertex.size=6)
    
    gg <- generate_g(xyData,cor.thr = cor.thr,edge.width = 1.5,node.size = 3)
    V(gg)$size=degree(gg)*size.factor
    ## check the degree of each node
    #identify those vertices part of less than three edges
    bad.vs<-V(gg)[degree(gg)<=degree.thr] 
    gg<-delete.vertices(gg, bad.vs) #exclude them from the graph
    message("V: ",vcount(gg))
    message("E: ",ecount(gg))
    plot(gg,layout=layout,vertex.label.cex=0.6,...)
    dev.off()
    message(resfig)
    if(showPlot){
        plot(gg,layout=layout,vertex.label.cex=0.6,...)
    }
    ## http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats
    gfile <- paste(para@outdir,"/",para@prefix,"-",paste(group,collapse = "-"),
                   "-network.",graph_format,sep="")
    message("Save the graph object to file: ",gfile)
    write_graph(gg,file=gfile,format = graph_format)
    return(gg)
    
}

##' @title Plot correlation heatmap
##' @description  This function plots correlation heatmap.
##' @param para A metaXpara object
##' @param valueID The name of the column that used for plot
##' @param samples Samples used for plot
##' @param label Label to show in figure
##' @param cor.method Method used for correlation
##' @param shownames A logical indicates whether show names when plot
##' @param width The width of the graphics region in inches. 
##' The default values are 6.
##' @param height The height of the graphics region in inches. 
##' The default values are 6.
##' @param anno A logical value indicates whether to plot heatmap with 
##' annotating class information
##' @param cluster A logical value indicates whether to do the cluster when 
##' anno is TRUE 
##' @param ... Additional parameter
##' @return The fig name
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
##' plotCorHeatmap(para,valueID="value",samples=NULL,width=6,anno=TRUE)
plotCorHeatmap=function(para,valueID="value",samples=NA,label="order",width=6,
                        cor.method="spearman",
                        height=6,anno=FALSE,cluster=FALSE,shownames=FALSE,...){
    peaksData <- para@peaksData
    samList <- read.delim(para@sampleListFile)
    
    if(is.null(samples)){
        message("use all the samples...")
    }else if(any(is.na(samples))){
        ## plot for QC samples
        peaksData <- dplyr::filter(peaksData,is.na(class))
        samList <- dplyr::filter(samList,is.na(class))
    }else{
        peaksData <- dplyr::filter(peaksData,class %in% samples)
        samList <- dplyr::filter(samList,class %in% samples)
    }
    
    if(label=="order"){
        xyData <- dcast(peaksData,ID~order,value.var = valueID)
        xyData <- dplyr::select(xyData,-ID)
    }else{
        xyData <- dcast(peaksData,ID~sample,value.var = valueID)
        xyData <- dplyr::select(xyData,-ID)
        
    }
    
    
    
    samList <- samList[order(samList$class,samList$order),]
    samList$order <- as.character(samList$order)
    if(label=="order"){
        diffs <- setdiff(samList$order,names(xyData))
        if(length(diffs)>=1){
            warning("The following samples are not found in the samples!")
            warning(paste(diffs,collapse = ","))
        }
        samList <- dplyr::filter(samList,!order %in% diffs)
        xyData <- xyData[,samList$order]
    }else{
        diffs <- setdiff(samList$sample,names(xyData))
        if(length(diffs)>=1){
            warning("The following samples are not found in the samples!")
            warning(paste(diffs,collapse = ","))
        }
        samList <- dplyr::filter(samList,!sample %in% diffs)
        xyData <- xyData[,samList$sample]
    }
    corres <- cor(xyData,method = cor.method,...)
    if(is.null(samples)){
        prefix = "ALL"
    }else{
        prefix <- ifelse(any(is.na(samples)),"QC",paste(samples,collapse = "_"))    
    }
    
    figpdf <- paste(para@outdir,"/",para@prefix,"-",prefix,"-corHeatmap.pdf",
                 sep="")
    figpng <- gsub(pattern = "pdf$",replacement = "png",x = figpdf)
    
    if(anno==FALSE){
        new.palette=colorRampPalette(c("black","red","yellow","white"),
                                     space="rgb")
        
        pdf(figpdf,width = width,height =height)
        gg <- levelplot(corres,col.regions=new.palette(40),cuts=30,
                        xlab="",ylab="",
                        scales=list(x=list(rot=90)))
        print(gg)
        dev.off()
        
        png(figpng,width = width,height = height,units = "in",res = 120)
        gg <- levelplot(corres,col.regions=new.palette(40),cuts=30,
                        xlab="",ylab="",
                        scales=list(x=list(rot=90)))
        print(gg)
        dev.off()
    }else{
        ## use pheatmap
        samList$class <- as.character(samList$class)
        samList$class[is.na(samList$class)] <- "QC"
        annotationList <- data.frame(Class=samList$class,
                                     Batch=as.character(samList$batch))
        if(label=="order"){
            row.names(annotationList) <- samList$order
        }else{
            row.names(annotationList) <- samList$sample
        }
        
        pdf(figpdf,width = width,height =height)
        pheatmap(corres,annotation=annotationList,border_color = NA,
                 cluster_rows=cluster,cluster_cols = cluster,
                 show_colnames = shownames,show_rownames=shownames,...)
        dev.off()
        
        png(figpng,width = width,height = height,units = "in",res = 120)
        pheatmap(corres,annotation=annotationList,border_color = NA,
                 cluster_rows=cluster,cluster_cols = cluster,
                 show_colnames = shownames,show_rownames=shownames,...)
        dev.off()
        
    }
    res <- list(fig=figpng,highfig=figpdf)
    return(res)
}
