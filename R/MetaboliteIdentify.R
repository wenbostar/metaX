

##' @title Metabolite identification
##' @description Metabolite identification
##' @param para An object of \code{metaXpara},use the information from the 
##' rawPeaks to do the identification
##' @param db The file name of database
##' @param delta The delta
##' @param mode The mode of data, positive or negative,1:positive, 2:negative
##' @param unit The unit of the delta,1:ppm,2:da. Default is ppm
##' @param dbType 1=HMDB,2=PubChem,3=KEGG,4=MassBank,5=LipidMaps,6=PlantCyc
##' @return The name of output file
##' @export
metaboliteAnnotation=function(para,db,delta=10,mode=1,unit="ppm",dbType=1){
    
    rawPeaks <- para@rawPeaks
    ## whether there is a column named "mz" in file
    if(sum(names(para@rawPeaks)=="mz")==0){
        ## xcms output?
        if(sum(names(para@rawPeaks)=="mzmed")==1){
            rawPeaks <- dplyr::rename(rawPeaks,mz=mzmed)
        }else{
            stop("please check your rawPeaks data, don't find valid mz column!\n")
        }
    }
    
    if(sum(names(para@rawPeaks)=="ID")==0){
        ## xcms output?
        if(sum(names(para@rawPeaks)=="name")==1){
            rawPeaks <- dplyr::rename(rawPeaks,ID=name)
        }else{
            stop("please check your rawPeaks data, don't find valid ID column!\n")
        }
    }
    
    
    ## save peak result to a file which will be used for identification
    savePeak2file <- paste(para@outdir,"/",para@prefix,
                           "-peak2identification.txt",sep="")
    message("Save data to file: ",savePeak2file)
    write.table(rawPeaks,file=savePeak2file,col.names = TRUE,
                row.names = FALSE,quote=FALSE,sep="\t")
    
    ## begin identification
    ## get the java sofware 
    
    javabin <- system.file("tool/MetaboliteIdentify.jar",package = "metaX")
    runcmd <- paste("java -jar ",paste("\"",javabin,"\"",sep=""),
                    " -i ",savePeak2file," -db ",db,
                    " -delta ", delta, " -unit ",unit,
                    " -prefix ", para@prefix ,
                    " -dbType ", dbType ,
                    " -o ",para@outdir," -mode ",mode,sep="")

    system(runcmd)
    outfile <- paste(para@outdir,"/",para@prefix,"-identify.txt",sep="")
    return(outfile)    
}
