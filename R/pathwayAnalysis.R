
##' @title Pathway analysis
##' @description Pathway analysis
##' @param id A vector of metabolite IDs
##' @param id.type The type of metabolite ID type, default is hmdb.
##' @param outfile The output file name
##' @return A data.frame object
##' @export
##' @examples 
##' \dontrun{
##' res <- pathwayAnalysis(id=c("HMDB00060","HMDB00056","HMDB00064"),
##'     outfile="pathway.csv")
##' head(res)
##' }
pathwayAnalysis=function(id,id.type="hmdb",outfile){

    myHttpheader<- c(
        "Connection"="keep-alive",
        "Host"="impala.molgen.mpg.de",
        "User-Agent"="Mozilla/5.0 (Windows NT 6.2; WOW64; rv:34.0) Gecko/20100101 Firefox/34.0",
        "Accept"="text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "Accept-Language"="zh-cn,zh;q=0.8,en-us;q=0.5,en;q=0.3",
        "Content-Type"="multipart/form-data; boundary=----WebKitFormBoundaryAt5WkZAJ0FIYcbhn",
        #"Accept-Encoding"="gzip,deflate"#,
        "Cookie"="_ZopeId=57367910A65do644nmY",
        "Referer"="http://impala.molgen.mpg.de/impala/"
    )
    d = debugGatherer()
    cHandle2 <- getCurlHandle(httpheader=myHttpheader,followlocation=1,
                              debugfunction=d$update,verbose=TRUE)
    postForm("http://impala.molgen.mpg.de/impala/impala/runAnalysis",
             geneacctype="invalid",
             compounds=paste(id,collapse = "\n"),
             compoundacctype=id.type,
             analysis="ora",
             curl = cHandle2,style = "httppost")
    a <- getURL("http://impala.molgen.mpg.de/impala/impala/downloadResults",
                curl=cHandle2 ,.encoding = 'UTF-8',.opts=list())
    message("Save result to file: ",outfile)
    write(a,file=outfile)
    res <- read.csv(outfile)
    return(res)
}