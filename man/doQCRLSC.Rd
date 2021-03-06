% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaX.R
\docType{methods}
\name{doQCRLSC}
\alias{doQCRLSC}
\title{Using the QC samples to do the quality control-robust spline signal 
correction}
\usage{
doQCRLSC(para, impute = TRUE, cpu = 0, ...)
}
\arguments{
\item{para}{An object of metaXpara}

\item{impute}{A logical indicates whether impute the result}

\item{cpu}{The number of cpu used for processing}

\item{...}{Additional parameters}
}
\value{
A list object
}
\description{
Using the QC samples to do the quality control-robust spline 
signal correction.
}
\details{
The smoothing parameter is optimised using leave-one-out cross validation 
to avoid overfitting.
}
\examples{
para <- new("metaXpara")
pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
rawPeaks(para) <- read.delim(pfile,check.names = FALSE)[1:20,]
sampleListFile(para) <- sfile
para <- reSetPeaksData(para)
para <- missingValueImpute(para)
res <- doQCRLSC(para,cpu=1)
}
\seealso{
\code{\link{plotQCRLSC}}
}
\author{
Bo Wen \email{wenbostar@gmail.com}
}
