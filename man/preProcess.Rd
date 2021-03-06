% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep.R
\name{preProcess}
\alias{preProcess}
\title{Pre-Processing}
\usage{
preProcess(para, t = 0, scale = c("auto", "range", "pareto", "vast",
  "level", "power", "none"), center = TRUE, valueID = "valueNorm")
}
\arguments{
\item{para}{An metaX object}

\item{t}{The method for transformation, 0=none, 1=log, 2=Cube root, 3=glog}

\item{scale}{The method of scaling: "auto", "range", "pareto", "vast", "level", "power","none"}

\item{center}{Centering}

\item{valueID}{The name of column used for transformation}

\item{...}{Additional parameter}
}
\value{
An new metaX object
}
\description{
Pre-Processing
}
\examples{
para <- new("metaXpara")
pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
sampleListFile(para) <- sfile
para <- reSetPeaksData(para)
para <- missingValueImpute(para)
para <- preProcess(para,valueID = "value",scale="uv")
}
\author{
Bo Wen \email{wenbostar@gmail.com}
}
