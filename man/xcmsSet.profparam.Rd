% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaX.R
\docType{methods}
\name{xcmsSet.profparam<-}
\alias{xcmsSet.profparam<-}
\title{xcmsSet.profparam}
\usage{
xcmsSet.profparam(para) <- value
}
\arguments{
\item{para}{An object of metaXpara}

\item{value}{value}
}
\value{
An object of metaXpara
}
\description{
xcmsSet.profparam
}
\examples{
library(faahKO)
para <- new("metaXpara")
xcmsSet.profparam(para) <- list(step=0.005)
}
\author{
Bo Wen \email{wenbostar@gmail.com}
}
