test_hasQC=function(){
    para <- new("metaXpara")
	pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
	sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
	para@rawPeaks <- read.delim(pfile,check.names = FALSE)
	para@sampleListFile <- sfile
	para <- reSetPeaksData(para)
	checkEquals(hasQC(para),TRUE)
}
