# countReadsByBin.R
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
library(GenomicRanges)
library(Rsamtools)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test.doCount()
    
} # runTests
#------------------------------------------------------------------------------------------------------------------------
doCount <- function(bamFile, chromosome, start, end, binSize=1e6)
{
   stopifnot(file.exists(bamFile))
   starts <- seq(from=start, to=(end-binSize), by=binSize)
   ends   <- seq(from=(start + binSize), to=(end), by=binSize)
   gr <- GRanges(seqnames="chr19", ranges=IRanges(start=starts, end=ends))
   regions <- ScanBamParam(which=gr)

  indexFile <- sprintf("%s.bai", bamFile)
  if(!file.exists(indexFile))
     indexBam(bamFile)

   tbl <- countBam(bamFile, param=regions)
   return(tbl)

} # doCount
#------------------------------------------------------------------------------------------------------------------------
test.doCount <- function()
{
   tbl.counts <- doCount("../data/ENCSR000DBY.19.chr.bam", "chr19", 0, 5e6, 1e6)
   checkEquals(dim(tbl.counts), c(5, 7))
   checkEquals(colnames(tbl.counts), c("space", "start", "end", "width", "file", "records", "nucleotides"))
   checkEqualsNumeric(mean(tbl.counts$records), 18487, tol=10) # a very rough check
   checkEqualsNumeric(sd(tbl.counts$records), 2028, tol=10)  


} # test.doCount
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
