# countReadsByBinCory.R
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

#------------------------------------------------------------------------------------------------------------------------

# my attempt at reading in the wellington footprints into an R data frame for the purpose of comparing them to the number of reads in the corresponding bam file
# should end up with a function

countFP <- function(bedFile)
{
        stopifnot(file.exists(bedFile))
        #read in the footprints
        tbl.temp <- read.table(bedFile, as.is=TRUE, sep = "\t", strip.white =TRUE)
        #combine the chr, start and stop and add the strand
        tbl.fp <- cbind.data.frame(tbl.temp$V1, tbl.temp$V2, tbl.temp$V3)
        # change colnames of dataframe
        colnames(tbl.fp) <- c("chrom", "start", "endpos")

        gr.fp <- with(tbl.fp, GRanges(seqnames=chrom, IRanges(start=start, end=endpos)))
        return(gr.fp)
}

#------------------------------------------------------------------------------------------------------------------------
# not sure how to do this. supposed to be a test to see if countFP works (it does, but the test doesn't)
test.doFootprints <- function()
{
	fp19 <- countFP("/local/Cory/resources/github/BDDS/footprints/functionalTests/output/wellington/ENCSR000DBY.19.chr.bam.chr19.peaks.400.bed.WellingtonFootprints.FDR.0.01.bed")
	
}

tbl.counts <- doCount("../data/ENCSR000DBY.19.chr.bam", "chr19", 0, 58617616, 1e5)
gr.fp <- countFP("/local/Cory/resources/github/BDDS/footprints/functionalTests/output/wellington/ENCSR000DBY.19.chr.bam.chr19.peaks.400.bed.WellingtonFootprints.FDR.0.01.bed")
gr.cs   <- with(tbl.counts, GRanges(seqnames=space, IRanges(start=start, end=end)))
tbl.overlaps <- as.data.frame(findOverlaps(gr.fp, gr.cs, type="any"))

