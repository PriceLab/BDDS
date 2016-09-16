library(RPostgreSQL)
library(GenomicRanges)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test.combineFootprintsAndFimo()
  test.combineFootprintsAndDatabasedFimo()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
explore <- function()
{
  tbl.fimo22 <- read.table("chr22.fimo.txt", sep="\t", as.is=TRUE, nrow=-1, comment.char='',
                           row.names=NULL, header=TRUE)

  colnames(tbl.fimo22) <- c("motif", "chrom", "motif.start", "motif.end", "motif.strand", "motif.score",
                            "motif.pvalue", "motif.sequence")
  tbl.fimo22$chrom <- paste("chr", tbl.fimo22$chrom, sep="")
  
  tbl.hint <- read.table("ENCSR000DBY.hint.bed", sep="\t", as.is=TRUE, nrow=-1)[, c(1,2,3,5)]
  colnames(tbl.hint) <- c("chrom", "footprint.start", "footprint.end", "footprint.score")
  
  tbl.hint.22 <- subset(tbl.hint, chrom=="chr22")
  
  gr.fimo <- with(tbl.fimo22, GRanges(seqnames=chrom, IRanges(start=motif.start, end=motif.end)))
  gr.hint <- with(tbl.hint.22,   GRanges(seqnames=chrom, IRanges(start=footprint.start, end=footprint.end)))
  tbl.overlaps <- as.data.frame(findOverlaps(gr.fimo, gr.hint))  # 83280 x 2
  
  # identify a range of hits (footprint indices) where there are some hits, some misses
  # footprints 25-28 inclusive include two footprints without motifs (26, 28) and
  # two footprints with multiple motifs (2 for #25, 4 for #7)
  #
  # subset(tbl.overlaps, subjectHits >= 25 & subjectHits <= 28)
  #
  #       queryHits subjectHits
  # 15984   1090321          27
  # 17312   1188165          27
  # 32431   2404130          27
  # 42068   3210467          25
  # 42069   3210468          25
  # 69735   5508105          27
  
  tbl.hint.22.sub <- tbl.hint.22[25:28,]   # 4 footprints
  tbl.overlaps.sub <- subset(tbl.overlaps, subjectHits >= 25 & subjectHits <= 28)   # 6 hits in two footprints
  length(tbl.overlaps.sub$subjectHits) 
  # in 4 footprints.  2 footprints have no motif. a good case study.  how to merge?
  as.data.frame(table(tbl.overlaps.sub$subjectHits))
  #  Var1 Freq
  # 1   25    2
  # 2   27    4
  
  tbl.empty.fimo <- data.frame(motif=NA_character_,
                               motif.start=NA_integer_,
                               motif.end=NA_integer_,
                               motif.strand=NA_character_,
                               motif.score=NA_real_,
                               motif.pvalue=NA_real_,
                               motif.sequence=NA_character_,
                               motif.footprint.overlap=NA_integer_,
                               stringsAsFactors=FALSE)
  
  fimo.indices <- tbl.overlaps.sub$queryHits
  hint.indices <- tbl.overlaps.sub$subjectHits
  hint.withMotifs.indices <- unique(tbl.overlaps.sub$subjectHits)
  hint.withoutMotifs.indices <- setdiff(25:28, hint.withMotifs.indices)
  
  tbl.mfp <- cbind(tbl.hint.22[hint.indices,], tbl.fimo22[fimo.indices, -grep("chrom", colnames(tbl.fimo22))])
  tbl.mfp$motif.footprint.overlap <- unlist(lapply(1:nrow(tbl.mfp),
      function(row) length(intersect(tbl.mfp$footprint.start[row]:tbl.mfp$footprint.end[row],
                                      tbl.mfp$motif.start[row]:tbl.mfp$motif.end[row]))))
  
  tbl.hintNoMotifs <- tbl.hint.22[hint.withoutMotifs.indices,]
  tbl.hintNoMotifs <- cbind(tbl.hintNoMotifs, tbl.empty.fimo)
  
  tbl.out <- rbind(tbl.mfp, tbl.hintNoMotifs)
  tbl.out <- tbl.out[order(tbl.out$footprint.start),]
  row.names(tbl.out) <- NULL

}
#------------------------------------------------------------------------------------------------------------------------
combineFootprintsAndFimo <- function(tbl.footprints, tbl.fimoMotifs, chromosome)
{
   stopifnot(colnames(tbl.footprints) == c("chrom", "footprint.start", "footprint.end", "footprint.score"))
   stopifnot(all (c("motif", "chrom", "motif.start", "motif.end", "motif.strand", "motif.score",
                    "motif.pvalue", "motif.sequence") %in% colnames(tbl.fimoMotifs)))

   tbl.fp <- subset(tbl.footprints, chrom==chromosome)
   stopifnot(nrow(tbl.fp) > 0)
   
   tbl.fimo <- subset(tbl.fimoMotifs, chrom==chromosome)
   stopifnot(nrow(tbl.fimo) > 0)

   gr.fimo <- with(tbl.fimo,  GRanges(seqnames=chrom, IRanges(start=motif.start, end=motif.end)))
   gr.fp <- with(tbl.fp,  GRanges(seqnames=chrom, IRanges(start=footprint.start, end=footprint.end)))
   tbl.overlaps <- as.data.frame(findOverlaps(gr.fimo, gr.fp))  # 83280 x 2

   tbl.empty.fimo <- data.frame(motif=NA_character_,
                                motif.start=NA_integer_,
                                motif.end=NA_integer_,
                                motif.strand=NA_character_,
                                motif.score=NA_real_,
                                motif.pvalue=NA_real_,
                                motif.sequence=NA_character_,
                                motif.footprint.overlap=NA_integer_,
                                stringsAsFactors=FALSE)
   
   fimo.indices <- tbl.overlaps$queryHits
   fp.indices <- tbl.overlaps$subjectHits

   fp.withMotifs.indices <- unique(tbl.overlaps$subjectHits)
   fp.withoutMotifs.indices <- setdiff(1:nrow(tbl.fp), fp.withMotifs.indices)
   stopifnot(length(fp.withMotifs.indices) + length(fp.withoutMotifs.indices) == nrow(tbl.fp))
   
   tbl.mfp <- cbind(tbl.fp[fp.indices,], tbl.fimo[fimo.indices, -grep("chrom", colnames(tbl.fimo))])

   printf("--- calculating motif.footprint.overlap")
   tbl.mfp$motif.footprint.overlap <- unlist(lapply(1:nrow(tbl.mfp),
             function(row) length(intersect(tbl.mfp$footprint.start[row]:tbl.mfp$footprint.end[row],
                                  tbl.mfp$motif.start[row]:tbl.mfp$motif.end[row]))))

   tbl.fpNoMotifs <- tbl.fp[fp.withoutMotifs.indices,]
   tbl.fpNoMotifs <- cbind(tbl.fpNoMotifs, tbl.empty.fimo)

   tbl.out <- rbind(tbl.mfp, tbl.fpNoMotifs)
   tbl.out <- tbl.out[order(tbl.out$footprint.start),]
   tbl.out$motif.footprint.overlap[is.na(tbl.out$motif.start)] <- NA_integer_
   row.names(tbl.out) <- NULL

   invisible(tbl.out)
   
} # combineFootprintsAndFimo
#------------------------------------------------------------------------------------------------------------------------
test.combineFootprintsAndFimo <- function()
{
   printf("--- test.combineFootprintsAndFimo")

   if(!exists("tbl.fimo22")){
      if(file.exists("test.RData")){
          load("test.RData", envir=.GlobalEnv)
      } else {
        tbl.fimo22 <- read.table("chr22.fimo.txt", sep="\t", as.is=TRUE, nrow=-1, comment.char='',
                                 row.names=NULL, header=TRUE)
        colnames(tbl.fimo22) <- c("motif", "chrom", "motif.start", "motif.end", "motif.strand", "motif.score",
                                  "motif.pvalue", "motif.sequence")
        tbl.fimo22$chrom <- paste("chr", tbl.fimo22$chrom, sep="")
        tbl.hint <- read.table("ENCSR000DBY.hint.bed", sep="\t", as.is=TRUE, nrow=-1)[, c(1,2,3,5)]
        colnames(tbl.hint) <- c("chrom", "footprint.start", "footprint.end", "footprint.score")
        tbl.hint.22 <- subset(tbl.hint, chrom=="chr22")
        save(tbl.fimo22, tbl.hint.22, file="test.RData")
        } # else: read in data afresh, save serialized
      } # variable not found: load serialized, or create from scratch

   tbl.out <- combineFootprintsAndFimo(tbl.hint.22[1:10, ], tbl.fimo22, "chr22")
   checkEquals(colnames(tbl.out), c("chrom", "footprint.start", "footprint.end", "footprint.score", 
                                    "motif", "motif.start", "motif.end", "motif.strand", 
                                    "motif.score", "motif.pvalue", "motif.sequence", "motif.footprint.overlap"))

   checkEquals(dim(tbl.out), c(36, 12))
   checkEquals(length(unique(tbl.out$footprint.start)), 10)
    
} # test.combineFootprintsAndFimo
#------------------------------------------------------------------------------------------------------------------------
run.combineFootprintsAndFimo <- function()
{
   printf("--- run.combineFootprintsAndFimo")

   if(!exists("tbl.fimo22")){
      if(file.exists("test.RData")){
          load("test.RData", envir=.GlobalEnv)
      } else {
        tbl.fimo22 <- read.table("chr22.fimo.txt", sep="\t", as.is=TRUE, nrow=-1, comment.char='',
                                 row.names=NULL, header=TRUE)
        colnames(tbl.fimo22) <- c("motif", "chrom", "motif.start", "motif.end", "motif.strand", "motif.score",
                                  "motif.pvalue", "motif.sequence")
        tbl.fimo22$chrom <- paste("chr", tbl.fimo22$chrom, sep="")
        tbl.hint <- read.table("ENCSR000DBY.hint.bed", sep="\t", as.is=TRUE, nrow=-1)[, c(1,2,3,5)]
        colnames(tbl.hint) <- c("chrom", "footprint.start", "footprint.end", "footprint.score")
        tbl.hint.22 <- subset(tbl.hint, chrom=="chr22")
        save(tbl.fimo22, tbl.hint.22, file="test.RData")
        } # else: read in data afresh, save serialized
      } # variable not found: load serialized, or create from scratch

   tbl.out <- combineFootprintsAndFimo(tbl.hint.22, tbl.fimo22, "chr22")
   checkEquals(colnames(tbl.out), c("chrom", "footprint.start", "footprint.end", "footprint.score", 
                                    "motif", "motif.start", "motif.end", "motif.strand", 
                                    "motif.score", "motif.pvalue", "motif.sequence", "motif.footprint.overlap"))

   checkEquals(dim(tbl.out), c(85509, 12))
   checkEquals(length(unique(tbl.out$footprint.start)), nrow(tbl.hint.22))
    
} # run.combineFootprintsAndFimo
#------------------------------------------------------------------------------------------------------------------------
test.combineFootprintsAndDatabasedFimo <- function()
{
   printf("--- test.combineFootprintsAndDatabasedFimo")
   if(!exists("db"))
      db <<- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="fimo", host="whovian")

   tbl.fimo22 <- dbGetQuery(db, "select * from hg38 where chr='22'")
   colnames(tbl.fimo22) <- c("motif", "chrom", "motif.start", "motif.end", "motif.strand", "motif.score",
                             "motif.pvalue", "empty", "motif.sequence")
   tbl.fimo22 <- tbl.fimo22[, -grep("empty", colnames(tbl.fimo22))]
   tbl.fimo22$chrom <- paste("chr", tbl.fimo22$chrom, sep="")
   tbl.hint <- read.table("ENCSR000DBY.hint.bed", sep="\t", as.is=TRUE, nrow=-1)[, c(1,2,3,5)]
   colnames(tbl.hint) <- c("chrom", "footprint.start", "footprint.end", "footprint.score")
   tbl.hint.22 <- subset(tbl.hint, chrom=="chr22")
   save(tbl.fimo22, tbl.hint.22, file="test.RData")

   tbl.out <- combineFootprintsAndFimo(tbl.hint.22[1:10, ], tbl.fimo22, "chr22")
   checkEquals(colnames(tbl.out), c("chrom", "footprint.start", "footprint.end", "footprint.score", 
                                    "motif", "motif.start", "motif.end", "motif.strand", 
                                    "motif.score", "motif.pvalue", "motif.sequence", "motif.footprint.overlap"))
   checkEquals(dim(tbl.out), c(36, 12))
   checkEquals(length(unique(tbl.out$footprint.start)), 10)

} # test.combineFootprintsAndDatabasedFimo
#------------------------------------------------------------------------------------------------------------------------
