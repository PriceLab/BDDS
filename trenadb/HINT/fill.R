library(RPostgreSQL)
library(GenomicRanges)
library(RUnit)
source("../regionAndHitsSchemas.R")
#------------------------------------------------------------------------------------------------------------------------
hint.path <- "/local/Cory/for_Paul/hint_sample"
#------------------------------------------------------------------------------------------------------------------------
if(!exists("db.hint"))
   db.hint <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hint", host="whovian")

if(!exists("db.fimo"))
   db.fimo <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="trena", host="whovian")

knownLocs <- new.env(parent=emptyenv())
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test.readHintTable()
  test.mergeFootprintsWithFimo()
  test.splitTableIntoRegionsAndHintHits()
  test.fill.to.database()
    
  #test.combineFootprintsAndFimo()
  #test.combineFootprintsAndDatabasedFimo()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
readHintTable <- function(directory, filename, nrows=NA, chromosome=NA)
{
   full.path <- file.path(directory, filename)
 
   if(!file.exists(full.path))
       return(data.frame)
   
   tbl <- read.table(full.path, sep="\t", as.is=TRUE)[, 1:5]
   colnames(tbl) <- c("chrom", "footprint.start", "footprint.end", "loc", "score")
   if(!is.na(chromosome))
      tbl <- subset(tbl, chrom==chromosome)

   if(!is.na(nrows))
      tbl <- tbl[1:nrows,]

   invisible(tbl)

} # readHintTable
#------------------------------------------------------------------------------------------------------------------------
test.readHintTable <- function()
{
   printf("--- test.readHintTable")

   tbl <- readHintTable(hint.path, "ENCSR000EMR.hint.bed", 5, "chr21")
   checkEquals(dim(tbl), c(5,5))
   checkEquals(colnames(tbl), c("chrom", "footprint.start", "footprint.end", "loc", "score"))
   checkEquals(unique(tbl$chrom), "chr21")

     # now read without chrom or nrow retraints
   tbl <- readHintTable(hint.path, "ENCSR000EMR.hint.bed")
   checkEquals(ncol(tbl), 5)
   checkTrue(nrow(tbl) > 600000)
   checkEquals(colnames(tbl), c("chrom", "footprint.start", "footprint.end", "loc", "score"))
   checkEquals(head(sort(unique(tbl$chrom))), c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14"))
   invisible(tbl)
   
} # test.readHintTable
#------------------------------------------------------------------------------------------------------------------------
mergeFimoWithFootprints <- function(tbl.fp, sampleID)
{
  chromosome <- unique(tbl.fp$chrom)
  stopifnot(length(chromosome) == 1)
  chromosome <- sub("chr", "", chromosome)  # fimo database lists chromosome without 'chr' prefix
  min.pos <- min(tbl.fp$footprint.start)
  max.pos <- max(tbl.fp$footprint.end)
  
  query <- sprintf("select * from fimo_hg38 where chrom='%s' and start >= %d and endpos <= %d",
                   chromosome, min.pos, max.pos)

  tbl.fimo <- dbGetQuery(db.fimo, query)
  colnames(tbl.fimo) <- c("motif", "chrom", "motif.start", "motif.end", "motif.strand", "fimo.score",
                          "fimo.pvalue", "empty", "motif.sequence")
  tbl.fimo <- tbl.fimo[, -grep("empty", colnames(tbl.fimo))]
  tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")

  gr.fimo <- with(tbl.fimo, GRanges(seqnames=chrom, IRanges(start=motif.start, end=motif.end)))

    # --- get some footprints

  gr.hint <- with(tbl.fp,   GRanges(seqnames=chrom, IRanges(start=footprint.start, end=footprint.end)))
  tbl.overlaps <- as.data.frame(findOverlaps(gr.fimo, gr.hint, type="within"))

  tbl.fimo$loc <- with(tbl.fimo, sprintf("%s:%d-%d", chrom, motif.start, motif.end))
  tbl.fimo$method <- "HINT"
  tbl.fimo$sample_id <- sampleID
  tbl.regions <- tbl.fimo[tbl.overlaps$queryHits,]

  tbl.regions <- cbind(tbl.regions, hint.score=tbl.fp[tbl.overlaps$subjectHits, "score"])
  invisible(tbl.regions)

} # mergeFimoWithFootprints
#------------------------------------------------------------------------------------------------------------------------
test.mergeFootprintsWithFimo <- function()
{
   printf("--- test.mergeFootprintsWithFimo")
   tbl.fp <- readHintTable(hint.path, "ENCSR000EMR.hint.bed")
   tbl.fp <- head(subset(tbl.fp, chrom=="chr21"), n=3)
   tbl <- mergeFimoWithFootprints(tbl.fp, "ENCSR00EMR")
   checkEquals(ncol(tbl), 12)
   checkEquals(sort(colnames(tbl)),
               c("chrom", "fimo.pvalue", "fimo.score", "hint.score", "loc", "method",
                 "motif", "motif.end", "motif.sequence", "motif.start", "motif.strand", "sample_id"))
   checkTrue(nrow(tbl) >= 9)
   checkEquals(length(grep("chr21:5014930-5014939", tbl$loc)), 5)
     # 5 distinct motifs mapped to these region
   checkEquals(sort(subset(tbl, loc == "chr21:5014930-5014939")$motif),
               c("MA0522.2", "MA0820.1", "MA0824.1", "MA0830.1", "MESP1_DBD"))
   invisible(tbl)

} # test.mergeFootprintsWithFimo
#------------------------------------------------------------------------------------------------------------------------
splitTableIntoRegionsAndHintHits <- function(tbl, minid)
{
   tbl.regions <- unique(tbl[, c("loc", "chrom", "motif.start", "motif.end")])
   colnames(tbl.regions) <- region.schema() # 29

   new.locs <- setdiff(tbl.regions$loc, names(knownLocs))
     # enter these new.locs into the hash
   lapply(new.locs, function(loc) knownLocs[[loc]] <- 0)   
   printf("novel locs: %d new/%d possible (%d now known, %d dups)",
          length(new.locs), nrow(tbl.regions), length(names(knownLocs)), nrow(tbl.regions) - length(new.locs))

   if(length(new.locs) > 0){
      keepers <- match(new.locs, tbl.regions$loc)
      tbl.regions <- tbl.regions[keepers,]
    } else {
       tbl.regions <- data.frame()
       }

      # now start working on the hits subtable.  
      #   add length.
      #   extract columns of interest.
      #      loc -> loc
      #      motif -> name
      #      motif.strand -> strand
      #      hint.score -> score1
      #      fimo.score -> score2
      #      fimo.pvalue -> score3
      #      
      #   make some table-wide assignments(provenance(minid), unused scores, 
      #                type(motif.in.footprint), method
      #   
      # hit.schema
      #  [1] "loc"        "type"       "name"       "length"     "strand"    
      #  [6] "sample_id"  "method"     "provenance" "score1"     "score2"    
      # [11] "score3"     "score4"     "score5"     "score6"    

   tbl$length <- with(tbl, 1 + motif.end - motif.start)
   tbl$provenance <- minid
   tbl$score4 <- NA
   tbl$score5 <- NA
   tbl$score6 <- NA
   tbl$type <- "motif.in.footprint"
   tbl.hits <- tbl[, c("loc", "type", "motif", "length", "motif.strand", "sample_id", 
                       "method", "provenance", "hint.score", "fimo.score", "fimo.pvalue",
                       "score4", "score5", "score6")]

   colnames(tbl.hits) <- hit.schema()
   invisible(list(regions=tbl.regions, hits=tbl.hits))
    
} # splitTableIntoRegionsAndHintHits    
#------------------------------------------------------------------------------------------------------------------------
test.splitTableIntoRegionsAndHintHits <- function(tbl)
{
   printf("--- test.splitTableIntoRegionsAndHintHits")
   tbl <- test.mergeFootprintsWithFimo()
   x <- splitTableIntoRegionsAndHintHits(tbl, "minid")
   checkEquals(sort(names(x)), c("hits", "regions"))
   checkEquals(colnames(x$regions), region.schema()) #
   checkEquals(colnames(x$hits), hit.schema())

} # test.splitTableIntoRegionsAndHintHits    
#------------------------------------------------------------------------------------------------------------------------
fill.to.database <- function(tbl.regions, tbl.hits)
{
   appendToRegionsTable(tbl.regions)
   appendToHitsTable(tbl.hits)
    
} # fill.to.database
#------------------------------------------------------------------------------------------------------------------------
test.fill.to.database <- function()
{
   printf("--- test.fill.to.database")
   createEmptyDatabaseTables()
   knownLocs <<- new.env(parent=emptyenv())

   tbl <- test.mergeFootprintsWithFimo()
   x <- splitTableIntoRegionsAndHintHits(tbl, "minid")

   fill.to.database(x$regions, x$hits)
   checkEquals(sort(dbListTables(db.hint)), c("hits", "regions"))

   tbl.regions <- dbGetQuery(db.hint, "select * from regions")
   checkEquals(dim(tbl.regions), c(5, 4))
   checkTrue(all(region.schema() == colnames(tbl.regions)))

   tbl.hits <- dbGetQuery(db.hint, "select * from hits")
   checkTrue(all(hit.schema() == colnames(tbl.hits)))
   checkEquals(dim(tbl.hits), c(9, 13))
   checkEquals(unique(tbl.hits$type), "motif.in.footprint")
   checkTrue(all(tbl.hits$loc %in% tbl.regions$loc))

} # test.fill.to.database
#------------------------------------------------------------------------------------------------------------------------
databaseSummary <- function()
{
    region.count <- dbGetQuery(db.hint, "select count(*) from regions")[1,1]
    hit.count <- dbGetQuery(db.hint, "select count(*) from hits")[1,1]
    printf("%d hits in %d regions", hit.count, region.count)

} # databaseSummary
#------------------------------------------------------------------------------------------------------------------------
createEmptyDatabaseTables <- function()
{
   system("/bin/psql -f createTables.sql")
    
} # appendToRegionsTable
#------------------------------------------------------------------------------------------------------------------------
appendToRegionsTable <- function(tbl)
{
   write.table(tbl, file="regions.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", na="NULL")
   system("/bin/psql -f fillRegions.sql")
    
} # appendToRegionsTable
#------------------------------------------------------------------------------------------------------------------------
appendToHitsTable <- function(tbl)
{
   write.table(tbl, file="hits.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", na="NULL")
   system("/bin/psql -f fillHits.sql")

} # appendToHitsTable
#------------------------------------------------------------------------------------------------------------------------
fill.all.samples.by.chromosome <- function(chromosome)
{
   # createEmptyDatabaseTables()
   knownLocs <<- new.env(parent=emptyenv())

   all.files <- list.files(hint.path, pattern=".*hint.bed$")
   minid <- "corys.hint.minid"

   for(file in all.files){
      printf("---- %s (%s) (%d/%d)", file, chromosome, grep(file, all.files), length(all.files))
      sample.id <- sub(".hint.bed", "", file)
      tbl.hint <- readHintTable(hint.path, file, NA, chromosome)
      tbl <- mergeFimoWithFootprints(tbl.hint, sample.id)
      x <- splitTableIntoRegionsAndHintHits(tbl, minid)
      printf("filling %d regions, %d hits for %s", nrow(x$regions), nrow(x$hits), sample.id)
      fill.to.database(x$regions, x$hits)
      databaseSummary()
      } # for file

} # fill.all.samples.by.chromosome
#------------------------------------------------------------------------------------------------------------------------
# tss of chr21 gene app produces only one motif hit in hint, but 1188 in piq
# try to figure that out
# explore piq motif hit chr21:25879610-25879624 for MA0518.1
explore.paucity.of.results <- function()
{
  tbl <- readHintTable(hint.path, "ENCSR000DBY.hint.bed", nrows=NA, chromosome="chr21")

  APP.tss <- 25880550
  loc.chrom <- "chr21"
  loc.start <- APP.tss - 1000
  loc.stop  <- APP.tss + 1000
                    
  tbl.app <- subset(tbl, footprint.start >= loc.start & footprint.end <= loc.stop)

} # explore.paucity.of.results
#------------------------------------------------------------------------------------------------------------------------
if(!interactive()){
   chromosomes <- paste("chr", 1:22, sep="")
   for(chromosome in chromosomes)
      fill.all.samples.by.chromosome(chromosome)
   }
