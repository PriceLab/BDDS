library(RPostgreSQL)
library(GenomicRanges)
library(RUnit)
source("../../src/regionAndHitsSchemas.R")
#------------------------------------------------------------------------------------------------------------------------
wellington.path <- "/local/Ben/BDDS/footprints/functionalTests/output/wellington"
test.sampleID <- "ENCSR000DBY"

#------------------------------------------------------------------------------------------------------------------------
if(!exists("db.wellington"))
   db.wellington <- dbConnect(PostgreSQL(), user= "trenatest", password="trenatest", dbname="testwellington", host="whovian")

if(!exists("db.fimo"))
   db.fimo <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="fimo", host="whovian")

knownLocs <- new.env(parent=emptyenv())

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test.readWellingtonTable()
  test.mergeFootprintsWithFimo()
  test.splitTableIntoRegionsAndWellingtonHits()
  #test.fill.to.database()
    
  #test.combineFootprintsAndFimo()
  #test.combineFootprintsAndDatabasedFimo()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
readWellingtonTable <- function(directory, sampleID, nrows=NA, chromosome=NA)
{
   filename <- grep(sampleID, list.files(directory), v=TRUE)
   full.path <- file.path(directory, filename)

   if(!file.exists(full.path))
       return(data.frame)
   
   tbl <- read.table(full.path, sep="\t", as.is=TRUE)
   colnames(tbl) <- c("chrom", "start", "end", "name", "score", "strand")
   tbl$chrom <- paste("chr", tbl$chrom, sep="")
   if(!is.na(chromosome))
      tbl <- subset(tbl, chrom==chromosome)

   if(!is.na(nrows))
      tbl <- tbl[1:nrows,]

   invisible(tbl)

} # readWellingtonTable
#------------------------------------------------------------------------------------------------------------------------
test.readWellingtonTable <- function()
{
   printf("--- test.readWellingtonTable")

   tbl <- readWellingtonTable(wellington.path, test.sampleID, 5, "chr21")
   checkEquals(dim(tbl), c(5,6))
   checkEquals(colnames(tbl), c("chrom", "start", "end", "name", "score", "strand"))
   checkEquals(unique(tbl$chrom), "chr21")

     # now read without chrom or nrow constraints
   tbl <- readWellingtonTable(wellington.path, test.sampleID)
   checkEquals(ncol(tbl), 6)
   checkTrue(nrow(tbl) > 30000)
   checkEquals(colnames(tbl), c("chrom", "start", "end", "name", "score", "strand"))
   checkEquals(head(sort(unique(tbl$chrom))), c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14"))
   
} # test.readWellingtonTable
#------------------------------------------------------------------------------------------------------------------------
mergeFimoWithFootprints <- function(tbl.fp, sampleID)
{
  chromosome <- unique(tbl.fp$chrom)
    # enforce treatment of just one chromosome at a time
  stopifnot(length(chromosome) == 1)
  min.pos <- min(tbl.fp$start)
  max.pos <- max(tbl.fp$end)
  
  fimo.chromosome <- sub("chr", "", chromosome)
  query <- sprintf("select * from hg38 where chr='%s' and start >= %d and endpos <= %d",
                   fimo.chromosome, min.pos, max.pos)

  tbl.fimo <- dbGetQuery(db.fimo, query)
  colnames(tbl.fimo) <- c("motif", "chrom", "motif.start", "motif.end", "motif.strand", "fimo.score",
                          "fimo.pvalue", "empty", "motif.sequence")
  tbl.fimo <- tbl.fimo[, -grep("empty", colnames(tbl.fimo))]
  tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")

  gr.fimo <- with(tbl.fimo, GRanges(seqnames=chrom, IRanges(start=motif.start, end=motif.end)))

    # --- get some footprints

  gr.wellington <- with(tbl.fp,   GRanges(seqnames=chrom, IRanges(start=start, end=end)))
  tbl.overlaps <- as.data.frame(findOverlaps(gr.fimo, gr.wellington, type="within"))

  tbl.fimo$loc <- with(tbl.fimo, sprintf("%s:%d-%d", chrom, motif.start, motif.end))
  tbl.fimo$method <- "WELLINGTON"
  tbl.fimo$sample_id <- sampleID
  tbl.regions <- tbl.fimo[tbl.overlaps$queryHits,]

  tbl.regions <- cbind(tbl.regions, wellington.score=tbl.fp[tbl.overlaps$subjectHits, "score"])
  invisible(tbl.regions)

} # mergeFimoWithFootprints
#------------------------------------------------------------------------------------------------------------------------
test.mergeFootprintsWithFimo <- function()
{
   printf("--- test.mergeFootprintsWithFimo")
   tbl.fp <- readWellingtonTable(wellington.path, test.sampleID, nrow=3, "chr21")
   tbl <- mergeFimoWithFootprints(tbl.fp, test.sampleID)
   checkEquals(ncol(tbl), 12)
   checkEquals(sort(colnames(tbl)),
               c("chrom", "fimo.pvalue", "fimo.score", "loc", "method",
                 "motif", "motif.end", "motif.sequence", "motif.start", "motif.strand", "sample_id",
                 "wellington.score"))
   checkTrue(nrow(tbl) >= 8)
   duplicated.loc <- "chr21:5290098-5290107"
   checkEquals(length(grep(duplicated.loc, tbl$loc)), 6)
     #  3 distinct motifs mapped to this region
   checkEquals(sort(subset(tbl, loc == duplicated.loc)$motif),
               c("MA0156.2", "MA0474.2", "MA0475.2", "MA0624.1", "MA0625.1", "MA0764.1"))
   invisible(tbl)

} # test.mergeFootprintsWithFimo
#------------------------------------------------------------------------------------------------------------------------
splitTableIntoRegionsAndWellingtonHits <- function(tbl, minid)
{
   tbl.regions <- unique(tbl[, c("loc", "chrom", "motif.start", "motif.end")])
   colnames(tbl.regions) <- region.schema() # 29
    # c("loc", "chrom", "motif_start", "motif_end")

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

   tbl.hits <- tbl[, c("loc", "motif", "motif.strand", "sample_id", "method", "wellington.score",
                       "fimo.score", "fimo.pvalue")]
   tbl.hits$length <- with(tbl, 1 + motif.end - motif.start)
   tbl.hits$provenance <- minid
   tbl.hits$score4 <- NA
   tbl.hits$score5 <- NA
   tbl.hits$score6 <- NA
   tbl.hits$type <- "motif.in.footprint"
   coi <- c("loc", "type", "motif", "length", "motif.strand", "sample_id", "method", "provenance", 
            "wellington.score", "fimo.score", "fimo.pvalue", "score4", "score5", "score6")
   tbl.hits <- tbl.hits[, coi]
   colnames(tbl.hits) <- hit.schema()
   invisible(list(regions=tbl.regions, hits=tbl.hits))
    
} # splitTableIntoRegionsAndWellingtonHits    
#------------------------------------------------------------------------------------------------------------------------
test.splitTableIntoRegionsAndWellingtonHits <- function(tbl)
{
   printf("--- test.splitTableIntoRegionsAndWellingtonHits")
   tbl.fp <- readWellingtonTable(wellington.path, test.sampleID, nrow=3, "chr21")
   tbl <- mergeFimoWithFootprints(tbl.fp, test.sampleID)

   x <- splitTableIntoRegionsAndWellingtonHits(tbl, "minid")
   checkEquals(sort(names(x)), c("hits", "regions"))
   checkEquals(colnames(x$regions), region.schema()) #
   checkEquals(colnames(x$hits), hit.schema())

} # test.splitTableIntoRegionsAndWellingtonHits    
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

   tbl.fp <- readWellingtonTable(wellington.path, test.sampleID, nrow=3, "chr21")
   tbl <- mergeFimoWithFootprints(tbl.fp, test.sampleID)
   x <- splitTableIntoRegionsAndWellingtonHits(tbl, "minid")

   fill.to.database(x$regions, x$hits)
   checkEquals(sort(dbListTables(db.wellington)), c("hits", "regions"))

   tbl.regions <- dbGetQuery(db.wellington, "select * from regions")
   checkEquals(dim(tbl.regions), c(2, 4))
   checkTrue(all(region.schema() == colnames(tbl.regions)))

   tbl.hits <- dbGetQuery(db.wellington, "select * from hits")
   checkTrue(all(hit.schema() == colnames(tbl.hits)))
   checkEquals(dim(tbl.hits), c(8, 14))
   checkEquals(unique(tbl.hits$type), "motif.in.footprint")
   checkTrue(all(tbl.hits$loc %in% tbl.regions$loc))

} # test.fill.to.database
#------------------------------------------------------------------------------------------------------------------------
databaseSummary <- function()
{
    region.count <- dbGetQuery(db.wellington, "select count(*) from regions")[1,1]
    hit.count <- dbGetQuery(db.wellington, "select count(*) from hits")[1,1]
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
   knownLocs <<- new.env(parent=emptyenv())

   all.sampleIDs <- unlist(lapply(strsplit(list.files(wellington.path, "ENCSR.*.bed"), ".", fixed=TRUE), "[", 1))
   minid <- "corys.wellington.minid"

   for(sampleID in all.sampleIDs){
      printf("---- %s (%s) (%d/%d)", sampleID, chromosome, grep(sampleID, all.sampleIDs), length(all.sampleIDs))
      tbl.wellington <- readWellingtonTable(wellington.path, sampleID, NA, chromosome)
      tbl <- mergeFimoWithFootprints(tbl.wellington, sampleID)
      x <- splitTableIntoRegionsAndWellingtonHits(tbl, minid)
      printf("filling %d regions, %d hits for %s", nrow(x$regions), nrow(x$hits), sampleID)
      fill.to.database(x$regions, x$hits)
      databaseSummary()
      } # for file

} # fill.all.samples.by.chromosome
#------------------------------------------------------------------------------------------------------------------------
examine.region <- function(chromosome, start, end, sampleIDs)
{
   knownLocs <<- new.env(parent=emptyenv())

   all.sampleIDs <- unlist(lapply(strsplit(list.files(wellington.path, "*.bed"), ".", fixed=TRUE), "[", 1))
   sampleIDs <- intersect(sampleIDs, all.sampleIDs)

   minid <- "corys.wellington.minid"

   for(sampleID in sampleIDs){
      printf("---- %s (%s) (%d/%d)", sampleID, chromosome, grep(sampleID, all.sampleIDs), length(all.sampleIDs))
      tbl.wellington <- readWellingtonTable(wellington.path, sampleID, NA, chromosome)
      tbl <- mergeFimoWithFootprints(tbl.wellington, sampleID)
      x <- splitTableIntoRegionsAndWellingtonHits(tbl, minid)
      printf("filling %d regions, %d hits for %s", nrow(x$regions), nrow(x$hits), sampleID)
      #fill.to.database(x$regions, x$hits)
      #databaseSummary()
      } # for file

} # examine.region
#------------------------------------------------------------------------------------------------------------------------
test_examine.region <- function()
{
   printf("--- test_examine.region")
   sampleIDs <- c("ENCSR000EJB")
   chrom <- "chr19"
   start <- 44903773
   end <- 44903785
   examine.region(chrom, start, end, sampleIDs)

} # test_examine.region
#------------------------------------------------------------------------------------------------------------------------
if(!interactive()){
    chromosomes <- paste("chr", c(1:18, 20:22), sep="")
    #chromosomes <- paste("chr", c(19), sep="")
    for(chromosome in chromosomes)
        fill.all.samples.by.chromosome(chromosome)
    }

