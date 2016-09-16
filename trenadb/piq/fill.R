library(RPostgreSQL)
library(GenomicRanges)
library(RUnit)
source("../regionAndHitsSchemas.R")
#------------------------------------------------------------------------------------------------------------------------
piq.path <- "/local/Cory/for_Paul/piq_complete"
#------------------------------------------------------------------------------------------------------------------------
#  createdb -U pshannon piq
#------------------------------------------------------------------------------------------------------------------------
if(!exists("db.piq"))
   db.piq <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="piqTest", host="whovian")

knownLocs <- new.env(parent=emptyenv())
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test.createMotifNameMap()
  test.processFile()
  test.splitTableIntoRegionsAndHits()
  #test.fill.to.database()
   
} # runTests
#------------------------------------------------------------------------------------------------------------------------
# create a list to translate piq's prefereed motif name (eg, "MA00061AhrArnt") to what we use ("MA0006.1")
createMotifNameMap <- function()
{
   x <- scan("all_motifs.meme", what=character(0), sep="\n")
   motif.lines <- x[grepl("^MOTIF", x)]
   motif.lines <- sub("MOTIF ", "", motif.lines)
   motif.name.tokens <- strsplit(motif.lines, " ")
   motif.names <- as.list(unlist(lapply(motif.name.tokens, "[", 1)))
   compressed.names.0 <- gsub(".", "", fixed=TRUE, motif.lines)  
   compressed.names.1 <- gsub(" ", "", fixed=TRUE, compressed.names.0)
   compressed.names.2 <- gsub(":", "", compressed.names.1)
   compressed.names.3 <- gsub("(", "", compressed.names.2, fixed=TRUE)
   compressed.names.4 <- gsub(")", "", compressed.names.3, fixed=TRUE)
   compressed.names   <- gsub("-", "", compressed.names.4, fixed=TRUE)
   names(motif.names) <- compressed.names

   motif.names

} # createMotifNameMap
#------------------------------------------------------------------------------------------------------------------------
test.createMotifNameMap <- function()
{
   printf("--- test.createMotifNameMap")
   checkEquals(map.motifName[["MA00061AhrArnt"]], "MA0006.1")

   filenames <- list.files(piq.path, pattern="*.bed")
   all.piq.motif.tf.names <- unique(unlist(lapply(strsplit(filenames, ".", fixed=TRUE), "[", 1)))
   motif.names <- unlist(lapply(all.piq.motif.tf.names, function(name) map.motifName[[name]]))
      # simple-minded check: every motif is of the form "MA0914.1" and extracted from MA09141ISL2
      # thus there should be a dot in every motif.name
   checkEquals(length(grep(".", motif.names, fixed=TRUE)), length(motif.names))

} # test.createMotifNameMap
#------------------------------------------------------------------------------------------------------------------------
# nrow only for testing
processFile <- function(fullpath, sampleName, chromosome=NA, nrows=NA)
{
  if(is.na(nrows))
     nrows <- -1

   tbl <- read.table(fullpath, sep="\t", header=FALSE, as.is=TRUE, nrows=nrows)
   colnames(tbl) <- c("chrom", "motifStart", "motifEnd", "motifName", "strand", "score1", "score2", "score3", "score4")
   if(!is.na(chromosome))
       tbl <- subset(tbl, chrom==chromosome)
   
   tbl$score5 <- rep(NA, nrow(tbl))
   tbl$score6 <- rep(NA, nrow(tbl))
   tbl$motifLength <- with(tbl,  1 + (motifEnd - motifStart))
   tbl$footprintStart <- rep(NA, nrow(tbl))
   tbl$footprintEnd <- rep(NA, nrow(tbl))
   tbl$footprintLength <- rep(NA, nrow(tbl))
   motifName <- unlist(map.motifName[tbl$motifName], use.names=FALSE)
   if(is.null(motifName))
     return(NA)

   tbl$motifName <- unlist(map.motifName[tbl$motifName], use.names=FALSE)
   tbl$method <- rep("piq", nrow(tbl))
   tbl$provenance <- rep("minid.XXXX", nrow(tbl))
   tmp <- tbl$score1 + tbl$score2 + tbl$score3 + tbl$score4
   tmp <- tmp - min(tmp)
   tbl$normalizedScore <- tmp/max(tmp)
   
   tbl$sampleID <- rep(sampleName, nrow(tbl))
   
   new.column.order <- c("chrom", "motifStart", "motifEnd", "motifName", "normalizedScore", "strand", 
                         "footprintStart", "footprintEnd",
                         "motifLength", "footprintLength",
                         "sampleID", "method", "provenance",
                         "score1", "score2", "score3", "score4", "score5", "score6")

   tbl <- tbl[, new.column.order]
   tbl

} # processFile
#------------------------------------------------------------------------------------------------------------------------
test.processFile <- function()
{
   printf("--- test.processFile")
   file <- "MA00022RUNX1.ENCSR000DBY.bed"
   full.path <- file.path(piq.path, file)
   tbl <- processFile(full.path, "ENCSR000DBY", chromosome=NA, nrows=2)

   checkEquals(colnames(tbl), c("chrom", "motifStart", "motifEnd", "motifName", "normalizedScore", "strand", 
                                "footprintStart", "footprintEnd", "motifLength", "footprintLength", "sampleID",
                                "method", "provenance", "score1", "score2", "score3", "score4", "score5", 
                                "score6"))

    checkEquals(unique(tbl$chrom), "chr1")
    checkEquals(unique(tbl$strand), "-")
    checkEquals(tbl$motifStart, c(15050, 44200))
    checkEquals(tbl$motifEnd, c(15061, 44211))
    checkTrue(all(is.na(tbl$footprintStart)))
    checkTrue(all(is.na(tbl$footprintEnd)))

} # test_processFile
#------------------------------------------------------------------------------------------------------------------------
splitTableIntoRegionsAndHits <- function(tbl, minid)
{
   tbl$loc <- with(tbl, sprintf("%s:%d-%d", chrom, motifStart, motifEnd))
    
   tbl.regions <- unique(tbl[, c("loc", "chrom", "motifStart", "motifEnd")])
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

   tbl.hits <- tbl[, c("loc", "motifName", "strand", "sampleID", "method", "score1", "score2", "score3",
                       "score4")]

      # add piq specific columns
   #browser()
   tbl.hits$length <- with(tbl, 1 + motifEnd - motifStart)
   tbl.hits$method <- "piq"
   tbl.hits$provenance <- minid
   tbl.hits$score5 <- NA
   tbl.hits$score6 <- NA
   tbl.hits$type <- "motif.in.footprint"

      # change some names

   colnames(tbl.hits)[grep("motifName", colnames(tbl.hits))] <- "name"
   colnames(tbl.hits)[grep("sampleID", colnames(tbl.hits))] <- "sample_id"

   stopifnot(all(colnames(tbl.hits) %in% hit.schema()))

      # now rearrange the columns to the standard form, as expected by the database
   tbl.hits <- tbl.hits[, hit.schema()]
   invisible(list(regions=tbl.regions, hits=tbl.hits))
    
} # splitTableIntoRegionsAndHits    
#------------------------------------------------------------------------------------------------------------------------
test.splitTableIntoRegionsAndHits <- function()
{
   printf("--- test.splitTableIntoRegionsAndHits")
    
   knownLocs <<- new.env(parent=emptyenv())
   file <- "MA00022RUNX1.ENCSR000DBY.bed"
   full.path <- file.path(piq.path, file)
   tbl <- processFile(full.path, "ENCSR000DBY", chromosome=NA, nrows=10)
   minid <- "piq.minid"
   x <- splitTableIntoRegionsAndHits(tbl, minid)
   checkEquals(dim(x$regions), c(10, 4))
   checkEquals(dim(x$hits), c(10, 14))
   checkEquals(length(names(knownLocs)), 10)

} # test.splitTableIntoRegionsAndHits
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

   file <- "MA00022RUNX1.ENCSR000DBY.bed"
   full.path <- file.path(piq.path, file)
   tbl <- processFile(full.path, "ENCSR000DBY", nrows=10)
   minid <- "piq.minid"
   x <- splitTableIntoRegionsAndHits(tbl, minid)
   fill.to.database(x$regions, x$hits)
   checkEquals(sort(dbListTables(db.piq)), c("hits", "regions"))

   tbl.regions <- dbGetQuery(db.piq, "select * from regions")
   checkEquals(dim(tbl.regions), c(10, 4))
   checkTrue(all(region.schema() == colnames(tbl.regions)))

   tbl.hits <- dbGetQuery(db.piq, "select * from hits")
   checkTrue(all(hit.schema() == colnames(tbl.hits)))
   checkEquals(dim(tbl.hits), c(10, 13))
   checkEquals(unique(tbl.hits$type), "motif.in.footprint")
   checkTrue(all(tbl.hits$loc %in% tbl.regions$loc))

} # test.fill.to.database
#------------------------------------------------------------------------------------------------------------------------
databaseSummary <- function()
{
    region.count <- dbGetQuery(db.piq, "select count(*) from regions")[1,1]
    hit.count <- dbGetQuery(db.piq, "select count(*) from hits")[1,1]
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

   all.files <- list.files(piq.path, pattern="*.bed$")
   minid <- "piq.minid.tbd"

   for(file in all.files){
      printf("---- %s (%s) (%d/%d)", file, chromosome, grep(file, all.files), length(all.files))
      sampleID <- strsplit(file, ".", fixed=TRUE)[[1]][2]
      full.path <- file.path(piq.path, file)
      tbl <- processFile(full.path, sampleID, chromosome) #, nrows=10)
      x <- splitTableIntoRegionsAndHits(tbl, minid)
      fill.to.database(x$regions, x$hits)
      databaseSummary()
      } # for file

} # fill.all.samples.by.chromosome
#------------------------------------------------------------------------------------------------------------------------
map.motifName <- createMotifNameMap()

if(!interactive()){
   fill.all.samples.by.chromosome("chr19")
   }
