library(RPostgreSQL)
library(GenomicRanges)
library(RUnit)
source("../regionAndHitsSchemas.R")
db.fimo <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="fimo", host="whovian")
#------------------------------------------------------------------------------------------------------------------------
chipseq.path <- "/proj/price1/sament/lymphoblast_trn/known_tfbs/hg38"
#------------------------------------------------------------------------------------------------------------------------
#  createdb -U pshannon chipseq
#------------------------------------------------------------------------------------------------------------------------
if(!exists("db"))
   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="chipseq", host="whovian")

knownLocs <- new.env(parent=emptyenv())
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test.createMotifNameMap()
  test.processFile()
  test.splitTableIntoRegionsAndHits()
   
} # runTests
#------------------------------------------------------------------------------------------------------------------------
# lots of multiple assignments, from 1 tf to multiple motifs:
#  tbl[grep("^ATF2$", tbl$tfs),]$motif
#    "MA0089.1" "MA0150.2" "MA0490.1" "MA0491.1" "MA0492.1" "MA0591.1" "MA0834.1" "MA0840.1"
# won't use this now.  instead, the name of the region will be the gene (tf) name
createMotifNameMap <- function()
{
  tbl <- read.table("../motifGeneMap/motif_to_tf_mappings_with_tfclass_include_multiple.csv", sep=",",
                    header=TRUE, stringsAsFactors=FALSE)


} # createMotifNameMap
#------------------------------------------------------------------------------------------------------------------------
test.createMotifNameMap <- function()
{
   printf("--- test.createMotifNameMap")
   tbl.motifGene <<- createMotifNameMap()
   checkEquals(dim(tbl.motifGene), c(9289, 2))
   checkEquals(tbl.motifGene[grep("^ATF2$", tbl.motifGene$tfs),]$motif,
               c("MA0089.1", "MA0150.2", "MA0490.1", "MA0491.1", "MA0492.1", "MA0591.1", "MA0834.1", "MA0840.1"))

} # test.createMotifNameMap
#------------------------------------------------------------------------------------------------------------------------
# nrow only for testing
processFile <- function(fullpath, tfName, chromosome=NA, nrows=NA)
{
   tbl <- read.table(fullpath, sep="\t", header=FALSE, as.is=TRUE)
   colnames(tbl) <- c("chrom", "motifStart", "motifEnd", "genes", "score1", "category")
   tbl$motifName <- tfName

   if(!is.na(chromosome))
       tbl <- subset(tbl, chrom==chromosome)

   if(nrow(tbl) == 0)
       return (data.frame())

   if(!is.na(nrows))
     tbl <- tbl[1:nrows,]
   
   tbl$score2 <- rep(NA, nrow(tbl))
   tbl$score3 <- rep(NA, nrow(tbl))
   tbl$score4 <- rep(NA, nrow(tbl))
   tbl$score5 <- rep(NA, nrow(tbl))
   tbl$score6 <- rep(NA, nrow(tbl))
   tbl$strand <- rep(NA, nrow(tbl))

   tbl$motifLength <- with(tbl,  1 + (motifEnd - motifStart))
   tbl$footprintStart <- rep(NA, nrow(tbl))
   tbl$footprintEnd <- rep(NA, nrow(tbl))
   tbl$footprintLength <- rep(NA, nrow(tbl))

   tbl$method <- rep("ChIP-SEQ", nrow(tbl))
   tbl$provenance <- rep("ChIP-minid.tbd", nrow(tbl))
   tmp <- tbl$score1
   tmp <- tmp - min(tmp)
   tbl$normalizedScore <- tmp/max(tmp)
   
   tbl$sampleID <- rep("pooled", nrow(tbl))
   
   new.column.order <- c("chrom", "motifStart", "motifEnd", "motifName", "normalizedScore", "strand", 
                         "footprintStart", "footprintEnd",
                         "motifLength", "footprintLength",
                         "sampleID", "method", "provenance",
                         "score1", "score2", "score3", "score4", "score5", "score6")
   missing <- setdiff(new.column.order, colnames(tbl))
   tbl[, new.column.order]

} # processFile
#------------------------------------------------------------------------------------------------------------------------
test.processFile <- function()
{
   printf("--- test.processFile")
   file <- list.files(chipseq.path)[1]
   tf.name <- strsplit(file, "_")[[1]][1]
   full.path <- file.path(chipseq.path, file)
   
   tbl <- processFile(full.path, tf.name, chromosome="chr21", nrows=2)

   checkEquals(colnames(tbl), c("chrom", "motifStart", "motifEnd", "motifName", "normalizedScore", "strand", 
                                "footprintStart", "footprintEnd", "motifLength", "footprintLength", "sampleID",
                                "method", "provenance", "score1", "score2", "score3", "score4", "score5", 
                                "score6"))

    checkEquals(unique(tbl$chrom), "chr21")
    checkTrue(all(tbl$sampleID == "pooled"))
    checkTrue(all(is.na(tbl$strand)))
    checkEquals(tbl$motifStart, c(14086719, 14096239))
    checkEquals(tbl$motifEnd,   c(14086869, 14096389))
    checkTrue(all(is.na(tbl$footprintStart)))
    checkTrue(all(is.na(tbl$footprintEnd)))

      # now test a file with now chr21 peaks
   file <- list.files(chipseq.path)[21]
   tf.name <- strsplit(file, "_")[[1]][1]
   full.path <- file.path(chipseq.path, file)
   chromosomes.read <- unique(read.table(full.path, sep="\t", header=FALSE, as.is=TRUE, stringsAsFactors=FALSE))
   checkTrue(! "chr21" %in% chromosomes.read)
   
   tbl <- processFile(full.path, tf.name, chromosome="chr21", nrows=2)
   checkEquals(dim(tbl), c(0,0))

} # test_processFile
#------------------------------------------------------------------------------------------------------------------------
splitTableIntoRegionsAndHits <- function(tbl, minid)
{
   tbl$loc <- with(tbl, sprintf("%s:%d-%d", chrom, motifStart, motifEnd))
   tbl$length <- with(tbl, 1 + motifEnd-motifStart)
   
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
   selected.columns <- c("loc", "motifName", "strand", "sampleID", "method", "score1", "score2", "score3", "score4")
   stopifnot(all(selected.columns %in% colnames(tbl)))
   tbl.hits <- tbl[, ]

      # add piq specific columns
   tbl.hits$method <- "cusanovitch"
   tbl.hits$provenance <- minid
   tbl.hits$score5 <- NA
   tbl.hits$score6 <- NA
   tbl.hits$type <- "chipseq.peak"

      # change some names

   colnames(tbl.hits)[grep("motifName", colnames(tbl.hits))] <- "name"
   colnames(tbl.hits)[grep("sampleID", colnames(tbl.hits))] <- "sample_id"

   #browser()
   stopifnot(all(hit.schema() %in% colnames(tbl.hits)))

      # now rearrange the columns to the standard form, as expected by the database
   tbl.hits <- tbl.hits[, hit.schema()]
   invisible(list(regions=tbl.regions, hits=tbl.hits))
    
} # splitTableIntoRegionsAndHits    
#------------------------------------------------------------------------------------------------------------------------
test.splitTableIntoRegionsAndHits <- function()
{
   printf("--- test.splitTableIntoRegionsAndHits")
    
   knownLocs <<- new.env(parent=emptyenv())
   file <- "ATF2_lymphoblast_binding_sites.bed.hg38.bed"
   full.path <- file.path(chipseq.path, file)
   tbl <- processFile(full.path, "ATF2", chromosome=NA, nrows=10)
   minid <- "chipseq.minid.tbd"
   x <- splitTableIntoRegionsAndHits(tbl, minid)

   checkEquals(sort(names(x)), c("hits", "regions"))
   checkEquals(dim(x$regions), c(10, 4))

   checkEquals(dim(x$hits), c(10, 14))
   checkTrue(all(x$hits$name == "ATF2"))
   checkTrue(all(x$hits$type == "chipseq.peak"))
   checkEqualsNumeric(mean(x$hits$score1), 134.2, tol=0.1)
   checkEquals(length(names(knownLocs)), 10)

   checkTrue(all(x$hits$loc %in% x$regions$loc))

} # test.splitTableIntoRegionsAndHits
#------------------------------------------------------------------------------------------------------------------------
fill.to.database <- function(tbl.regions, tbl.hits)
{
   appendToRegionsTable(tbl.regions)
   appendToHitsTable(tbl.hits)
    
} # fill.to.database
#------------------------------------------------------------------------------------------------------------------------
disabled.test.fill.to.database <- function()
{
   printf("--- test.fill.to.database")
   createEmptyDatabaseTables()
   knownLocs <<- new.env(parent=emptyenv())

   file <- "ATF2_lymphoblast_binding_sites.bed.hg38.bed"
   full.path <- file.path(chipseq.path, file)
   tbl <- processFile(full.path, "ATF2", nrows=10)
   minid <- "chipseq.minid.tbd"
   x <- splitTableIntoRegionsAndHits(tbl, minid)
   fill.to.database(x$regions, x$hits)
   checkEquals(sort(dbListTables(db)), c("hits", "regions"))

   tbl.regions <- dbGetQuery(db, "select * from regions")
   checkEquals(dim(tbl.regions), c(10, 4))
   checkTrue(all(region.schema() == colnames(tbl.regions)))

   tbl.hits <- dbGetQuery(db, "select * from hits")
   checkTrue(all(hit.schema() == colnames(tbl.hits)))
   checkEquals(dim(tbl.hits), c(10, 13))
   checkEquals(unique(tbl.hits$type), "chipseq.peak")
   checkTrue(all(tbl.hits$loc %in% tbl.regions$loc))

} # disabled.test.fill.to.database
#------------------------------------------------------------------------------------------------------------------------
databaseSummary <- function()
{
    region.count <- dbGetQuery(db, "select count(*) from regions")[1,1]
    hit.count <- dbGetQuery(db, "select count(*) from hits")[1,1]
    printf("%d hits in %d regions", hit.count, region.count)

} # databaseSummary
#------------------------------------------------------------------------------------------------------------------------
#createEmptyDatabaseTables <- function()
#{
#   system("/bin/psql -f createTables.sql")
#    
#} # appendToRegionsTable
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

   all.files <- list.files(chipseq.path, pattern="*.bed$")
   minid <- "chipseq.minid.tbd"

   for(file in all.files){
      printf("---- %s (%s) (%d/%d)", file, chromosome, grep(file, all.files), length(all.files))
      tf.name <- strsplit(file, "_")[[1]][1]
      full.path <- file.path(chipseq.path, file)
      tbl <- processFile(full.path, tf.name, chromosome) #, nrows=10)
      if(nrow(tbl) == 0)
          next;
      x <- splitTableIntoRegionsAndHits(tbl, minid)
      fill.to.database(x$regions, x$hits)
      databaseSummary()
      } # for file

} # fill.all.samples.by.chromosome
#------------------------------------------------------------------------------------------------------------------------
map.motifName <- createMotifNameMap()

if(!interactive()){
   chromosomes <- paste("chr", 1:22, sep="")
   for(chromosome in chromosomes)
      fill.all.samples.by.chromosome(chromosome)
   }
