source("../src/dependencies.R")
source("../src/DBFunctions.R")
source("../src/TableParsing.R")
source("../src/Tests.R")

#-------------------------------------------------------------------------------
# path to output of the makefile based tests for wellington on CHR 19
wellington.path <- "../../functionalTests/output/wellington"
test.sampleID <- "ENCSR000DBY"

#-------------------------------------------------------------------------------
if(!exists("db.wellington"))
   db.wellington <- getDBConnection(dbname="testwellington")

if(!exists("db.fimo"))
   db.fimo <- getDBConnection(user= "trena", 
                              password="trena", 
                              dbname="fimo", 
                              host="whovian")

knownLocs <- new.env(parent=emptyenv())
#-------------------------------------------------------------------------------
#TODO: fix fill.to.database call
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
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
if(!interactive()){
    chromosomes <- paste("chr", c(1:18, 20:22), sep="")
    #chromosomes <- paste("chr", c(19), sep="")
    for(chromosome in chromosomes)
        fill.all.samples.by.chromosome(chromosome)
    }