source("../src/dependencies.R")
source("../src/dbFunctions.R")
source("../src/tableParsing.R")
source("../src/tests.R")

#-------------------------------------------------------------------------------
# path to output of the makefile based tests for wellington on CHR 19
wellington.path <- "/local/Ben/BDDS/footprints/functionalTests/output/wellington"
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
fill.all.samples.by.chromosome <- function(chromosome = "chr19", minid = "temp.filler.minid")
{
   knownLocs <<- new.env(parent=emptyenv())

   all.sampleIDs <- unlist(lapply(strsplit(list.files(wellington.path, "ENCSR.*.bed"), ".", fixed=TRUE), "[", 1))

   for(sampleID in all.sampleIDs){
      printf("---- %s (%s) (%d/%d)", sampleID, chromosome, grep(sampleID, all.sampleIDs), length(all.sampleIDs))
      tbl.wellington <- readWellingtonTable(wellington.path, sampleID, NA, chromosome)
      print("Wellington table read. Merging with Fimo...")
      tbl <- mergeFimoWithFootprints(tbl.wellington, sampleID)
      print("Merged. Now splitting table to regions and hits...")
      x <- splitTableIntoRegionsAndWellingtonHits(tbl, minid)
      printf("filling %d regions, %d hits for %s", nrow(x$regions), nrow(x$hits), sampleID)
      fill.to.database(x$regions, x$hits, db.wellington)
      databaseSummary()
      } # for file

} # fill.all.samples.by.chromosome
#-------------------------------------------------------------------------------
if(!interactive()){
    #chromosomes <- paste("chr", c(1:18, 20:22), sep="")
    chromosomes <- paste("chr", c(19), sep="")
    for(chromosome in chromosomes)
        fill.all.samples.by.chromosome(chromosome)
    }