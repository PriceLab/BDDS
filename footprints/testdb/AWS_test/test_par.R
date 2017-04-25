# users can either step through this file, or call this file with 
# r -f example.R

# THIS ASSUMES THAT THE TESTHINT DATABASE EXISTS. The recipe for building that
# database is in ../dbInitialization/createHintTest.sql

# THIS EXAMPLE USES THE HINT OUTPUT MADE BY RUNNING make hint at /BDDS/footprints/functionalTests/ 

#-------------------------------------------------------------------------------
# load functions, libraries and dependencies
source("../src/dependencies.R")
source("../src/dbFunctions.R")
source("../src/tableParsing.R")
source("../src/tests.R")
source("../src/main.R")

#-------------------------------------------------------------------------------
# set path to hint output 

# for makefile based tests for hint on CHR 19 on whovian, use this:
data.path <- "/scratch/data/test_set"
#test.sampleID <- "ENCSR000DBY"

# on globus genomics machines, use this:
# hint.path <- "/local/Ben/BDDS/footprints/functionalTests/output/hint"
# test.sampleID <- "ENCSR000DBY"

#-------------------------------------------------------------------------------
# establish database connections:
# for whovian, use this:

if(!exists("db.wellington"))
#   db.wellington <- getDBConnection("testwellington_localhost")
    db.wellington <- "testwellington_localhost"

if(!exists("db.fimo"))
#   db.fimo <- getDBConnection("fimo_localhost")
    db.fimo <- "fimo_localhost"
# for bdds-rds-globusgenomics.org, use:
# if(!exists("db.hint"))
#   db.hint <- getDBConnection("testhint")
# 
# if(!exists("db.fimo"))
#   db.fimo <- getDBConnection("fimo")
#-------------------------------------------------------------------------------
#if(!interactive()){
#    chromosomes <- paste("chr", c(1:22), sep="")
    chromosomes <- paste("chr", c(21,22), sep="")
#    for(chromosome in chromosomes)

  # Create parallel structure here
  library(foreach); library(doParallel)
  cores <- detectCores()
  cl <- makeCluster(cores[1] - 1)
  registerDoParallel(cl)
  
  #give it the readDataTable function and other variables
  clusterExport(cl, varlist = c("readDataTable","data.path",
				"db.fimo", "db.wellington",
				"mergeFimoWithFootprints",
				"fillToDatabase", "databaseSummary",
				"splitTableIntoRegionsAndHits",
				"getDBConnection","region.schema",
				"printf","hit.schema",
				"appendToRegionsTable","appendToHitsTable"),
				envir = environment())
  junk <- clusterEvalQ(cl, library(GenomicRanges))
  junk <- clusterEvalQ(cl, library(RPostgreSQL))
#  junk <- clusterEvalQ(cl, source("./fillFunctions.R"))

  foreach(i=1:length(chromosomes)) %dopar% {
    fillAllSamplesByChromosome(chromosome = chromosomes[[i]],
                                   dbConnection = db.wellington,
                                   fimo = db.fimo,
                                   minid = "testhint_par.minid",
                                   #dbUser = "ben",
                                   dbUser = "trena",
                                   dbTable = "testwellington",
                                   sourcePath = data.path,
                                   isTest = FALSE,
                                   method = "HINT")
				   }				  
#    }
