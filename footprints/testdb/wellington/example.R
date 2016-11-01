# users can either step through this file, or call this file with 
# r -f example.R

#-------------------------------------------------------------------------------
# load functions and dependencies
source("../src/dependencies.R")
source("../src/dbFunctions.R")
source("../src/tableParsing.R")
source("../src/tests.R")
source("../src/main.R")

#-------------------------------------------------------------------------------
# set path to wellington output 

# for makefile based tests for wellington on CHR 19 on whovian, use this:
wellington.path <- "/local/Ben/BDDS/footprints/functionalTests/output/wellington"
test.sampleID <- "ENCSR000DBY"

# on globus genomics machines, use this:
# wellington.path <- "/local/Ben/BDDS/footprints/functionalTests/output/wellington"
# test.sampleID <- "ENCSR000DBY"

#-------------------------------------------------------------------------------
# establish database connections:
# for whovian, use this:

# if(!exists("db.wellington"))
#    db.wellington <- getDBConnection("testwellington_whovian")
# 
# if(!exists("db.fimo"))
#    db.fimo <- getDBConnection("fimo_whovian")
 
# for bdds-rds-globusgenomics.org, use:
if(!exists("db.wellington"))
  db.wellington <- getDBConnection("testwellington")

if(!exists("db.fimo"))
  db.fimo <- getDBConnection("fimo")
#-------------------------------------------------------------------------------

if(!interactive()){
    #chromosomes <- paste("chr", c(1:18, 20:22), sep="")
    chromosomes <- paste("chr", c(19), sep="")
    for(chromosome in chromosomes)
        fillAllSamplesByChromosome(chromosome)
    }