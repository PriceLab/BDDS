# users can either step through this file, or call this file with 
# r -f example.R

# THIS ASSUMES THAT THE TESTHINT DATABASE EXISTS. The recipe for building that
# database is in ../dbInitialization/createHintTest.sql

# THIS EXAMPLE USES THE HINT OUTPUT MADE BY RUNNING make hint at /BDDS/footprints/functionalTests/ 

#-------------------------------------------------------------------------------
# load functions and dependencies
source("../src/dependencies.R")
source("../src/dbFunctions.R")
source("../src/tableParsing.R")
source("../src/tests.R")
source("../src/main.R")

#-------------------------------------------------------------------------------
# set path to hint output 

# for makefile based tests for hint on CHR 19 on whovian, use this:
hint.path <- "/local/Ben/BDDS/footprints/functionalTests/output/hint"
#test.sampleID <- "ENCSR000DBY"

# on globus genomics machines, use this:
# hint.path <- "/local/Ben/BDDS/footprints/functionalTests/output/hint"
# test.sampleID <- "ENCSR000DBY"

#-------------------------------------------------------------------------------
# establish database connections:
# for whovian, use this:

if(!exists("db.hint"))
   db.hint <- getDBConnection("testhint_whovian")

if(!exists("db.fimo"))
   db.fimo <- getDBConnection("fimo_whovian")
 
# for bdds-rds-globusgenomics.org, use:
# if(!exists("db.hint"))
#   db.hint <- getDBConnection("testhint")
# 
# if(!exists("db.fimo"))
#   db.fimo <- getDBConnection("fimo")
#-------------------------------------------------------------------------------

if(!interactive()){
    #chromosomes <- paste("chr", c(1:18, 20:22), sep="")
    chromosomes <- paste("chr", c(19), sep="")
    for(chromosome in chromosomes)
        fillAllSamplesByChromosome(chromosome = chromosome,
                                   dbConnection = db.hint,
                                   fimo = db.fimo,
                                   minid = "testexample.filler.minid",
                                   #dbUser = "ben",
                                   dbUser = "trena",
                                   dbTable = "testhint",
                                   sourcePath = hint.path,
                                   isTest = FALSE,
                                   method = "HINT")
    }
