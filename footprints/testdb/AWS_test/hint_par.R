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
source("../src/main_SJ.R")

#-------------------------------------------------------------------------------
# set path to hint output 
printf = function(...) print(noquote(sprintf(...)))
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
    # pass a string; don't connect yet
    db.wellington <- "testwellington_localhost"

if(!exists("db.fimo"))
#   db.fimo <- getDBConnection("fimo_localhost")
    # pass a string; don't connect yet
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
    chromosomes <- paste("chr", c(13:14), sep="")
    for(chromosome in chromosomes)
        fillAllSamplesByChromosome(chromosome = chromosome,
                                   dbConnection = db.wellington,
                                   fimo = db.fimo,
                                   minid = "testhint_par.minid",
                                   #dbUser = "ben",
                                   dbUser = "trena",
                                   dbTable = "testwellington",
                                   sourcePath = data.path,
                                   isTest = FALSE,
                                   method = "HINT")
#    }
