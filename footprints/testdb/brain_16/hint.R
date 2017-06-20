# users can either step through this file, or call this file with 
# r -f example.R

# THIS ASSUMES THAT THE TESTHINT DATABASE EXISTS. The recipe for building that
# database is in ../dbInitialization/createHintTest.sql

# THIS EXAMPLE USES THE BRAIN HINT OUTPUT MADE BY RUNNING make hint at /scratch/data/footprints

print(date())
#-------------------------------------------------------------------------------
# set path to hint output 
data.path <- "/scratch/data/footprints/seed16/brain/hint"
#-------------------------------------------------------------------------------
# establish database connections:

if(!exists("db.hint"))
    db.hint <- "brain_hint_16_localhost"

if(!exists("db.fimo"))
    db.fimo <- "fimo_localhost"
#-------------------------------------------------------------------------------
if(!interactive()){    
    chromosomes <- paste("chr", c(1:22,"X","Y","MT"), sep="")
    
    # Create parallel structure here    
    library(foreach); library(doParallel)    
    cores <- detectCores()    
    cl <- makeCluster(cores[1] - 1)    
    registerDoParallel(cl)      

    # Pass path variables and source files
    clusterExport(cl, varlist = c("data.path","db.fimo", "db.hint"),
                  envir = environment())
    
    junk <- clusterEvalQ(cl, source("../src/dependencies.R"))
    junk <- clusterEvalQ(cl, source("../src/dbFunctions.R"))
    junk <- clusterEvalQ(cl, source("../src/tableParsing.R"))
    junk <- clusterEvalQ(cl, source("../src/tests.R"))
    junk <- clusterEvalQ(cl, source("../src/main.R"))

    # Run on all 24 possible chromosomes at once
    foreach(i=1:length(chromosomes)) %dopar% {
        fillAllSamplesByChromosome(chromosome = chromosomes[[i]],
                                   dbConnection = db.hint,
                                   fimo = db.fimo,
                                   minid = "brain_hint_16.minid",
                                   dbUser = "trena",
                                   dbTable = "brain_hint_16",
                                   sourcePath = data.path,
                                   isTest = FALSE,
                                   method = "HINT")
				   }				  
}

print("Database fill complete; creating indices")

# Index the database
source("../src/dbFunctions.R")
source("../src/dependencies.R")
dbConnection <- getDBConnection(db.hint)
dbSendQuery(dbConnection, "create index regions_index on regions (loc, start, endpos);")
dbSendQuery(dbConnection, "create index hits_index on hits (loc);")
dbDisconnect(dbConnection)

print(date())
