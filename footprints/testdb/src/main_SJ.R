#-------------------------------------------------------------------------------
fillAllSamplesByChromosome <- function(dbConnection = db.wellington,
                                       fimo = db.fimo,
                                       chromosome = "chr19",
                                       minid = "temp.filler.minid",
                                       dbUser = "ben",
                                       dbTable = "testwellington",
                                       sourcePath = wellington.path,
                                       isTest = True,
                                       method = "DEFAULT")
{
  knownLocs <<- new.env(parent=emptyenv())
  
  all.sampleIDs <- unlist(lapply(strsplit(list.files(sourcePath, 
                                                     "ENCSR.*.bed$"), 
                                          ".", fixed=TRUE), "[", 1))


  # Create parallel structure here
  library(foreach); library(doParallel)
  cores <- detectCores()
  cl <- makeCluster(cores[1] - 1)
  registerDoParallel(cl)
  
  #give it the readDataTable function and other variables
  clusterExport(cl, varlist = c("readDataTable", "all.sampleIDs",
                                "sourcePath", "chromosome",
				"fimo", "dbConnection",
				"method", "minid",
				"dbUser", "dbTable",
				"mergeFimoWithFootprints",
				"fillToDatabase", "databaseSummary",
				"splitTableIntoRegionsAndHits",
				"getDBConnection","region.schema",
				"knownLocs","printf","hit.schema",
				"appendToRegionsTable","appendToHitsTable"),
				envir = environment())
  junk <- clusterEvalQ(cl, library(GenomicRanges))
  junk <- clusterEvalQ(cl, library(RPostgreSQL))

  foreach(i=1:length(all.sampleIDs)) %dopar% {
    if (isTest) {
      # nrow set for testing
      tbl.wellington <- readDataTable(sourcePath, all.sampleIDs[[i]], nrow = 10, 
                                          chromosome)
    } else {
      tbl.wellington <- readDataTable(sourcePath, all.sampleIDs[[i]], NA, 
                                          chromosome)
    }
    print("Data table read. Merging with Fimo...")
    # Connect to fimo
    fimo <- getDBConnection(fimo)
    tbl <- mergeFimoWithFootprints(tbl.wellington, all.sampleIDs[[i]], 
                                   dbConnection = fimo,
                                   method)
    dbDisconnect(fimo)
    
    print("Merged. Now splitting table to regions and hits...")
    x <- splitTableIntoRegionsAndHits(tbl, minid)
    #printf("filling %d regions, %d hits for %s", nrow(x$regions), 
    #       nrow(x$hits), all.sampleIDs[i])

    # Connect to footprint DB
    dbConnection <- getDBConnection(dbConnection)
    fillToDatabase(x$regions, x$hits, dbConnection, dbUser, dbTable)
    databaseSummary(dbConnection)
    dbDisconnect(dbConnection)
  } # for sampleID
} # fill.all.samples.by.chromosome
#-------------------------------------------------------------------------------
