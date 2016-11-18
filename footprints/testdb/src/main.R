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
  
  for(sampleID in all.sampleIDs){
    printf("---- %s (%s) (%d/%d)", sampleID, chromosome, 
           grep(sampleID, all.sampleIDs), length(all.sampleIDs))

    if (isTest) {
      # nrow set for testing
      tbl.wellington <- readDataTable(sourcePath, sampleID, nrow = 10, 
                                          chromosome)
    } else {
      tbl.wellington <- readDataTable(sourcePath, sampleID, NA, 
                                          chromosome)
    }
    print("Data table read. Merging with Fimo...")
    tbl <- mergeFimoWithFootprints(tbl.wellington, sampleID, 
                                   dbConnection = fimo,
                                   method)
    print("Merged. Now splitting table to regions and hits...")
    x <- splitTableIntoRegionsAndHits(tbl, minid)
    printf("filling %d regions, %d hits for %s", nrow(x$regions), 
           nrow(x$hits), sampleID)
    fillToDatabase(x$regions, x$hits, dbConnection, dbUser, dbTable)
    databaseSummary(dbConnection)
  } # for sampleID
} # fill.all.samples.by.chromosome
#-------------------------------------------------------------------------------
