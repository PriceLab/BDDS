#-------------------------------------------------------------------------------
fillAllSamplesByChromosome <- function(dbConnection = "db.wellington",
                                       fimo = "db.fimo",
                                       chromosome = "chr19",
                                       minid = "temp.filler.minid")
{
  knownLocs <<- new.env(parent=emptyenv())
  
  all.sampleIDs <- unlist(lapply(strsplit(list.files(wellington.path, 
                                                     "ENCSR.*.bed$"), 
                                          ".", fixed=TRUE), "[", 1))
  
  for(sampleID in all.sampleIDs){
    printf("---- %s (%s) (%d/%d)", sampleID, chromosome, 
           grep(sampleID, all.sampleIDs), length(all.sampleIDs))
    tbl.wellington <- readWellingtonTable(wellington.path, sampleID, NA, 
                                          chromosome)
    print("Wellington table read. Merging with Fimo...")
    tbl <- mergeFimoWithFootprints(tbl.wellington, sampleID, 
                                   dbConnection = db.fimo)
    print("Merged. Now splitting table to regions and hits...")
    x <- splitTableIntoRegionsAndWellingtonHits(tbl, minid)
    printf("filling %d regions, %d hits for %s", nrow(x$regions), 
           nrow(x$hits), sampleID)
    fillToDatabase(x$regions, x$hits, dbConnection)
    databaseSummary(dbConnection)
  } # for file
  
} # fill.all.samples.by.chromosome
#-------------------------------------------------------------------------------
