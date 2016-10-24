

#-------------------------------------------------------------------------------
runTests <- function()
{
  test.getDBConnection()
  test.readWellingtonTable()
  test.mergeFootprintsWithFimo()
  test.splitTableIntoRegionsAndWellingtonHits()
  #test.fill.to.database()
  
  #test.combineFootprintsAndFimo()
  #test.combineFootprintsAndDatabasedFimo()
  
} # runTests
#-------------------------------------------------------------------------------
test.readWellingtonTable <- function()
{
  printf("--- test.readWellingtonTable")
  
  tbl <- readWellingtonTable(wellington.path, test.sampleID, 5, "chr19")
  checkEquals(dim(tbl), c(5,6))
  checkEquals(colnames(tbl), c("chrom", "start", "end", "name", "score", "strand"))
  checkEquals(unique(tbl$chrom), "chr19")
  
  # now read without chrom or nrow constraints
  tbl <- readWellingtonTable(wellington.path, test.sampleID)
  checkEquals(ncol(tbl), 6)
  checkTrue(nrow(tbl) > 300)
  checkEquals(colnames(tbl), c("chrom", "start", "end", "name", "score", "strand"))
  checkEquals(head(sort(unique(tbl$chrom))), "chr19")
  
} # test.readWellingtonTable
#-------------------------------------------------------------------------------
test.mergeFootprintsWithFimo <- function()
{
  printf("--- test.mergeFootprintsWithFimo")
  tbl.fp <- readWellingtonTable(wellington.path, test.sampleID, nrow=3, "chr19")
  tbl <- mergeFimoWithFootprints(tbl.fp, test.sampleID)
  checkEquals(ncol(tbl), 12)
  checkEquals(sort(colnames(tbl)),
              c("chrom", "fimo.pvalue", "fimo.score", "loc", "method",
                "motif", "motif.end", "motif.sequence", "motif.start", "motif.strand", "sample_id",
                "wellington.score"))
  checkTrue(nrow(tbl) >= 8)
  
  duplicated.loc <- "chr19:636760-636770"
  checkEquals(length(grep(duplicated.loc, tbl$loc)), 5)
  #  3 distinct motifs mapped to this region
  checkEquals(sort(subset(tbl, loc == duplicated.loc)$motif),
              c("MA0003.3", "MA0812.1", "MA0812.1", "MA0814.1", "MA0814.1"))
  
  invisible(tbl)
  
} # test.mergeFootprintsWithFimo
#-------------------------------------------------------------------------------
test.splitTableIntoRegionsAndWellingtonHits <- function(tbl)
{
  printf("--- test.splitTableIntoRegionsAndWellingtonHits")
  tbl.fp <- readWellingtonTable(wellington.path, test.sampleID, nrow=3, "chr19")
  tbl <- mergeFimoWithFootprints(tbl.fp, test.sampleID)
  
  x <- splitTableIntoRegionsAndWellingtonHits(tbl, "minid")
  checkEquals(sort(names(x)), c("hits", "regions"))
  checkEquals(colnames(x$regions), region.schema()) #
  checkEquals(colnames(x$hits), hit.schema())
  
} # test.splitTableIntoRegionsAndWellingtonHits    
#-------------------------------------------------------------------------------
test.fill.to.database <- function()
{
  printf("--- test.fill.to.database")
  if(!exists("db.wellington.test"))
    db.wellington.test <- 
      getDBConnection(user="trenatest", 
                      password="trenatest", 
                      dbname="testwellington")
  
  createEmptyDatabaseTables('trenatest', 'testwellington')
  knownLocs <<- new.env(parent=emptyenv())
  
  tbl.fp <- readWellingtonTable(wellington.path, test.sampleID, nrow=3, "chr19")
  tbl <- mergeFimoWithFootprints(tbl.fp, test.sampleID)
  x <- splitTableIntoRegionsAndWellingtonHits(tbl, "minid")
  
  fill.to.database(x$regions, x$hits, db.wellington.test)
  checkEquals(sort(dbListTables(db.wellington.test)), c("hits", "regions"))
  
  tbl.regions <- dbGetQuery(db.wellington.test, "select * from regions")
  checkEquals(dim(tbl.regions), c(2, 4))
  checkTrue(all(region.schema() == colnames(tbl.regions)))
  
  tbl.hits <- dbGetQuery(db.wellington.test, "select * from hits")
  checkTrue(all(hit.schema() == colnames(tbl.hits)))
  checkEquals(dim(tbl.hits), c(8, 14))
  checkEquals(unique(tbl.hits$type), "motif.in.footprint")
  checkTrue(all(tbl.hits$loc %in% tbl.regions$loc))
  
} # test.fill.to.database
#-------------------------------------------------------------------------------
test.getDBConnection <- function()
{
  printf("--- test.getDBConnection")
  
  db.testConnection <- dbConnect(
    PostgreSQL(), 
    user= "trenatest", 
    password="trenatest", 
    dbname="trenatest", 
    host="whovian")
  
  # test the class of db.testConnection
  checkEquals(class(db.testConnection)[1], "PostgreSQLConnection")
  
  # test that we can read tables from db.testConnection
  checkEquals(dbExistsTable(db.testConnection, "testtable"), TRUE)
  
  # test that we can write to db.testConnection
  checkEquals(
    dbWriteTable(conn = db.testConnection,
                 name = "testtable",
                 value=data.frame("teststring"), 
                 append = TRUE, 
                 row.names = FALSE),
    TRUE)
  
  # test that we can read from db.testConnection
  checkEquals(
    is.data.frame(
      dbGetQuery(db.testConnection, "select * from testtable")
    ), TRUE)
  
  # clean up and disconnect
  dbGetQuery(db.testConnection, "DELETE FROM testtable")
  dbDisconnect(db.testConnection)
} # test.getDBConnection
#-------------------------------------------------------------------------------
test_examine.region <- function()
{
  printf("--- test_examine.region")
  sampleIDs <- c("ENCSR000EJB")
  chrom <- "chr19"
  start <- 44903773
  end <- 44903785
  examine.region(chrom, start, end, sampleIDs)
  
} # test_examine.region
#-------------------------------------------------------------------------------