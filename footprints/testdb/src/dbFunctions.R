#-------------------------------------------------------------------------------
getDBConnection <- function(database) {
  port = "5432"
  driver = dbDriver("PostgreSQL")

  if (database == "trenatest") {
    user= "ben"
    password="ben_PASS"
    dbname="trenatest" 
    host="bdds-rds.globusgenomics.org"
    
  } else if (database == "fimo") {
    user= "ben"
    password="ben_PASS"
    dbname="fimo"
    host="bdds-rds.globusgenomics.org"
    
  } else if (database == "testwellington") {
    user = "ben"
    password = "ben_PASS"
    dbname = "testwellington"
    host = "bdds-rds.globusgenomics.org"
    
  } else if (database == "wellington") {
    user= "ben"
    password="ben_PASS"
    dbname="wellington" 
    host="bdds-rds.globusgenomics.org"
    
  } else if (database == "trenatest_whovian") {
    user= "trenatest"
    password="trenatest"
    dbname="trenatest" 
    host="whovian"
    
  } else if (database == "fimo_whovian") {
    user= "trena" 
    password="trena" 
    dbname="fimo"
    host="whovian"
    
  } else if (database == "testwellington_whovian") {
    user = "trenatest"
    password = "trenatest"
    dbname = "testwellington"
    host = "whovian"
    
  } else if (database == "wellington_whovian") {
    user = "trena"
    password = "trena"
    dbname = "wholeBrain-wellington"
    host = "whovian"
  }
  
  dbConnect(drv=driver, user=user, password=password, dbname=dbname, host=host, 
            port=port)
} # getDBConnection

#-------------------------------------------------------------------------------
region.schema <- function()
{
  c("loc", "chrom", "start", "endpos")
} # region.schema
#-------------------------------------------------------------------------------
hit.schema <- function()
{
  c("loc", "type", "name", "length", 
    "strand", "sample_id", "method", "provenance",
    "score1", "score2", "score3", "score4", "score5", "score6")
} # hit.schema
#-------------------------------------------------------------------------------
databaseSummary <- function(dbConnection = "db.wellington.test")
{
  region.count <- dbGetQuery(dbConnection, "select count(*) from regions")[1,1]
  hit.count <- dbGetQuery(dbConnection, "select count(*) from hits")[1,1]
  printf("%d hits in %d regions", hit.count, region.count)
} # databaseSummary
#-------------------------------------------------------------------------------
createEmptyDatabaseTables <- function(dbUser, dbName, dbConnection= "db.wellington.test")
{
  sql_command <- paste('drop table regions;
  drop table hits;
  
  create table regions(loc varchar primary key,
  chrom varchar,
  start int,
  endpos int);
  
  grant all on table "regions" to ', dbUser, ';', '
  
  create table hits(loc varchar,
  type varchar,
  name varchar,
  length int,
  strand char(1),
  sample_id varchar,
  method varchar,
  provenance varchar,
  score1 real,
  score2 real,
  score3 real,
  score4 real,
  score5 real,
  score6 real);
  
  grant all on table "hits" to ', dbUser, ';', sep="")

  dbGetQuery(dbConnection, sql_command)
} # createEmptyDatabaseTables
#-------------------------------------------------------------------------------
appendToRegionsTable <- function(tbl, dbConnection="db.wellington.test")
{
  dbWriteTable(dbConnection, "regions", tbl, row.names=FALSE, append=TRUE)
} # appendToRegionsTable
#-------------------------------------------------------------------------------
appendToHitsTable <- function(tbl, dbConnection="db.wellington.test")
{
  dbWriteTable(dbConnection, "hits", tbl, row.names=FALSE, append=TRUE)
} # appendToHitsTable
#-------------------------------------------------------------------------------
fillToDatabase <- function(tbl.regions, tbl.hits, dbConnection="db.wellington.test")
{
  appendToRegionsTable(tbl.regions, dbConnection)
  appendToHitsTable(tbl.hits, dbConnection)
}
