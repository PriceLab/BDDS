#-------------------------------------------------------------------------------
getDBConnection <- function(driver = dbDriver("PostgreSQL"), user= "trenatest", 
                            password="trenatest", dbname="trenatest", 
                            host="whovian") {
  dbConnect(driver, user=user, password=password, dbname=dbname, host=host)
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
databaseSummary <- function()
{
  region.count <- dbGetQuery(db.wellington, "select count(*) from regions")[1,1]
  hit.count <- dbGetQuery(db.wellington, "select count(*) from hits")[1,1]
  printf("%d hits in %d regions", hit.count, region.count)
  
} # databaseSummary
#-------------------------------------------------------------------------------
## NOTE: original sql command started with \connect wholeBrain-wellington - needed?
createEmptyDatabaseTables <- function(dbuser, dbConnection)
{
  sql_command <- '
  drop table regions;
  drop table hits;
  
  create table regions(loc varchar primary key,
  chrom varchar,
  start int,
  endpos int);
  
  grant all on table "regions" to trena;
  
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
  
  grant all on table "hits" to ' + dbuser + ';'

  dbGetQuery(dbConnection, sql_command)
} # createEmptyDatabaseTables
#-------------------------------------------------------------------------------
## NOTE: original sql command started with \connect wholeBrain-wellington - needed?
appendToRegionsTable <- function(tbl, dbConnection)
{
  write.table(tbl, file="regions.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", na="NULL")
  sql_command <- "
  \copy regions from 'regions.tsv' delimiter E'\t' CSV NULL as 'NULL';"
  dbGetQuery(dbConnection, sql_command)
  
} # appendToRegionsTable
#-------------------------------------------------------------------------------
## NOTE: original sql command started with \connect wholeBrain-wellington - needed?
appendToHitsTable <- function(tbl, dbConnection)
{
  write.table(tbl, file="hits.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", na="NULL")
  sql_command <- "
  \copy hits from 'hits.tsv' delimiterE'\t' CSV NULL as 'NULL'"
  dbGetQuery(dbConnection, sql_command)
  
} # appendToHitsTable
#-------------------------------------------------------------------------------
fill.to.database <- function(tbl.regions, tbl.hits, dbConnection)
{
  appendToRegionsTable(tbl.regions, dbConnection)
  appendToHitsTable(tbl.hits, dbConnection)
}