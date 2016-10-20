library(RPostgreSQL)
library(RUnit)

printf <- function(...) print(noquote(sprintf(...)))

#-------------------------------------------------------------------------------
runTests <- function()
{
  test.getDBConnection()
} # runTests

#-------------------------------------------------------------------------------

getDBConnection <- function(driver = dbDriver("PostgreSQL"), user= "trenatest", 
                            password="trenatest", dbname="trenatest", 
                            host="whovian") {
  dbConnect(driver, user=user, password=password, dbname=dbname, host=host)
} # getDBConnection

#-------------------------------------------------------------------------------
test.getDBConnection <- function()
{
  printf("--- test.getDBConnection")
  
  db.testConnection <- dbConnect(PostgreSQL(), user= "trenatest", password="trenatest", dbname="trenatest", host="whovian")
  
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