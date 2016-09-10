library(RPostgreSQL)
library(FimoClient)
library(getDNAClient)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("db.hint"))
   db.hint <- dbConnect(PostgreSQL(), user="pshannon", dbname="hintTest2")

if(!exists("db.fimo"))
   db.fimo <- dbConnect(PostgreSQL(), user="pshannon", dbname="hg38")

if(!exists("db.chipSeq"))
   db.chipSeq <- dbConnect(PostgreSQL(), user="pshannon", dbname="chipseqTest")

if(!exists("dna.service"))
    dna.service <- getDNAClient("hg38")

if(!exists("fimo.service"))
   fimo.service <-  FimoClient("whovian", 5558)

#------------------------------------------------------------------------------------------------------------------------
getHits <- function(db, chrom, start, stop)
{
   query.p0 <- "select loc, chrom, start, stop from regions"
   query.p1 <- sprintf("where chrom='%s' and start > %d and stop < %d", chrom, start, stop)
   query.regions <- paste(query.p0, query.p1)
   browser()
   tbl.regions <- dbGetQuery(db, query.regions)
   if(nrow(tbl.regions) == 0)
       return(data.frame())
   loc.set <- sprintf("('%s')", paste(tbl.regions$loc, collapse="','"))
   query.hits <- sprintf("select * from hits where loc in %s", loc.set)
   tbl.hits <- dbGetQuery(db, query.hits)
   merge(tbl.regions, tbl.hits, on="loc")

} # getHits
#------------------------------------------------------------------------------------------------------------------------
getFimoHits <- function(chrom, start, stop)
{
   chrom <- sub("^chr", "", chrom)
   query <- sprintf("select * from hg38 where chr='%s' and start > %d and endpos < %d", chrom, start, stop)
   dbGetQuery(db.fimo, query)

} # getFimoHits
#------------------------------------------------------------------------------------------------------------------------

chrom <- "chr21"
loc.start <- 25881625
loc.end   <- 25882287

tbl.cs <- getHits(db.chipSeq, chrom, loc.start, loc.end)
tbl.cs.ctcf <- subset(tbl.cs, name=="CTCF")
seq <- getSequenceByLoc(dna.service, chrom, loc.start, loc.end)
tbl.fimo1 <- requestMatch(fimo.service, list(seq))
tbl.fimo2 <- getFimoHits(chrom, loc.start, loc.end)
