library(igvR)
library(RPostgreSQL)
source("../utils.R")
igv <- igvR()

db.w <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="skin_wellington", host="whovian")
db.h <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="skin_hint", host="whovian")

chrom <- "chr17"
tss <- 50201632
start <- tss - 1000
end   <- tss + 1000

tbl.w <- getHits(db.w, chrom, start, end)
tbl.bed.w <- unique(tbl.w[, c("chrom", "start", "endpos", "sample_id",  "score2")])
displayBedTable(igv, tbl.bed.w, "wellington")

tbl.h <- getHits(db.h, chrom, start, end)
tbl.bed.h <- unique(tbl.h[, c("chrom", "start", "endpos", "sample_id",  "score2")])
displayBedTable(igv, tbl.bed.h, "hint")
                
