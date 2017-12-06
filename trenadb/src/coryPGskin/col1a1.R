library(igvR)
library(RPostgreSQL)
source("../utils.R")
igv <- igvR()

#db.w <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="skin_wellington", host="whovian")
#db.h <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="skin_hint", host="whovian")

db.w.lympho <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="lymphoblast_wellington", host="whovian")
db.h.lympho <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="lymphoblast_hint", host="whovian")

db.w.brain <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="brain_wellington", host="whovian")
db.h.brain <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="brain_hint", host="whovian")


chrom <- "chr17"
tss <- 50201632
start <- tss - 5000
end   <- tss + 5000

#tbl.w <- getHits(db.w, chrom, start, end)
#tbl.bed.w <- unique(tbl.w[, c("chrom", "start", "endpos", "name",  "score2")])
#displayBedTable(igv, tbl.bed.w, "wellington")

#tbl.h <- getHits(db.h, chrom, start, end)
#tbl.bed.h <- unique(tbl.h[, c("chrom", "start", "endpos", "name",  "score2")])
#displayBedTable(igv, tbl.bed.h, "hint")
 
#tbl.w.lympho <- getHits(db.w.lympho, chrom, start, end)
#tbl.bed.w.lympho <- unique(tbl.w.lympho[, c("chrom", "start", "endpos", "name",  "score2")])
#displayBedTable(igv, tbl.bed.w.lympho, "lympho_wellington")

#tbl.h.lympho <- getHits(db.h.lympho, chrom, start, end)
#tbl.bed.h.lympho <- unique(tbl.h.lympho[, c("chrom", "start", "endpos", "name",  "score2")])
#displayBedTable(igv, tbl.bed.h.lympho, "lympho_hint")

tbl.w.brain <- getHits(db.w.brain, chrom, start, end)
tbl.bed.w.brain <- unique(tbl.w.brain[, c("chrom", "start", "endpos", "name",  "score2")])
displayBedTable(igv, tbl.bed.w.brain, "brain_wellington")

tbl.h.brain <- getHits(db.h.brain, chrom, start, end)
tbl.bed.h.brain <- unique(tbl.h.brain[, c("chrom", "start", "endpos", "name",  "score2")])
displayBedTable(igv, tbl.bed.h.lympho, "brain_hint")

               
