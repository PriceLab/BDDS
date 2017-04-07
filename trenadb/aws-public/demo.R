library(RPostgreSQL)
db.hg38 <- dbConnect(PostgreSQL(), user= "galaxy", password="bdds_postgres_rds_pass*word", dbname="hg38",
                    host="bddsrds.globusgenomics.org")
dbListTables(db.hg38)
query <- "select * from gtf where moleculetype='gene' and gene_biotype='protein_coding'"
system.time(tbl.gtf <- dbGetQuery(db.hg38, query));     # 19797 x 30, 13 seconds

db.brainHint <- dbConnect(PostgreSQL(), user= "galaxy", password="bdds_postgres_rds_pass*word", dbname="brain_hint",
                    host="bddsrds.globusgenomics.org")
dbListTables(db.brainHint)
# [1] "hits"      "hits_1"    "hits_2"    "hits_3"    "regions"   "regions_1"
# [7] "regions_2" "regions_3"

