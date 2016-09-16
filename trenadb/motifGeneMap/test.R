library(RPostgreSQL)
tbl.map <- read.table("motif_to_tf_mappings_with_tfclass_include_multiple.csv", sep=",", as.is=TRUE, header=TRUE)
dim(subset(tbl.map, tfs=="ELF1"))  # [1] 25  2
write.table(tbl.map, file="genesMotifs.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
db <- dbConnect(PostgreSQL(), user="pshannon", dbname="trena")

dbGetQuery(db, "select * from tfmotifs limit 3")

