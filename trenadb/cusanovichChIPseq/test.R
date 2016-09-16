library(RPostgreSQL)
#------------------------------------------------------------------------------------------------------------------------
chrom <- "chr13"
loc.start <- 41019200
loc.end   <- 41019360

loc.string <- sprintf("%s:d-d", chrom, loc.start, loc.end)

db.cs <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="chipseq", host="whovian")
dbGetQuery(db.cs, "select * from hits limit 3")

query.regions <- sprintf("select * from regions where chrom='%s' and start > %d and endpos < %d", chrom, loc.start, loc.end)
system.time(tbl.regions <- dbGetQuery(db.cs, query.regions))  # 0.064 seconds
dim(tbl.regions)

query.hits <- sprintf("select * from hits where loc='%s'", tbl[1, "loc"])
system.time(tbl.hits <- dbGetQuery(db.cs, query.hits))     # 0.35 seconds
dim(tbl.hits)

tbl.out <- merge(tbl.regions, tbl.hits, on="loc")
preferred.column.order <- c("chrom", "start", "endpos", "name", "strand", "score1",
                            "type", "length", "sample_id", "method", "provenance",
                            "score2", "score3", "score4", "score5", "score6")
tbl.out <- tbl.out[, preferred.column.order]
head(tbl.out[, 1:10])

# now find the motif's mapped by fimo in this region
db.trena <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="trena", host="whovian")
# take a look
dbGetQuery(db.trena, "select * from fimo_hg38 limit 3")
query.fimo <- sprintf("select * from fimo_hg38 where chrom='%s' and start >= %d and endpos <= %d",
                      "13", loc.start, loc.end)
system.time(tbl.fimo <- dbGetQuery(db.trena, query.fimo))  # 5.17 seconds; index needed?

# what genes (tfs) have been mapped into the chr13:41019200-4109360 by cusanovich ChIPseq?
# which have motifs?
genes.tfs <- sort(unique(tbl.out$name))                                             # 36
genes.allMapped <- dbGetQuery(db.trena, "select distinct gene from tfMotifs")[,1]   # 847
genes.tfs.withMotifs <- intersect(genes.tfs, genes.allMapped)                       # 27/36

# what motifs are associated with each of these allegedly bound tfs?
motifs.csGenes <- lapply(genes.tfs.withMotifs, function(gene)
          dbGetQuery(db, sprintf("select motif from tfmotifs where gene = '%s'", gene))[,1])
names(motifs.csGenes) <- genes.tfs.withMotifs

# filter this list, keeping only motifs actually mapped by fimo in the target region
motifs.fimo <- sort(unique(tbl.fimo$motifname))
for(tf in names(motifs.csGenes)){
   found.by.fimo <- intersect(motifs.csGenes[[tf]], motifs.fimo)
   printf("%8s: %s", tf, paste(found.by.fimo, collapse=","))
   }
