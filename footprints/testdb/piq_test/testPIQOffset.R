
library(RPostgreSQL)
library(dplyr)

# Read in the top 10000 lines of PIQ and sample 1000
tbl <- read.table(file = "/scratch/data/footprints/brain_piq_16/ENCSR000DBW.bed",
                  sep="\t", as.is=TRUE, nrows = 10000)
# tbl <- tbl[sample(nrow(tbl), 1000),]
colnames(tbl) <- c("chrom", "motif.start", "motif.end", "motif", "motif.strand",
                           "score1", "score2", "score3", "score4");
tbl <- select(tbl,1:4)

# Fix the motifs by removing the RC and add a column for the fimo end
tbl$motif <- sub('(MA\\d{4}).*','\\1',tbl$motif)
tbl$fimo.end <- NA

# Create the fimo db connection
fimo <- dbConnect(PostgreSQL(),
                  user = 'trena',
                  password = 'trena',
                  dbname = 'fimo',
                  host = 'localhost')

# For each entry in the table, query the database to get the fimo end point and save it
for(i in 1:length(tbl$chrom)){

    my.query <- paste("select endpos from fimo_hg38 where motifname like '%",
                              tbl$motif[i],"%' and start = ", tbl$motif.start[i],";",
                              sep="")

    fimo.end <- dbGetQuery(fimo, my.query)
    
    if(length(fimo.end$endpos) > 0){
        tbl$fimo.end[i] <- fimo.end$endpos[1]}
}

tbl %>% mutate(diff = motif.end - fimo.end) %>% summarise(avg = mean(diff,na.rm=TRUE),
                                                          total = n(),
                                                          nas = sum(is.na(diff)))
tbl %>% filter(is.na(fimo.end)) %>% summarise(diff.motifs = n_distinct(motif))
                  


