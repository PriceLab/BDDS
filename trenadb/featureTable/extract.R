library(RPostgreSQL)
library(GenomicRanges)
library(RUnit)
source("../regionAndHitsSchemas.R")
#------------------------------------------------------------------------------------------------------------------------
if(!exists("db.chipseq"))
   db.chipseq <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="chipseq", host="whovian")

if(!exists("db.gtf"))
   db.gtf <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="gtf", host="whovian")

if(!exists("apoe")){
   tbl.tmp <- dbGetQuery(db.gtf, "select * from hg38human where gene_name='APOE' and moleculetype='gene'")
   apoe <- list(chrom=tbl.tmp[1, "chr"], start=tbl.tmp[1, "start"])
   }

if(!exists("db.hint"))
   db.hint <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="hint", host="whovian")

if(!exists("db.trena"))
  db.trena <- dbConnect(PostgreSQL(), user="pshannon", dbname="trena")
      
if(!exists("tbl.genesmotifs"))
    tbl.genesmotifs <- dbGetQuery(db.trena, "select * from tfmotifs")

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test.addFimoRegions()
  test.toFeatureTable()
  test.toFeatureTable.big()
  test.hintToFeatureTable()
  
} # runTests
#------------------------------------------------------------------------------------------------------------------------
getHits <- function(db, chrom, start, end)
{
   query.p0 <- "select loc, chrom, start, endpos from regions"
   query.p1 <- sprintf("where chrom='%s' and start > %d and endpos < %d", chrom, start, end)
   query.regions <- paste(query.p0, query.p1)
   tbl.regions <- dbGetQuery(db, query.regions)
   if(nrow(tbl.regions) == 0)
       return(data.frame())
   loc.set <- sprintf("('%s')", paste(tbl.regions$loc, collapse="','"))
   query.hits <- sprintf("select * from hits where loc in %s", loc.set)
   tbl.hits <- dbGetQuery(db, query.hits)
   merge(tbl.regions, tbl.hits, on="loc")

} # getHits
#------------------------------------------------------------------------------------------------------------------------
getFimoHits <- function(chrom, start, end)
{
   chrom <- sub("^chr", "", chrom)
   query <- sprintf("select * from fimo_hg38 where chrom='%s' and start > %d and endpos < %d", chrom, start, end)
   dbGetQuery(db.trena, query)

} # getFimoHits
#------------------------------------------------------------------------------------------------------------------------
# find all overlaps between the 151 base pair chipseq regions, and short motif-based fimo regions
# then expand the tbl.cs by joining it with all the tbl.fimo regions which overlap with each chipseq region
# then whittle those down, keeping only those fimo/chipseq rows where the fimo motif is in fact
# associated with the transcript factor reported in tbl.cs
addFimoRegions <- function(tbl.cs, tbl.fimo)
{
   printf("--- entering addFimoRegions, %d cs, %d fimo", nrow(tbl.cs), nrow(tbl.fimo))
    
   gr.fimo <- with(tbl.fimo, GRanges(seqnames=chrom, IRanges(start=start, end=endpos)))
   gr.cs   <- with(tbl.cs,   GRanges(seqnames=chrom, IRanges(start=start, end=endpos)))
   tbl.overlaps <- as.data.frame(findOverlaps(gr.fimo, gr.cs, type="any"))
   tbl.combined <- cbind(tbl.fimo[tbl.overlaps$queryHits, ], tbl.cs[tbl.overlaps$subjectHits, c("type", "name", "score1")])
   tbl.combined <- unique(tbl.combined)

   sharedMotif <- function(tf, motif){
      motifs.for.this.gene <- subset(tbl.genesmotifs, gene==tf)$motif
      found <- motif %in% motifs.for.this.gene
      return(found)
      }

   tbl.combined$bindingSite <- sapply(1:nrow(tbl.combined),
                                   function(i) sharedMotif(tbl.combined$name[i], tbl.combined$motifname[i]))

   rownames(tbl.combined) <- NULL

   printf("tbl.combined: %d %d", nrow(tbl.combined), ncol(tbl.combined))
   
      # discard the rows where the fimo motif is not associated with the tf
   tbl.out <- subset(tbl.combined, bindingSite==TRUE)

      # add a loc & length field
   tbl.out$loc <- with(tbl.out, sprintf("%s:%d-%d", chrom, start, endpos))
   tbl.out$length <- with(tbl.out, 1 + endpos - start)

   printf("tbl.out: %d %d", nrow(tbl.out), ncol(tbl.out))

   invisible(tbl.out)

} # addFimoRegions
#------------------------------------------------------------------------------------------------------------------------
test.addFimoRegions <- function()
{
   printf("--- test.addFimoRegions")
   shoulder <- 1000
   tbl.apoe.chipseq <- getHits(db.chipseq, apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
   tbl.fimo <- getFimoHits(apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
   tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")
   tbl.expanded <- addFimoRegions(tbl.apoe.chipseq, tbl.fimo)
   expected.colnames <- c("motifname", "chrom", "start", "endpos", "strand", "motifscore", "pval",
                          "empty", "sequence", "type", "name", "score1", "bindingSite", "loc", "length")
   checkTrue(all(expected.colnames %in% colnames(tbl.expanded)))

} # test.addFimoRegions
#------------------------------------------------------------------------------------------------------------------------
# annotateWithMotifs <- function(tbl.cs, tbl.fimo)
# {
#    gr.fimo <- with(tbl.fimo, GRanges(seqnames=chrom, IRanges(start=start, end=endpos)))
#    gr.cs   <- with(tbl.cs,   GRanges(seqnames=chrom, IRanges(start=start, end=endpos)))
#    tbl.overlaps <- as.data.frame(findOverlaps(gr.fimo, gr.cs, type="any"))
#    tbl.combined <- cbind(tbl.cs[tbl.overlaps$subjectHits,],
#                          tbl.fimo[tbl.overlaps$queryHits,])[, c("loc", "name", "motifname")]
#    tbl.combined <- unique(tbl.combined)
# 
#    sharedMotif <- function(gene, motif){
#       motifs.for.this.gene <- dbGetQuery(db.trena, sprintf("select motif from tfmotifs where gene='%s'", gene))$motif
#       found <- motif %in% motifs.for.this.gene
#       #printf("gene: %s   has motifs: %s  includes? %s: %s", gene, paste(motifs.for.this.gene, collapse=","), motif, found
#      return(found)
#      }
#
#   tbl.combined$bindingSite <- sapply(1:nrow(tbl.combined),
#                                   function(i) sharedMotif(tbl.combined$name[i], tbl.combined$motif[i]))
#
#    rownames(tbl.combined) <- NULL
#    subset(tbl.combined, bindingSite==TRUE)
#    tbl.status <- subset(tbl.combined, bindingSite==TRUE)
# 
#    failed.tfs <- setdiff(tbl.cs$name, tbl.status$name)
# 
#    if(length(failed.tfs) > 0){
#       tbl.failed <-  subset(tbl.cs, name %in% failed.tfs)[, c("loc", "name")]
#       tbl.failed$motifname <- "noMotif"
#       tbl.failed$bindingSite <- FALSE
#       tbl.status <- rbind(tbl.status, tbl.failed)
#       }
# 
#      # now add bindingSite info
# 
#    tbl.out <- merge(tbl.cs, tbl.status, by=c("loc", "name"))
#    tbl.out[, c(hit.schema(), "bindingSite", "motifname")]
# 
#    invisible(tbl.out)
#    
# } # annotateWithMotifs
#------------------------------------------------------------------------------------------------------------------------
test.annotateWithMotifs <- function()
{
   printf("--- test.annotateWithMotifs")
   shoulder <- 2000
   tbl.apoe <- getHits(db.chipseq, apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
   tbl.fimo <- getFimoHits(apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
   tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")
   tbl.anno <- annotateWithMotifs(tbl.apoe, tbl.fimo)
   checkEquals(colnames(tbl.anno), c("loc", "name", "chrom", "start", "endpos", "type", "length", "strand", "sample_id",
                                     "method", "provenance", "score1", "score2", "score3", "score4", "score5", "score6",
                                     "motifname", "bindingSite"))
   checkTrue(nrow(tbl.anno) >= nrow(tbl.apoe))  # should be true of any chipseq table.
   
   shoulder <- 5000
   tbl.apoe <- getHits(db.chipseq, apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
   tbl.fimo <- getFimoHits(apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
   tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")
   tbl.anno <- annotateWithMotifs(tbl.apoe, tbl.fimo)
   checkEquals(colnames(tbl.anno), c("loc", "name", "chrom", "start", "endpos", "type", "length", "strand", "sample_id",
                                     "method", "provenance", "score1", "score2", "score3", "score4", "score5", "score6",
                                     "motifname", "bindingSite"))
   checkTrue(nrow(tbl.anno) >= nrow(tbl.apoe))  # should be true of any chipseq table.

} # test.annotateWithMotifs
#------------------------------------------------------------------------------------------------------------------------
toFeatureTable <- function(tbl.hits)
{
   printf("--- entering toFeatureTable, tbl.hits is (%d, %d)", nrow(tbl.hits), ncol(tbl.hits))

   motifNames <- paste("chipseq", sort(unique(tbl.hits$motifname)), sep="_")
   column.names <- c("uLoc", motifNames)
   uLocs <- unique(tbl.hits$loc)
   tbl <- data.frame(matrix(data=0, nrow=length(uLocs), ncol=length(column.names)), stringsAsFactors=FALSE)
   colnames(tbl) <- column.names
   tbl$uLoc <- uLocs
   for(r in 1:nrow(tbl.hits)){
      row <- tbl.hits[r, "loc"]
      col <- sprintf("chipseq_%s", tbl.hits[r, "motifname"])
      tbl[grep(row, tbl$uLoc), col] <- tbl[grep(row, tbl$uLoc), col] + 1
      }

   tbl

} # toFeatureTable
#------------------------------------------------------------------------------------------------------------------------
test.toFeatureTable <- function(shoulder=1000)
{
   printf("--- test.toFeatureTable")

       # chipseq data: one 151 base pair hit, claims that this is shared by 3 tfs: CTCF, RUNX3, PBX3
   tbl.apoe <- getHits(db.chipseq, apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
       # fimo identifies 511 binding sites, among 164 unique motifs
   tbl.fimo <- getFimoHits(apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
       # our current fimo database lists chromosome names as, e.g., "10".  standardize them
   tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")

   tbl.expanded <- addFimoRegions(tbl.apoe, tbl.fimo)
   tbl <- toFeatureTable(tbl.expanded)

      # we expect one entry in the feature table for every row in tbl.expanded
   printf("found %d loc/motif hits from tbl.expanded of %d rows", sum(tbl[, -1]), nrow(tbl.expanded))
   checkEquals(nrow(tbl.expanded), sum(tbl[, -1]))

      # for every row, identify each column with a 1, make sure tbl.expanded[

   for(r in 1:nrow(tbl)){
      hits <- which(tbl[r,] == 1)
      this.loc <- tbl$uLoc[r]
      if(length(hits) > 0){
         column.names <- sub("chipseq_", "", colnames(tbl)[hits])
         #if(length(column.names) > 10) browser()
         #printf("loc: %s   column.names: %s", this.loc, paste(column.names, collapse=","))
         #printf("checked out? %s (%d,%d)", 
         #       nrow(subset(tbl.expanded, motifname %in% column.names & loc==this.loc)) == length(hits),
         #       nrow(subset(tbl.expanded, motifname %in% column.names & loc==this.loc)), length(hits))
                
         checkEquals(nrow(subset(tbl.expanded, motifname %in% column.names & loc==this.loc)), length(hits))
         } # if hits
      } # for r

   #browser();
   #x <- 99

} # test.toFeatureTable
#------------------------------------------------------------------------------------------------------------------------
test.toFeatureTable.big <- function(shoulder=15000)
{
   printf("--- test.toFeatureTable.big")

       # chipseq data: one 151 base pair hit, claims that this is shared by 3 tfs: CTCF, RUNX3, PBX3

   #browser()
   printf("   shoulder: %d", shoulder)
   tbl.apoe <- getHits(db.chipseq, apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
   printf("   tbl.apoe: %d rows", nrow(tbl.apoe))

       # fimo identifies 511 binding sites, among 164 unique motifs
   tbl.fimo <- getFimoHits(apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
   printf("   tbl.fimo: %d rows", nrow(tbl.fimo))

       # our current fimo database lists chromosome names as, e.g., "10".  standardize them
   tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")

   tbl.expanded <- addFimoRegions(tbl.apoe, tbl.fimo)
   printf("    testing freshly created 'tbl.expanded': %d x %d", nrow(tbl.expanded), ncol(tbl.expanded))
   tbl <- toFeatureTable(tbl.expanded)
   checkEquals(nrow(tbl.expanded), sum(tbl[, -1]))
   x <- 99

} # test.toFeatureTable.big
#------------------------------------------------------------------------------------------------------------------------
hintToFeatureTable <- function(tbl.hits)
{

   motifNames <- paste("hint", sort(unique(tbl.hits$name)), sep="_")
   column.names <- c("uLoc", motifNames)
   uLocs <- unique(tbl.hits$loc)
   tbl <- data.frame(matrix(data=0, nrow=length(uLocs), ncol=length(column.names)), stringsAsFactors=FALSE)
   colnames(tbl) <- column.names
   tbl$uLoc <- uLocs
   for(r in 1:nrow(tbl.hits)){
      row <- tbl.hits[r, "loc"]
      col <- sprintf("hint_%s", tbl.hits[r, "name"])
      tbl[grep(row, tbl$uLoc), col] <- tbl[grep(row, tbl$uLoc), col] + 1
      }

   invisible(tbl)

} # hintToFeatureTable
#------------------------------------------------------------------------------------------------------------------------
test.hintToFeatureTable <- function(shoulder=1000)
{
   #tbl.apoe <- getHits(db.hint, apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
   load("tbl.apoe.hint.RData")
   tbl.hits <- tbl.apoe
   printf("   tbl.hits: %d rows", nrow(tbl.hits))

   tbl <- hintToFeatureTable(tbl.hits)
   checkEquals(nrow(tbl.hits), sum(tbl[, -1]))

      # for every row, identify each column with a 1, make sure tbl.expanded[

   for(r in 1:nrow(tbl)){
      hits <- which(tbl[r,] == 1)
      this.loc <- tbl$uLoc[r]
      if(length(hits) > 0){
         column.names <- sub("hint_", "", colnames(tbl)[hits])
         #if(length(column.names) > 10) browser()
         #printf("loc: %s   column.names: %s", this.loc, paste(column.names, collapse=","))
         #printf("checked out? %s (%d,%d)", 
         #       nrow(subset(tbl.hits, name %in% column.names & loc==this.loc)) == length(hits),
         #       nrow(subset(tbl.hits, name %in% column.names & loc==this.loc)), length(hits))
         checkEquals(nrow(subset(tbl.hits, name %in% column.names & loc==this.loc)), length(hits))
         } # if hits
      } # for r


} # test.hintToFeatureTable
#------------------------------------------------------------------------------------------------------------------------
ensembl <- function(chrom, start, end)
{
   #browser()
   tbl.hits.cs <- getHits(db.chipseq, chrom, start, end)
   tbl.fimo.cs <- getFimoHits(chrom, start, end)
   tbl.fimo.cs$chrom <- paste("chr", tbl.fimo.cs$chrom, sep="")

   tbl.expanded.cs <- addFimoRegions(tbl.hits.cs, tbl.fimo.cs)
   tbl.feature.cs <- toFeatureTable(tbl.expanded.cs)

   tbl.hits.hint <- getHits(db.hint, chrom, start, end)
   tbl.feature.hint <- hintToFeatureTable(tbl.hits.hint)
   printf ("uLoc overlap: %d", length(intersect(tbl.feature.hint$uLoc, tbl.feature.cs$uLoc)))
   x <- tbl.feature.cs$uLoc; y <- tbl.feature.hint$uLoc;
   printf("chipseq regions: %d", length(x))
   printf("   HINT regions: %d", length(y))
   browser()
   x <- 99
    
} # ensembl
#------------------------------------------------------------------------------------------------------------------------
test.ensembl <- function()
{
   printf("--- test.ensembl")
   shoulder <- 1000
   tbl <- ensembl(apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
   
} # test.ensembl
#------------------------------------------------------------------------------------------------------------------------
