library(RPostgreSQL)
library(GenomicRanges)
library(RUnit)
#library(igvR)
#------------------------------------------------------------------------------------------------------------------------
source("../../regionAndHitsSchemas.R")
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

if(!exists("db.wellington"))
   db.wellington <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="wellington", host="whovian")

if(!exists("db.piq"))
   db.piq <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="piq", host="whovian")

if(!exists("db.trena"))
   db.trena <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="trena", host="whovian")
      
if(!exists("tbl.genesmotifs"))
    tbl.genesmotifs <- dbGetQuery(db.trena, "select * from tfmotifs")

if(!exists("db.fimo"))
    db.fimo <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="fimo", host="whovian")

if(!exists("tbl.chipseq"))
    db.chipseq <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="chipseq", host="whovian")

#if(!exists("igv"))
#    igv <- igvR()

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test.createWellingtonTable()
   test.locStringToBedTable()
   test.addFimoRegions()
   test.chipseqToFeatureTable()
   test.ensemble()
   test.cleanFeatureTable()
  
} # runTests
#------------------------------------------------------------------------------------------------------------------------
getHits <- function(db, chrom, start, end)
{
   query.p0 <- "select loc, chrom, start, endpos from regions"
   query.p1 <- sprintf("where chrom='%s' and start >= %d and endpos <= %d", chrom, start, end)
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
# motif != NA used only for testing, to pull out interesting edge cases
createWellingtonTable <- function(chrom, start, end, motif=NA, sampleID=NA)
{
   tbl.hitsw <- getHits(db.wellington, chrom, start, end) # 6 x 17
   browser()
   if(!is.na(motif))
      tbl.hitsw <- subset(tbl.hitsw, name==motif)
   if(!is.na(sampleID))
      tbl.hitsw <- subset(tbl.hitsw, sample_id==sampleID)
   printf("found %d wellington hits in %d bases", nrow(tbl.hitsw), 1 + end - start)
   #displayBedTable(igv, tbl.hitsw[, c("chrom", "start", "endpos", "name", "score2")], "wellington")
   tbl.std <- tbl.hitsw[, c("loc",  "name", "length", "score1", "score2", "score3", "sample_id")]
   tbl.collapsed <- as.data.frame(table(tbl.std$loc, tbl.std$name))
   tbl.collapsed <- subset(tbl.collapsed, Freq != 0)
   colnames(tbl.collapsed) <- c("loc", "motif.w", "sample.count.w")

      # preallocate and zero fill the summary columns for each row

   distinct.hits.count <- nrow(tbl.collapsed)  # unique loc/motif combinations
   length <- rep(0, distinct.hits.count)
   #score1.median <- rep(0.0, distinct.hits.count)
   score1       <-   rep(0.0, distinct.hits.count)
   #score1.best <-   rep(0.0, distinct.hits.count)
   score2.median <- rep(0.0, distinct.hits.count)
   score2.best  <-   rep(0.0, distinct.hits.count)
   score3.median <- rep(0.0, distinct.hits.count)
   score3.best <-   rep(0.0, distinct.hits.count)
   length <-        rep(0,   distinct.hits.count)

   for(r in 1:nrow(tbl.collapsed)){
       this.loc <- tbl.collapsed$loc[r]
       this.motif <- tbl.collapsed$motif.w[r]
       tbl.sub <- subset(tbl.std, loc==this.loc & name==this.motif)
       printf("count of duplicated samples: %d, %s, %s", length(which(duplicated(tbl.sub$sample_id))),
              this.loc, this.motif)
       browser()
       score1[r]        <- max(tbl.sub$score1)
       score2.median[r] <- median(tbl.sub$score2)
       score2.best[r]   <- max(tbl.sub$score2)
       score3.median[r] <- median(tbl.sub$score3)
       score3.best[r]   <- min(tbl.sub$score3)
       length[r] <- median(tbl.sub$length) # should all be identical, but this covers all situations
       #printf("  found %d rows for %s %s", nrow(tbl.sub), this.loc, this.motif)
       x <- 99
       } # for r

   tbl.out <- cbind(tbl.collapsed, length, score1, score2.median, score2.best, score3.median, score3.best)
   colnames(tbl.out) <- c("loc", "motif.w", "samplecount.w", "length.w", "score1",
                          "score2.w.median",  "score2.w.best", "score3.w.median", "score3.w.best")
   tbl.out

}  # createWellingtonTable
#------------------------------------------------------------------------------------------------------------------------
test.createWellingtonTable <- function()
{
   printf("--- test.createWellingtonTable")

     # try with the development test set, a region with 118 wellington footprints in only 6 locs, across just 19 bases

   chrom <- "chr19"
   start <- 44903772
   end   <- 44903790

      # MA0813.1 hits on both strands, provides chance to test reduction to one hit only
   tbl.w <- createWellingtonTable(chrom, start, end, motif="MA0813.1", sampleID="ENCSR000EJE")
   checkEquals(dim(tbl.w), c(12, 10))
   checkEquals(length(unique(tbl.w$loc)), 7)
   checkEquals(length(unique(tbl.w$motif.w)), 12)

      # expand upstream and downstream by an extra 10kb

   tbl.w <- createWellingtonTable(chrom, start-10000, end+10000)
   checkEquals(dim(tbl.w), c(179, 10))
   checkEquals(length(unique(tbl.w$loc)), 139)

   # make sure that not all score3 medians are equal to score3 best
   checkTrue(length(which(!tbl.w$score3.w.median == tbl.w$score3.w.best)) > 0)
   
   # does a motif-rich, high-sample footprint produce a properly reduced table?
   
   tbl.wx <- createWellingtonTable("chr19", 45423910, 45423921)
   checkEquals(dim(tbl.wx), c(17, 10))
   checkEquals(dim(unique(tbl.wx)), c(17, 10))
   checkEquals(dim(unique(tbl.wx[, c("loc", "motif.w")])), c(17, 2))

}  # test.createWellingtonTable
#------------------------------------------------------------------------------------------------------------------------
createHintTable <- function(chrom, start, end)
{
   tbl.hitsh <- getHits(db.hint, chrom, start, end) # 6 x 17
   printf("found %d hint hits in %d bases", nrow(tbl.hitsh), 1 + end - start)
   #displayBedTable(igv, tbl.hitsw[, c("chrom", "start", "endpos", "name", "score2")], "hint")
   tbl.std <- tbl.hitsh[, c("loc",  "name", "length", "sample_id", "score1", "score2", "score3")]
   tbl.collapsed <- as.data.frame(table(tbl.std$loc, tbl.std$name))
   tbl.collapsed <- subset(tbl.collapsed, Freq != 0)
   colnames(tbl.collapsed) <- c("loc", "motif.h", "sample.count.h")

      # preallocate and zero fill the summary columns for each row

   distinct.hits.count <- nrow(tbl.collapsed)  # unique loc/motif combinations
   length <- rep(0, distinct.hits.count)
   score1.median <- rep(0.0, distinct.hits.count)
   score1.best <-   rep(0.0, distinct.hits.count)
   score2.median <- rep(0.0, distinct.hits.count)
   score2.best <-   rep(0.0, distinct.hits.count)
   score3.median <- rep(0.0, distinct.hits.count)
   score3.best <-   rep(0.0, distinct.hits.count)
   length <-        rep(0,   distinct.hits.count)

    # score1: hint score, ranges from 2 to 10,160, mean of ~123, median 44
    # score2: fimo score, -18.7 to 30.9, mean and median both about 11.5
    # score3: fimo pval, min 5.9e-13, max 1e-4, mean 4.3e-5, median 4.01 e-5

   for(r in 1:nrow(tbl.collapsed)){
       this.loc <- tbl.collapsed$loc[r]
       this.motif <- tbl.collapsed$motif.h[r]
       tbl.sub <- subset(tbl.std, loc==this.loc & name==this.motif)
       score1.median[r] <- median(tbl.sub$score1)
       score1.best[r]   <- max(tbl.sub$score1)
       score2.median[r] <- median(tbl.sub$score2)
       score2.best[r]   <- max(tbl.sub$score2)
       score3.median[r] <- median(tbl.sub$score3)
       score3.best[r]   <- min(tbl.sub$score3)
       length[r] <- median(tbl.sub$length) # should all be identical, but this covers all situations
       #printf("  found %d rows for %s %s", nrow(tbl.sub), this.loc, this.motif)
       x <- 99
       } # for r

   tbl.out <- cbind(tbl.collapsed, length, score1.median, score1.best, score2.median, score2.best,
                    score3.median, score3.best)
   colnames(tbl.out) <- c("loc", "motif.h", "samplecount.h", "length.h", "score1.h.median",  "score1.h.best",
                          "score2.h.median",  "score2.h.best", "score3.h.median", "score3.h.best")
   tbl.out

}  # createHintTable
#------------------------------------------------------------------------------------------------------------------------
test.createHintTable <- function()
{
   printf("--- test.createHintTable")

     # try with the development test set,

   chrom <- "chr19"
   start <- 44903772
   end   <- 44903790
   tbl.h <- createHintTable(chrom, start, end)
   checkEquals(dim(tbl.h), c(5, 10))
   checkEquals(length(unique(tbl.h$loc)), 3)
   checkEquals(length(unique(tbl.h$motif.h)), 5)

      # expand upstream and downstream by an extra 10kb
      # "found 995 hint hits in 20019 bases"

   tbl.h <- createHintTable(chrom, start-10000, end+10000)
   checkEquals(dim(tbl.h), c(641, 10))
   checkEquals(length(unique(tbl.h$loc)), 497)

   checkEquals(range(tbl.h$score1.h.best), c(6, 318))

     # make sure that not all score3 medians are equal to score3 best

   checkTrue(length(which(!tbl.h$score3.h.median == tbl.h$score3.h.best)) > 0)
   checkEquals(round(range(tbl.h$score2.h.best)), c(-14, 23))
   checkEqualsNumeric(min(tbl.h$score3.h.best), 2.46e-08)
   checkEqualsNumeric(max(tbl.h$score3.h.best), 9.99e-05)
   
}  # test.createHintTable
#------------------------------------------------------------------------------------------------------------------------
explore <- function()
{
   chrom <- apoe$chrom
   start <- apoe$start - 2000
   end  <- apoe$start + 500

   tbl.hitsw <- getHits(db.wellington, chrom, start, end) # 6 x 17
   tbl.hitsh <- getHits(db.hint, chrom, start, end)       # 97 17

   #tbl.hitsp <- getHits(db.piq, chrom, start, end)
   #load("tbl.piq.hits.chr19.1100bases.RData")
   load("tbl.piq.hits.chr19.2500bases.RData")
   

   tbl.chipseq <- getHits(db.chipseq, chrom, start, end)
   tbl.fimo <- getFimoHits(chrom, start, end)
   tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")
   tbl.csWithFimoAll <- addFimoRegions(tbl.chipseq, tbl.fimo)   # 76 15
   tbl.csWithFimoBound <- subset(tbl.csWithFimoAll, bindingSite==TRUE)
   tbl.csWithFimoBound$fullname <- paste(tbl.csWithFimoBound$name, tbl.csWithFimoBound$motifname, sep="-") # 5 16

   displayBedTable(igv, tbl.hitsw[, c("chrom", "start", "endpos", "name", "score2")], "wellington")
   displayBedTable(igv, tbl.hitsh[, c("chrom", "start", "endpos", "name", "score2")], "hint")
   displayBedTable(igv, tbl.csWithFimoBound[, c("chrom", "start", "endpos", "fullname", "score1")], "chipseq-fimo")
   tbl.piqFiltered <- subset(tbl.hitsp, score4 > 0.7)
   displayBedTable(igv, tbl.piqFiltered[, c("chrom", "start", "endpos", "name", "score4")], "piqFiltered")

   tbl.wellington <- tbl.hitsw[, c("loc",  "name", "length", "score1", "score2", "score3")]
   colnames(tbl.wellington) <- c("loc", "motif.w", "length",  "well1", "well2", "well3")

   tbl.chipseq <- tbl.csWithFimoBound[, c("loc", "name", "length", "score1", "motifscore", "pval")]
   colnames(tbl.chipseq) <- c("loc", "tf", "length",  "chip1", "chip2", "chip3")

   tbl.hint <- tbl.hitsh[, c("loc",  "name", "length", "score1", "score2", "score3")]
   colnames(tbl.hint) <-  c("loc", "motif.h", "length",  "hint1", "hint2", "hint3")
   
   tbl.piq <- tbl.hitsp[, c("chrom", "start", "endpos", "loc", "name", "length", "score1", "score2", "score3", "score4")]
   colnames(tbl.piq) <- c("chrom", "start", "endpos", "loc", "motif.p", "length", "piq1", "piq2", "piq3", "piq4")
   tbl.piq2 <- tbl.piq[order (tbl.piq$loc, tbl.piq$motif, tbl.piq$piq4, decreasing=TRUE),]
   dups <- which(duplicated(tbl.piq2[, c("loc", "motif.p")]))
   tbl.piq3 <- tbl.piq2[-dups,]
   tbl.piq4 <- refimo(tbl.piq3)

   #tbl.csw <- merge(tbl.chipseq, tbl.wellington, by=c("loc", "name", "length"), all=TRUE)
   #tbl.cswh <- merge(tbl.csw, tbl.hint,  by=c("loc", "name", "length"), all=TRUE)
   #tbl.cswhp <- merge(tbl.cswh, tbl.piq3, by=c("loc", "name", "length"), all=TRUE)

   tbl.csw <- merge(tbl.chipseq, tbl.wellington, by=c("loc", "length"), all=TRUE)
   tbl.cswh <- merge(tbl.csw, tbl.hint,  by=c("loc",  "length"), all=TRUE)
   tbl.cswhp <- merge(tbl.cswh, tbl.piq4, by=c("loc",  "length"), all=TRUE)
   
} # explore
#------------------------------------------------------------------------------------------------------------------------
# we want wellington, hint, chipseq and piq to all use the same fimo regions
# in contrast to the other methods, where i do the fimo intersections, piq does its own
# and apparently does so on an obsolete motif library.  and it may eliminate second-best hits
# case in point:
#   tbl.cswhp
# 13  chr19:44904896-44904914     19  CTCF    69 13.14750 8.31e-06     <NA>       NA       NA       NA                 <NA>    NA        NA       NA     <NA>       NA          NA         NA       NA
# 14  chr19:44904896-44904915     20  <NA>    NA       NA       NA     <NA>       NA       NA       NA                 <NA>    NA        NA       NA MA0139.1  8.12448 -0.19550100   5.687290 0.876744
# 15  chr19:44904899-44904917     19  CTCF    69  9.70492 4.36e-05     <NA>       NA       NA       NA                 <NA>    NA        NA       NA     <NA>       NA          NA         NA       NA
#
# solution:
#  getFimoHits("chr19", 44904895, 44904915)
#   motifname chrom    start   endpos strand motifscore     pval empty            sequence
# 1  MA0139.1    19 44904896 44904914      +    13.1475 8.31e-06       TGGCAGCCAGGGGGAGGTG
# 2  MA0155.1    19 44904900 44904911      +    12.5816 2.35e-05              AGCCAGGGGGAG
# 3  GTF2I.p2    19 44904904 44904912      +    13.5610 1.50e-05                 AGGGGGAGG
# choose motif with minium pval, thus getting the best full-length character match
#
#  tbl.refimo <- getFimoHits("chr19", 44904896-2, 44904915+2)
#  tbl.refimo[which(tbl.refimo$pval == min(tbl.refimo$pval)),]
#   motifname chrom    start   endpos strand motifscore     pval empty            sequence
# 2  MA0139.1    19 44904896 44904914      +    13.1475 8.31e-06       TGGCAGCCAGGGGGAGGTG
# rewirte tbl with this corrected finding

refimo <- function(tbl, shoulder=1)
{
   for(r in 1:nrow(tbl)){
     chrom <- tbl$chrom[r]
     start <- tbl$start[r]
     end <-  tbl$endpos[r]
     tbl.fimo <- getFimoHits(chrom, start-shoulder, end+shoulder)
     if(nrow(tbl.fimo) == 0) next;
     replacement <-tbl.fimo[order(tbl.fimo$pval),][1,,drop=TRUE]
     if(replacement$motifname ==  tbl$motif[r]){
        printf("changing %s  %s", tbl$loc[r], tbl$motif[r])
        tbl$start[r] <- replacement$start
        tbl$endpos[r] <- replacement$endpos
        tbl$loc[r] <- sprintf("%s:%d-%d", tbl$chrom[r], tbl$start[r], tbl$endpos[r]);
        tbl$length[r] <- 1 + tbl$endpos[r] - tbl$start[r]
        printf("new loc: %s", tbl$loc[r]);
        } # if matched motif, possibly changed start:end
   } # for r
   
   tbl   

} # refimo
#------------------------------------------------------------------------------------------------------------------------
test.refimo <- function()
{
   printf("--- test.refimo")
   if(!exists("tbl.piq3"))
      load("tbl.piq3.RData", envir=.GlobalEnv)
   x <- refimo(tbl.piq3[35,,drop=FALSE])
   checkEquals(x$loc, "chr19:44904896-44904914")
   checkEquals(x$length, 19)
   tbl.changed <- refimo(tbl.piq3)
   checkEquals(dim(tbl.changed), dim(tbl.piq3))
   checkEquals(length(setdiff(tbl.changed$loc, tbl.changed$piq3)), 39)

} # test.refimo
#------------------------------------------------------------------------------------------------------------------------
# find all overlaps between the 151 base pair chipseq regions, and short motif-based fimo regions
# then expand the tbl.cs by joining it with all the tbl.fimo regions which overlap with each chipseq region
# then whittle those down, keeping only those fimo/chipseq rows where the fimo motif is in fact
# associated with the transcript factor reported in tbl.cs
addFimoRegions <- function(tbl.cs, tbl.fimo)
{
   #printf("--- entering addFimoRegions, %d cs, %d fimo", nrow(tbl.cs), nrow(tbl.fimo))
    
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

   #printf("tbl.combined: %d %d", nrow(tbl.combined), ncol(tbl.combined))
   
      # discard the rows where the fimo motif is not associated with the tf
   tbl.out <- tbl.combined
   #tbl.out <- subset(tbl.combined, bindingSite==TRUE)

      # add a loc & length field
   tbl.out$loc <- with(tbl.out, sprintf("%s:%d-%d", chrom, start, endpos))
   tbl.out$length <- with(tbl.out, 1 + endpos - start)

   # printf("tbl.out: %d %d", nrow(tbl.out), ncol(tbl.out))

   invisible(tbl.out)

} # addFimoRegions
#------------------------------------------------------------------------------------------------------------------------
test.addFimoRegions <- function()
{
   printf("--- test.addFimoRegions")
   chrom <- "chr19"
   start <- 44904800
   end <-   44905000
   tbl.chipseq <- getHits(db.chipseq, chrom, start, end)
   tbl.fimo <- getFimoHits(chrom, start, end)
   tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")
   tbl.expanded <- addFimoRegions(tbl.chipseq, tbl.fimo)
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
old.toFeatureTable <- function(tbl.hits, methodName)
{
   motifNames <- paste(methodName, sort(unique(tbl.hits$name)), sep="_")
   column.names <- c("uLoc", motifNames)
   uLocs <- unique(tbl.hits$loc)
   tbl <- data.frame(matrix(data=0, nrow=length(uLocs), ncol=length(column.names)), stringsAsFactors=FALSE)
   colnames(tbl) <- column.names
   tbl$uLoc <- uLocs

   for(r in 1:nrow(tbl.hits)){
      row <- tbl.hits[r, "loc"]
      col <- sprintf("%s_%s", methodName, tbl.hits[r, "name"])
      #printf("[%s, %s]", row, col)
      tbl[grep(row, tbl$uLoc), col] <- tbl[grep(row, tbl$uLoc), col] + 1
      }

   invisible(tbl)

} # old.toFeatureTable
#------------------------------------------------------------------------------------------------------------------------
toFeatureTable.v1 <- function(tbl.std)
{
   browser()
   #motifNames <- paste(methodName, sort(unique(tbl.hits$name)), sep="_")
   column.names <- c("uLoc", motifNames)
   uLocs <- unique(tbl.hits$loc)
   tbl <- data.frame(matrix(data=0, nrow=length(uLocs), ncol=length(column.names)), stringsAsFactors=FALSE)
   colnames(tbl) <- column.names
   tbl$uLoc <- uLocs

   for(r in 1:nrow(tbl.hits)){
      row <- tbl.hits[r, "loc"]
      col <- sprintf("%s_%s", methodName, tbl.hits[r, "name"])
      #printf("[%s, %s]", row, col)
      tbl[grep(row, tbl$uLoc), col] <- tbl[grep(row, tbl$uLoc), col] + 1
      }

   invisible(tbl)

} # toFeatureTable.v1
#------------------------------------------------------------------------------------------------------------------------
test.toFeatureTable.v1 <- function()
{
   printf("--- test.toFeatureTable.v1")
   load("hint.normalized.5rows.10hits.RData")
   ft <- toFeatureTable.v1(tbl.h)

} # test.toFeatureTable.v1
#------------------------------------------------------------------------------------------------------------------------
chipseqToFeatureTable <- function(tbl.hits, methodName)
{
   motifNames <- paste(methodName, sort(unique(tbl.hits$name)), sep="_")
   column.names <- c("uLoc", motifNames)
   uLocs <- unique(tbl.hits$loc)
   tbl <- data.frame(matrix(data=0, nrow=length(uLocs), ncol=length(column.names)), stringsAsFactors=FALSE)
   colnames(tbl) <- column.names
   tbl$uLoc <- uLocs


   for(r in 1:nrow(tbl.hits)){
      row <- tbl.hits[r, "loc"]
      col <- sprintf("%s_%s", methodName, tbl.hits[r, "name"])
      #printf("[%s, %s]", row, col)
      tbl[grep(row, tbl$uLoc), col] <- tbl[grep(row, tbl$uLoc), col] + 1
      }

   invisible(tbl)

} # chipseqToFeatureTable
#------------------------------------------------------------------------------------------------------------------------
test.chipseqToFeatureTable <- function()
{
   printf("--- test.chipseqToFeatureTable")
   chrom <- "chr19"
   start <- 44906493
   end <-  44906663
   tbl.chipseq <- getHits(db.chipseq, chrom, start, end)
   tbl.fimo <- getFimoHits(chrom, start, end)
   tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")
   tbl.csWithFimo <- addFimoRegions(tbl.chipseq, tbl.fimo)

   tbl.exact.features  <- chipseqToFeatureTable(subset(tbl.csWithFimo, bindingSite==TRUE), "chipseqFimoTfMatch")
   checkEquals(dim(tbl.exact.features), c(2, 2))
   checkEquals(colnames(tbl.exact.features), c("uLoc", "chipseqFimoTfMatch_PBX3"))
   tbl.all.features  <- chipseqToFeatureTable(tbl.csWithFimo, "chipseqFimoRegionMatch")
   checkEquals(colnames(tbl.all.features), c("uLoc", "chipseqFimoRegionMatch_PBX3"))
   checkEquals(sum(tbl.all.features[, -1]), nrow(tbl.csWithFimo))

}  # test.chipseqToFeatureTable
#------------------------------------------------------------------------------------------------------------------------
test.old.toFeatureTable <- function(shoulder=100)
{
   tbl.apoe <- getHits(db.hint, apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
     # chr19:44,904,793-44,904,958 includes two 151-bp chipseq hits about 800bp upstream of apoe tss
     # chr19:44,906,493-44,906,663 includes one 151-bp chipssq hit (PBX3) just downstream of the tss
   tbl.hhits <- getHits(db.hint, "chr19", 44906493, 44906663)
   tbl.whits <- getHits(db.wellington, "chr19", 44906493, 44906663)
   tbl.hft <- old.toFeatureTable(tbl.hhits, "hint")
   tbl.wft <- old.toFeatureTable(tbl.whits, "wellington")
   #checkEquals(nrow(tbl.hits), sum(tbl[, -1]))

      # for every row, identify each column with a 1, make sure tbl.expanded[

   for(r in 1:nrow(tbl)){
      hits <- which(tbl[r,] == 1)
      this.loc <- tbl$uLoc[r]
      if(length(hits) > 0){
         column.names <- sub("HINT_", "", colnames(tbl)[hits])
         #if(length(column.names) > 10) browser()
         #printf("loc: %s   column.names: %s", this.loc, paste(column.names, collapse=","))
         #printf("checked out? %s (%d,%d)", 
         #       nrow(subset(tbl.hits, name %in% column.names & loc==this.loc)) == length(hits),
         #       nrow(subset(tbl.hits, name %in% column.names & loc==this.loc)), length(hits))
         checkEquals(nrow(subset(tbl.hits, name %in% column.names & loc==this.loc)), length(hits))
         } # if hits
      } # for r


} # test.old.toFeatureTable
#------------------------------------------------------------------------------------------------------------------------
oldEnsembl <- function(chrom, start, end)
{
   #browser()
   tbl.hits.cs <- getHits(db.chipseq, chrom, start, end)
   tbl.fimo.cs <- getFimoHits(chrom, start, end)
   tbl.fimo.cs$chrom <- paste("chr", tbl.fimo.cs$chrom, sep="")

   tbl.expanded.cs <- addFimoRegions(tbl.hits.cs, tbl.fimo.cs)
   tbl.feature.cs <- old.toFeatureTable(tbl.expanded.cs)

   tbl.hits.hint <- getHits(db.hint, chrom, start, end)
   tbl.feature.hint <- hintToFeatureTable(tbl.hits.hint)
   #printf ("uLoc overlap: %d", length(intersect(tbl.feature.hint$uLoc, tbl.feature.cs$uLoc)))
   x <- tbl.feature.cs$uLoc; y <- tbl.feature.hint$uLoc;
   printf("chipseq regions: %d", length(x))
   printf("   HINT regions: %d", length(y))
   browser()
   x <- 99
    
} # oldEnsembl
#------------------------------------------------------------------------------------------------------------------------
test.oldEnsembl <- function()
{
   printf("--- test.ensembl")
   shoulder <- 1000
   tbl <- ensembl(apoe$chrom, apoe$start - shoulder, apoe$start + shoulder)
   
} # test.ensembl
#------------------------------------------------------------------------------------------------------------------------
ensemble <- function(chrom, start, end)
{
   tbl.w <- createWellingtonTable(chrom, start, end)
   tbl.h <- createHintTable(chrom, start, end)
   tbl.merged <- merge(tbl.w, tbl.h, by=c("loc"), all=TRUE)
   browser()
   
   tbl.chipseq <- getHits(db.chipseq, chrom, start, end)
   stopifnot(nrow(tbl.chipseq) > 0)
   tbl.fimo <- getFimoHits(chrom, start, end)
   tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")
   tbl.csWithFimo <- addFimoRegions(tbl.chipseq, tbl.fimo)
   tbl.csWithTrueFimo <- subset(tbl.csWithFimo, bindingSite==TRUE)
   
   tbl.csft <- tbl.csWithTrueFimo[, c("loc", "motifname", "name", "score1")]
   colnames(tbl.csft) <- c("loc", "csmotif", "csTF", "csscore")
   tbl.wfc <- merge(tbl.merged, tbl.csft, by="loc", all=TRUE)

   classInfo <- sapply(tbl.wfc, class)
   for(colname in names(classInfo)){
      if(classInfo[[colname]] == "factor"){
          printf("correctoring character class in column '%s'", colname)
          tbl.wfc[[colname]] <- as.character(tbl.wfc[[colname]])
          } # if factor
      } # for colname

   invisible(tbl.wfc)

} # ensemble
#------------------------------------------------------------------------------------------------------------------------
test.ensemble <- function()
{
   printf("--- test.ensemble")

   chrom <- "chr19"
     # 229 bases
   start <- 44903672
   end   <- 44903900

   ft <- ensemble(chrom, start, end)
   checkEquals(dim(ft), c(33, 22))
   tfHits <- ft[which(!is.na(ft$csTF)), "csTF"]
   checkEquals(length(tfHits), 1)
   checkEquals(tfHits, "CTCF")
    
     # test in one 12 base region, very rich in motifs, includes all most combinations of present & missing motifs
   chrom <- "chr19"
   target.loc <- "chr19:45423911-45423920"
   start <- 45423910
   end   <- 45423921
   tbl.h <- getHits(db.hint, chrom, start, end)
   tbl.w <- getHits(db.wellington, chrom, start, end)
   ft.x  <- ensemble(chrom, start-1000, end + 1000) # be sure to get a chipseq hit.  5400 hits
     # eliminate the off target rows (present only because the footprint is wider than our target.loc
   ft.xs <- ft.x # subset(ft.x, loc==target.loc)  # 5401
   ft.xs.1 <- subset(ft.xs, (motif.w == csmotif & is.na(motif.h)) |
                            (motif.h == csmotif & is.na(motif.w)) |
                            (motif.w == motif.h & motif.w == csmotif) |
                            (is.na(csTF)))  # 27 22

     # tbl <- ft.xs.1
     # save(tbl, file="shortRichEnsemblTestResult.RData")

     # make sure there are no factors!
   checkEquals(sort(unique(unlist(lapply(ft.xs.1, class), use.names=FALSE))), c("character",  "integer",  "numeric"))

     # make sure that we have rows of all kinds:  hint, wellington and chipseq and shared hits
   checkEquals(nrow(ft.xs.1), 319)

   checkEquals(nrow(subset(ft.xs.1, !is.na(motif.h)  &  is.na(motif.w)  &  is.na(csmotif))), 81)
   checkEquals(nrow(subset(ft.xs.1, !is.na(motif.h)  &  is.na(motif.w)  & !is.na(csmotif))), 20)
   checkEquals(nrow(subset(ft.xs.1, !is.na(motif.h)  & !is.na(motif.w)  &  is.na(csmotif))), 94)
   checkEquals(nrow(subset(ft.xs.1, !is.na(motif.h)  & !is.na(motif.w)  & !is.na(csmotif))), 105)

   checkEquals(nrow(subset(ft.xs.1, is.na(motif.h)  &  is.na(motif.w)  &  is.na(csmotif))),   0)
   checkEquals(nrow(subset(ft.xs.1, is.na(motif.h)  &  is.na(motif.w)  & !is.na(csmotif))),   0)
   checkEquals(nrow(subset(ft.xs.1, is.na(motif.h)  & !is.na(motif.w)  &  is.na(csmotif))),  16)
   checkEquals(nrow(subset(ft.xs.1, is.na(motif.h)  & !is.na(motif.w)  & !is.na(csmotif))),   3)
   checkEquals(81 + 20 + 94 + 105 + 0 + 0 + 16 + 3, 319)

     #--------- test a rich region: 89 hint hits of 6 motifs
     # MA0058.3 MA0104.3 MA0617.1 MA0622.1 MA0626.1 MA0825.1 
     #       14       15       16       16       14       14 
     #  65 wellington 
     #  MA0058.3 MA0104.3 MA0617.1 MA0622.1 MA0626.1 MA0825.1 
     #        10       11       12       12       10       10 




} # test.ensemble
#------------------------------------------------------------------------------------------------------------------------
locStringToBedTable <- function(locStrings)
{
   tokensList <- strsplit(locStrings, ":")
   chroms     <- unlist(lapply(tokensList, "[", 1))
   posPairStrings  <- unlist(lapply(tokensList, "[", 2))

   posPairs <- lapply(strsplit(posPairStrings, "-"), as.integer)
   tbl.startEnd <- data.frame(matrix(as.integer(unlist(posPairs)), ncol=2, byrow=TRUE))
   
   tbl <- cbind(chrom=chroms, tbl.startEnd, stringsAsFactors=FALSE)   
   colnames(tbl) <- c("chrom", "start", "end")
   tbl[order(tbl$chrom, tbl$start),]
   
}  # locStringToBedTable
#------------------------------------------------------------------------------------------------------------------------
test.locStringToBedTable <- function()
{
   printf("--- test.locStringToBedTable")

   s <- c("chr19:44906708-44906728", "chr19:44906711-44906731", "chr19:44906550-44906559")
   tbl <- locStringToBedTable(s)
   checkEquals(colnames(tbl), c("chrom", "start", "end"))
   checkEquals(as.list(tbl[1,]), list(chrom="chr19", start=44906550, end=44906559))
   checkEquals(as.list(tbl[3,]), list(chrom="chr19", start=44906711, end=44906731))
   checkEquals(as.list(tbl[2,]), list(chrom="chr19", start=44906708, end=44906728))

}  # test.locStringToBedTable
#------------------------------------------------------------------------------------------------------------------------
# fivenum(scores.hint)                 2          8         12         20     17836
# fivenum(scores.wellington)  -1000.0000   -36.8523   -18.9043   -12.9359   -10.0000
# make missing values for score1.h.median, score1.h.best 0
# make missing values for score1.w.median 0
# collapse motif.w and motif.h into one 'motif' column
cleanFeatureTable <- function(tbl)
{
   x <- tbl

   x$score1.w.median[is.na(x$score1.w.median)] <- 0
   x$score1.w.best[is.na(x$score1.w.best)] <- 0
   x$score1.h.median[is.na(x$score1.h.median)] <- 0
   x$score1.h.best[is.na(x$score1.h.best)] <- 0

   x$score2.w.median[is.na(x$score2.w.median)] <- -99
   x$score2.w.best[is.na(x$score2.w.best)] <- -99
   x$score2.h.median[is.na(x$score2.h.median)] <- -99
   x$score2.h.best[is.na(x$score2.h.best)] <- -99

   x$score3.w.median[is.na(x$score3.w.median)] <- 1.0
   x$score3.w.best[is.na(x$score3.w.best)] <- 1.0
   x$score3.h.median[is.na(x$score3.h.median)] <- 1.0
   x$score3.h.best[is.na(x$score3.h.best)] <- 1.0

   x$samplecount.w[is.na(x$samplecount.w)] <- 0
   x$length.w[is.na(x$length.w)] <- 0
   
   x$samplecount.h[is.na(x$samplecount.h)] <- 0
   x$length.h[is.na(x$length.h)] <- 0
   
   x$csmotif[is.na(x$csmotif)] <- "NA"
   x$csTF[is.na(x$csTF)] <- "NA"

      # collapse motifs from the two methods into one column
      # some will be identical, NAs will be complementary

   motif <- as.character(x$motif.w)
   nas <- which(is.na(motif)) 
   motif[nas] <- as.character(x$motif.h[nas])
   x$motif <- motif

   motif.columns.to.delete <- grep("motif.", colnames(x), fixed=TRUE)
   if(length(motif.columns.to.delete) > 0)
      x <- x[, -motif.columns.to.delete]

   motiflength <- x$length.w
   zeros <- which(motiflength == 0)
   motiflength[zeros] <- x$length.h[zeros]
   x$motiflength <- motiflength

   length.columns.to.delete <- grep("length.", colnames(x), fixed=TRUE)
   if(length(length.columns.to.delete) > 0)
      x <- x[, -length.columns.to.delete]

   x$csscore[which(is.na(x$csscore))] <- 0
   
   invisible(unique(x))
   
} # cleanFeatureTable
#------------------------------------------------------------------------------------------------------------------------
test.cleanFeatureTable <- function()
{
   printf("--- test.cleanFeatureTable")

   load ("shortRichEnsemblTestResult.RData")
   checkEquals(dim(tbl), c(319,22))

   ft <- cleanFeatureTable(tbl)

   checkTrue("motif" %in% colnames(ft))
   checkTrue(!"motif.h" %in% colnames(ft))
   checkTrue(!"motif.w" %in% colnames(ft))
   checkEquals(nrow(tbl), nrow(ft))

   checkTrue(all(ft$score1.w.median <= 0))
   checkTrue(all(ft$score1.w.best <= 0))
    
   checkTrue(all(ft$score2.w.median >= -99))
   checkTrue(all(ft$score2.w.best >= -99))

   checkTrue(all(ft$score3.w.median >= 0))
   checkTrue(all(ft$score3.w.median <= 1))

   checkTrue(all(ft$score3.w.best >= 0))
   checkTrue(all(ft$score3.w.best <= 1))

   checkTrue(all(ft$score1.h.median >= 0))
   checkTrue(all(ft$score1.h.best >= 0))
    
   checkTrue(all(ft$score2.h.median >= -99))
   checkTrue(all(ft$score2.h.best >= -99))

   checkTrue(all(ft$score3.h.median >= 0))
   checkTrue(all(ft$score3.h.median <= 1))

   checkTrue(all(ft$score3.h.best >= 0))
   checkTrue(all(ft$score3.h.best <= 1))

   checkTrue(all(ft$samplecount.w >= 0))
   checkTrue(all(ft$samplecount.h >= 0))

   checkTrue(all(ft$length.w >= 0))
   checkTrue(all(ft$length.h >= 0))

} # test.cleanFeatureTable
#------------------------------------------------------------------------------------------------------------------------
run <- function()
{
   chrom <- "chr19"
   start <- 42000000
   end   <- 47000000

   ft.0 <- ensemble(chrom, start, end)
   ft <- cleanFeatureTable(ft.0)
   filename <- sprintf("ft.%s.RData", format (Sys.time(), "%a.%b.%d.%Y-%H:%M:%S"))
   save(ft, file=filename)
   
} # run
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
    run()
