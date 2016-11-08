library(RPostgreSQL)
library(GenomicRanges)
library(RUnit)
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


#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
utils.runTests <- function()
{
   test.locStringToBedTable()
   test.getFimoHits()
   test.getHits()
   test.createHintTable()
   test.createHintTable_ignoreStrand()
   test.createWellingtonTable()
   test.createWellingtonTable_ignoreStrand()
   test.createPiqTable()
   
} # runTestUtils
#------------------------------------------------------------------------------------------------------------------------
getHits <- function(db, chrom, start, end, locs=NA, motifs=NA)
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
   tbl.out <- merge(tbl.regions, tbl.hits, on="loc")

   if(!any(is.na(motifs)))
       tbl.out <- subset(tbl.out, name %in% motifs)

   if(!any(is.na(locs)))
       tbl.out <- subset(tbl.out, loc %in% locs)

   tbl.out

} # getHits
#------------------------------------------------------------------------------------------------------------------------
test.getHits <- function()
{
   printf("--- test.getHits")

   chrom <- "chr19"
   start <- 44903772
   end   <- 44903790

   x <- getHits(db.hint, chrom, start, end)
   checkEquals(dim(x), c(10, 17))

   target.locs <- c("chr19:44903775-44903786")
   x2 <- getHits(db.hint, chrom, start, end, locs=target.locs)
   checkEquals(dim(x2), c(2, 17))
   checkEquals(unique(x2$loc), target.locs)

   target.motifs <- c("MA0696.1", "sci09.v1_Plagl1_0972")
   x3 <- getHits(db.hint, chrom, start, end, motifs=target.motifs)
   checkEquals(dim(x3), c(4, 17))
   checkEquals(sort(unique(x3$name)), sort(target.motifs))

   target.locs <- c("chr19:44890339-44890348", "chr19:44890339-44890348")
   target.motifs <- c("MA0612.1", "MA0644.1")
   x4 <- getHits(db.hint, chrom, start-20000, end+20000, locs=target.locs, motifs=target.motifs)
   x5 <- getHits(db.hint, chrom, start-20000, end+20000, locs=target.locs)
   x6 <- getHits(db.hint, chrom, start-20000, end+20000, motifs=target.motifs)

   checkEquals(dim(x4), c(8,  17))
   checkEquals(dim(x5), c(76, 17))
   checkEquals(dim(x6), c(16, 17))
    
} # test.getHits
#------------------------------------------------------------------------------------------------------------------------
getFimoHits <- function(chrom, start, end, exact=FALSE)
{
   chrom <- sub("^chr", "", chrom)

   if(exact)
       query <- sprintf("select * from fimo_hg38 where chrom='%s' and start = %d and endpos = %d", chrom, start, end)
   else
      query <- sprintf("select * from fimo_hg38 where chrom='%s' and start >= %d and endpos <= %d", chrom, start, end)

   dbGetQuery(db.trena, query)

} # getFimoHits
#------------------------------------------------------------------------------------------------------------------------
test.getFimoHits <- function()
{
   printf("--- test.getFimoHits")

   chrom <- "chr19"
   start <- 44903772
   end   <- 44903790

   x <- getFimoHits(chrom, start, end)
   checkEquals(dim(x), c(15, 9))
   expected <- c("motifname", "chrom", "start", "endpos", "strand", "motifscore", "pval", "sequence")
   checkTrue(all(expected %in% colnames(x)))
    
} # test.getFimoHits
#------------------------------------------------------------------------------------------------------------------------
# refactor this into a single createMethodTable, for both hint and wellington?
createWellingtonTable <- function(chrom, start, end, motifs=NA, locs=NA, collapseOnStrand=FALSE)
{
   tbl.hitsw <- getHits(db.wellington, chrom, start, end, motifs=motifs, locs=locs) # 6 x 17
   printf("found %d wellington hits in %d bases", nrow(tbl.hitsw), 1 + end - start)
   if(nrow(tbl.hitsw) == 0)
       return(data.frame())

   if(collapseOnStrand){
       with(tbl.hitsw, tbl.hitsw <- tbl.hitsw[order(loc, name, sample_id,  score3, decreasing=FALSE),])
       strand.duplications <- which(duplicated(tbl.hitsw[, c("loc", "name", "sample_id")]))
       if(length(strand.duplications) > 0){
           printf("eliminating %d double-stranded hits", length(strand.duplications))
           tbl.hitsw <- tbl.hitsw[-strand.duplications,]
           }  # some hits for loc & motif on both strands
       } # collapseOnStrand

   tbl.std <- tbl.hitsw[, c("loc",  "name", "length", "score1", "score2", "score3", "sample_id")]
   tbl.collapsed <- as.data.frame(table(tbl.std$loc, tbl.std$name))
   tbl.collapsed <- subset(tbl.collapsed, Freq != 0)
   colnames(tbl.collapsed) <- c("loc", "motif.w", "sample.count.w")

      # preallocate and zero fill the summary columns for each row

   distinct.hits.count <- nrow(tbl.collapsed)  # unique loc/motif combinations
   length <- rep(0, distinct.hits.count)
   score1.median <- rep(0.0, distinct.hits.count)
   score1.best <-   rep(0.0, distinct.hits.count)
   score2.median <- rep(0.0, distinct.hits.count)
   score2.best  <-   rep(0.0, distinct.hits.count)
   score3.median <- rep(0.0, distinct.hits.count)
   score3.best <-   rep(0.0, distinct.hits.count)
   length <-        rep(0,   distinct.hits.count)

   for(r in 1:nrow(tbl.collapsed)){
       this.loc <- tbl.collapsed$loc[r]
       this.motif <- tbl.collapsed$motif.w[r]
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
   colnames(tbl.out) <- c("loc", "motif.w", "samplecount.w", "length.w", "score1.w.median", "score1.w.best",
                          "score2.w.median",  "score2.w.best", "score3.w.median", "score3.w.best")
   tbl.out$motif.w <- as.character(tbl.out$motif.w)
   tbl.out$loc <- as.character(tbl.out$loc)
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

   #tbl.w <- createWellingtonTable(chrom, start, end, motif="GTF2A1,2.p2")
   #tbl.w <- createWellingtonTable(chrom, start, end, motif="MA0813.1", sampleID="ENCSR000EJE")
   tbl.w <- createWellingtonTable(chrom, start, end)
   checkTrue(! "factor" %in% as.character(lapply(tbl.w, class)))
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

   motif <- "MA0760.1"
   tbl.wx2 <- createWellingtonTable("chr19", 45423910, 45423921, motifs=motif)
   checkEquals(nrow(tbl.wx2), 1)


   start <- 44903772
   end   <- 44903790
   locs <- c("chr19:44903773-44903784", "chr19:44903773-44903784", "chr19:44903773-44903785")
   tbl.w <- createWellingtonTable(chrom, start, end); # , locs=locs)
   checkEquals(nrow(tbl.w), 12)

   tbl.w2 <- createWellingtonTable(chrom, start, end, locs=locs)
   checkEquals(nrow(tbl.w2), 5)
   checkTrue(all(tbl.w2$loc %in% locs))

   motifs <- c("MA0810.1", "MA0815.1")
   tbl.w3 <- createWellingtonTable(chrom, start, end, locs=locs, motifs=motifs)
   checkEquals(nrow(tbl.w3), 2)
   checkTrue(all(tbl.w3$loc %in% locs))
   checkTrue(all(tbl.w3$motif.w %in% motifs))
   
}  # test.createWellingtonTable
#------------------------------------------------------------------------------------------------------------------------
createHintTable <- function(chrom, start, end, motifs=NA, locs=NA, collapseOnStrand=FALSE)
{
   tbl.hitsh <- getHits(db.hint, chrom, start, end, motifs=motifs, locs=locs)
   printf("found %d hint hits in %d bases", nrow(tbl.hitsh), 1 + end - start)

   if(nrow(tbl.hitsh) == 0)
      return(data.frame())

   if(collapseOnStrand){
         # resort tbl.hitsh so that all rows with equal loc/name/sample_id are sorted
         # in descending order by score3, the fimo pval.  then, when any duplicates are eliminated
         # they are the lower scoring rows
       with(tbl.hitsh, tbl.hitsh <- tbl.hitsh[order(loc, name, sample_id,  score3, decreasing=FALSE),])
       strand.duplications <- which(duplicated(tbl.hitsh[, c("loc", "name", "sample_id")]))
       if(length(strand.duplications) > 0){
           printf("eliminating %d double-stranded hits", length(strand.duplications))
           tbl.hitsh <- tbl.hitsh[-strand.duplications,]
           }  # some hits for loc & motif on both strands
       } # collapseOnStrand

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
   tbl.out$motif.h <- as.character(tbl.out$motif.h)
   tbl.out$loc <- as.character(tbl.out$loc)

   tbl.out

}  # createHintTable
#------------------------------------------------------------------------------------------------------------------------
test.createHintTable <- function()
{
   printf("--- test.createHintTable")

      # extract a very small table with confusing scores
   chrom <- "chr19"
   start <- 45423500
   end   <- 45423600
   tbl <- createHintTable(chrom, start, end)
   checkEquals(nrow(tbl), 49)

   motif <- "MA0599.1"
   tbl <- createHintTable(chrom, start, end, motifs=motif)
   checkEquals(nrow(tbl), 2)
   checkEquals(tbl$motif.h, c(motif, motif))
   checkEquals(length(unique(tbl$loc)), 2)

   loc <- "chr19:45423532-45423541"
   tbl <- createHintTable(chrom, start, end, locs=loc)
   checkEquals(nrow(tbl), 2)
   checkEquals(tbl$loc, c(loc, loc))
   checkEquals(length(unique(tbl$motif.h)), 2)

   start <- 44903772
   end   <- 44903790
   tbl.h <- createHintTable(chrom, start, end)
   checkTrue(! "factor" %in% as.character(lapply(tbl.h, class)))

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
test.createHintTable_ignoreStrand <- function()
{
    printf("--- test.createHintTable_ignoreStrand")
    x <- createHintTable("chr19", 920283, 920294)

      #                   loc  motif.h samplecount.h length.h score1.h.median score1.h.best
      # 1 chr19:920283-920294 MA0524.2            11       12              48           424
      # 2 chr19:920283-920294 MA0810.1            11       12              48           424
      # 3 chr19:920283-920294 MA0811.1            22       12              48           424

    checkEquals(x$samplecount.h, c(11, 11, 22))

    x2 <- createHintTable("chr19", 920283, 920294, collapseOnStrand=TRUE)
    checkEquals(x2$samplecount.h, c(11, 11, 11))

} # test.createHintTable_ignoreStrand
#------------------------------------------------------------------------------------------------------------------------
test.createWellingtonTable_ignoreStrand <- function()
{
   printf("--- test.createWellingtonTable_ignoreStrand, deferred")
   loc <- "chr19:572471-572480"  #  MA0632.1   (empirically discovered)
   chrom <- "chr19"
   start <- 572471
   end <- 572480
   test.motif <- "MA0632.1"

   x <- createWellingtonTable(chrom, start, end, motifs=test.motif, collapseOnStrand=FALSE)
   checkEquals(x$samplecount.w, 26)

   y <- createWellingtonTable(chrom, start, end, motif=test.motif, collapseOnStrand=TRUE)
   checkEquals(y$samplecount.w, 13)
   
} # test.createWellingtonTable_ignoreStrand
#------------------------------------------------------------------------------------------------------------------------
createPiqTable <- function(chrom, start, end, motifs=NA, locs=NA, collapseOnStrand=FALSE)
{
   tbl.hits <- getHits(db.piq, chrom, start, end, motifs=motifs, locs=locs)
   printf("found %d piq hits in %d bases", nrow(tbl.hits), 1 + end - start)

   if(nrow(tbl.hits) == 0)
      return(data.frame())

   if(collapseOnStrand){
         # resort tbl.hits so that all rows with equal loc/name/sample_id are sorted
         # in descending order by score3, the fimo pval.  then, when any duplicates are eliminated
         # they are the lower scoring rows
       with(tbl.hits, tbl.hits <- tbl.hits[order(loc, name, sample_id,  score3, decreasing=FALSE),])
       strand.duplications <- which(duplicated(tbl.hits[, c("loc", "name", "sample_id")]))
       if(length(strand.duplications) > 0){
           printf("eliminating %d double-stranded hits", length(strand.duplications))
           tbl.hits <- tbl.hits[-strand.duplications,]
           }  # some hits for loc & motif on both strands
       } # collapseOnStrand

   tbl.std <- tbl.hits[, c("loc",  "name", "length", "sample_id", "score1", "score2", "score3", "score4")]
   tbl.collapsed <- as.data.frame(table(tbl.std$loc, tbl.std$name))
   tbl.collapsed <- subset(tbl.collapsed, Freq != 0)
   colnames(tbl.collapsed) <- c("loc", "motif.p", "sample.count.p")


      # preallocate and zero fill the summary columns for each row

   distinct.hits.count <- nrow(tbl.collapsed)  # unique loc/motif combinations
   length <- rep(0, distinct.hits.count)

   score1.median <- rep(0.0, distinct.hits.count)
   score1.best <-   rep(0.0, distinct.hits.count)

   score2.median <- rep(0.0, distinct.hits.count)
   score2.best <-   rep(0.0, distinct.hits.count)

   score3.median <- rep(0.0, distinct.hits.count)
   score3.best <-   rep(0.0, distinct.hits.count)

   score4.median <- rep(0.0, distinct.hits.count)
   score4.best <-   rep(0.0, distinct.hits.count)

   length <-        rep(0,   distinct.hits.count)

    # score1: hint score, ranges from 2 to 10,160, mean of ~123, median 44
    # score2: fimo score, -18.7 to 30.9, mean and median both about 11.5
    # score3: fimo pval, min 5.9e-13, max 1e-4, mean 4.3e-5, median 4.01 e-5

   for(r in 1:nrow(tbl.collapsed)){
       this.loc <- tbl.collapsed$loc[r]
       this.motif <- tbl.collapsed$motif.p[r]
       tbl.sub <- subset(tbl.std, loc==this.loc & name==this.motif)
       score1.median[r] <- median(tbl.sub$score1)
       score1.best[r]   <- max(tbl.sub$score1)
       score2.median[r] <- median(tbl.sub$score2)
       score2.best[r]   <- max(tbl.sub$score2)
       score3.median[r] <- median(tbl.sub$score3)
       score3.best[r]   <- max(tbl.sub$score3)
       score4.median[r] <- median(tbl.sub$score4)
       score4.best[r]   <- max(tbl.sub$score4)
       length[r] <- median(tbl.sub$length) # should all be identical, but this covers all situations
       #printf("  found %d rows for %s %s", nrow(tbl.sub), this.loc, this.motif)
       x <- 99
       } # for r

   tbl.out <- cbind(tbl.collapsed, length, score1.median, score1.best, score2.median, score2.best,
                    score3.median, score3.best, score4.median, score4.best)
   colnames(tbl.out) <- c("loc", "motif.p", "samplecount.p", "length.p", "score1.p.median",  "score1.p.best",
                          "score2.p.median",  "score2.p.best", "score3.p.median", "score3.p.best",
                          "score4.p.median", "score4.p.best")
   tbl.out$motif.p <- as.character(tbl.out$motif.p)
   tbl.out$loc <- as.character(tbl.out$loc)

   tbl.out

}  # createPiqTable
#------------------------------------------------------------------------------------------------------------------------
test.createPiqTable <- function()
{
   printf("--- test.createPiqTable")

      # extract a very small table with confusing scores
   chrom <- "chr19"
   start <- 45423500
   end   <- 45423600
   tbl <- createPiqTable(chrom, start, end)
   checkEquals(dim(tbl), c(19, 12))
   checkEquals(colnames(tbl), c("loc", "motif.p", "samplecount.p", "length.p", "score1.p.median", "score1.p.best",
                                "score2.p.median", "score2.p.best", "score3.p.median", "score3.p.best",
                                "score4.p.median", "score4.p.best"))
   checkEquals(length(unique(tbl$motif.p)), 12)
   checkTrue(all(tbl$samplecount.p == 18))

     # no negative or pvalue scores, so all "best" scores should be >= median scores
     # by direct inspection, what seems plausible is true for scores2-4, that best   
     # is > median
     # todo: recall the nature of these scores, and see if score1.median == score1. best is plausible

   checkTrue(with(tbl, all(score1.p.median <= score1.p.best)))
   checkTrue(with(tbl, all(score2.p.median < score2.p.best)))
   checkTrue(with(tbl, all(score3.p.median < score3.p.best)))
   checkTrue(with(tbl, all(score4.p.median < score4.p.best)))
   
}  # test.createPiqTable
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
#if(!interactive())
#   utils.runTests()
