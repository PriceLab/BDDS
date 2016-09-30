#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
runTestUtils <- function()
{
   test.createWellingtonTable()
   test.createHintTable()
   test.createHintTable_ignoreStrand()
   test.createWellingtonTable_ignoreStrand()
   
} # runTestUtils
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
   tbl.out <- merge(tbl.regions, tbl.hits, on="loc")
   tbl.out

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
createWellingtonTable <- function(chrom, start, end, motif=NA, sampleID=NA, collapseOnStrand=FALSE)
{
   tbl.hitsw <- getHits(db.wellington, chrom, start, end) # 6 x 17

   if(!is.na(motif))
      tbl.hitsw <- subset(tbl.hitsw, name==motif)
   if(!is.na(sampleID))
      tbl.hitsw <- subset(tbl.hitsw, sample_id==sampleID)

   if(collapseOnStrand){
       strand.duplications <- which(duplicated(tbl.hitsw[, c("loc", "name", "sample_id")]))
       if(length(strand.duplications) > 0){
           printf("eliminating %d double-stranded hits", length(strand.duplications))
           tbl.hitsw <- tbl.hitsw[-strand.duplications,]
           }  # some hits for loc & motif on both strands
       } # collapseOnStrand


   printf("found %d wellington hits in %d bases", nrow(tbl.hitsw), 1 + end - start)
   #displayBedTable(igv, tbl.hitsw[, c("chrom", "start", "endpos", "name", "score2")], "wellington")
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
       #printf("count of duplicated samples: %d, %s, %s", length(which(duplicated(tbl.sub$sample_id))),
       #       this.loc, this.motif)
       #browser()
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

   tbl.w <- createWellingtonTable(chrom, start, end, motif="GTF2A1,2.p2")
   tbl.w <- createWellingtonTable(chrom, start, end, motif="MA0813.1", sampleID="ENCSR000EJE")
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
   

}  # test.createWellingtonTable
#------------------------------------------------------------------------------------------------------------------------
createHintTable <- function(chrom, start, end, motif=NA, sample=NA, loc=NA, collapseOnStrand=FALSE)
{
   tbl.hitsh <- getHits(db.hint, chrom, start, end) # 6 x 17

   printf("found %d hint hits in %d bases", nrow(tbl.hitsh), 1 + end - start)

   if(collapseOnStrand){
       strand.duplications <- which(duplicated(tbl.hitsh[, c("loc", "name", "sample_id")]))
       if(length(strand.duplications) > 0){
           printf("eliminating %d double-stranded hits", length(strand.duplications))
           tbl.hitsh <- tbl.hitsh[-strand.duplications,]
           }  # some hits for loc & motif on both strands
       } # collapseOnStrand

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
   start <- 45423537
   end <- 45423550
   motif <- "MA0512.2"
   loc <- "chr19:45423537-45423550"
   tbl.h <- createHintTable(chrom, start, end, motif=motif, loc=loc)


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

   x <- createWellingtonTable(chrom, start, end, motif=test.motif, sampleID=NA, collapseOnStrand=FALSE)
   checkEquals(x$samplecount.w, 26)

   y <- createWellingtonTable(chrom, start, end, motif=test.motif, sampleID=NA, collapseOnStrand=TRUE)
   checkEquals(y$samplecount.w, 13)
   
} # test.createWellingtonTable_ignoreStrand
#------------------------------------------------------------------------------------------------------------------------
