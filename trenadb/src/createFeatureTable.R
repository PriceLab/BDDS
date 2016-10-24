library(RPostgreSQL)
library(GenomicRanges)
library(RUnit)
#library(igvR)
#------------------------------------------------------------------------------------------------------------------------
source("regionAndHitsSchemas.R")
source("utils.R")
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
   utils.runTests()
   test.locStringToBedTable()
   test.addFimoRegions()
   test.refimoForPiq()

   #test.chipseqToFeatureTable()
   #test.ensemble_hint_wellington_empty_chipseq()
   #test.ensemble_hint_wellington_chipseq()
   test.ensemble.flex.1.loc.only()
   test.ensemble.flex.all.locs()
   
   #test.run50k()
  
} # runTests
#------------------------------------------------------------------------------------------------------------------------
# note: this seems obsolete (21 oct 2016).
# and gets all motifs in that region.
#
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

refimo.old <- function(tbl, shoulder=1)
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

} # refimo.old
#------------------------------------------------------------------------------------------------------------------------
test.refimo.old <- function()
{
   printf("--- test.refimo.old")

   if(!exists("tbl.piq3"))
      load("tbl.piq3.RData", envir=.GlobalEnv)

   x <- refimo.old(tbl.piq3[35,,drop=FALSE])
   checkEquals(x$loc, "chr19:44904896-44904914")
   checkEquals(x$length, 19)
   tbl.changed <- refimo.old(tbl.piq3)
   checkEquals(dim(tbl.changed), dim(tbl.piq3))
   checkEquals(length(setdiff(tbl.changed$loc, tbl.changed$piq3)), 39)

} # test.refimo.old
#------------------------------------------------------------------------------------------------------------------------
# refimoForPiq - for tables acquired from db.piq: takes locs at face value, queries db.fimo
# gets all matching motifs
# to explain by example
# our current piq workflow produces
#    chr19:45423532-45423542 MA0039.2
#    chr19:45423532-45423542 MA0599.1
#
# wellington has
# 1507  chr19:45423532-45423542 KLF16_DBD             7       11        -28.5045      -10.3715    16.8163  1.04e-06
# 17131 chr19:45423532-45423542  MA0741.1             7       11        -28.5045      -10.3715    16.8621  1.04e-06
# 17689 chr19:45423532-45423542  MA0746.1             7       11        -28.5045      -10.3715    15.6748  1.04e-06
#
# and fimo has
# getFimoHits("chr19", 45423530, 45423545)
#   motifname chrom    start   endpos strand motifscore     pval empty       sequence
# 1  MA0039.2    19 45423532 45423541      -    15.8061 5.23e-06           GGGGCGTGGC
# 2  MA0493.1    19 45423531 45423541      +    12.4242 3.07e-05          CGCCACGCCCC
# 3  MA0599.1    19 45423532 45423541      +    15.2653 4.42e-06           GCCACGCCCC
# 4  MA0740.1    19 45423531 45423544      +    18.8276 3.30e-07       CGCCACGCCCCTTT
# 5  MA0741.1    19 45423532 45423542      +    16.8621 1.04e-06          GCCACGCCCCT
# 6  MA0746.1    19 45423532 45423542      +    15.6748 1.04e-06          GCCACGCCCCT
# 7  MA0747.1    19 45423532 45423543      +    16.6579 8.72e-07         GCCACGCCCCTT
# 8 KLF16_DBD    19 45423532 45423542      +    16.8163 1.04e-06          GCCACGCCCCT
# 
# 
# but fimo says
# getFimoHits("chr19", 45423530, 45423545)
#   motifname chrom    start   endpos strand motifscore     pval empty       sequence
# 1  MA0039.2    19 45423532 45423541      -    15.8061 5.23e-06           GGGGCGTGGC
# 2  MA0493.1    19 45423531 45423541      +    12.4242 3.07e-05          CGCCACGCCCC
# 3  MA0599.1    19 45423532 45423541      +    15.2653 4.42e-06           GCCACGCCCC
# 4  MA0740.1    19 45423531 45423544      +    18.8276 3.30e-07       CGCCACGCCCCTTT
# 5  MA0741.1    19 45423532 45423542      +    16.8621 1.04e-06          GCCACGCCCCT
# 6  MA0746.1    19 45423532 45423542      +    15.6748 1.04e-06          GCCACGCCCCT
# 7  MA0747.1    19 45423532 45423543      +    16.6579 8.72e-07         GCCACGCCCCTT
# 8 KLF16_DBD    19 45423532 45423542      +    16.8163 1.04e-06          GCCACGCCCCT
# 
# the function below, refimo.old, changes (shrinks) the region for the supplied motif.
# see "refimo(" for an approach which takes the reported region, takes it at face value
refimoForPiq <- function(tbl)
{
    dups <- which(duplicated(tbl[, c("loc", "samplecount.p")]))
    if(length(dups) > 0)
        tbl <- tbl[-dups,]
    
    locStrings <- unique(tbl$loc)
    tbl.locs <-locStringToBedTable(locStrings)
    tbl.out <- data.frame()

    for(r in 1:nrow(tbl.locs)){
        tbl.fimo <- getFimoHits(tbl.locs$chrom[r], tbl.locs$start[r], tbl.locs$end[r], exact=TRUE)
        tbl.out <- rbind(tbl.out, tbl.fimo)
        }
   tbl.out$loc <- sprintf("chr%s:%d-%d", tbl.out$chrom, tbl.out$start, tbl.out$endpos)

   old.cols.to.keep <- c("loc",   "samplecount.p", "length.p",
                         "score1.p.median", "score1.p.best",
                         "score2.p.median",  "score2.p.best",
                         "score3.p.median", "score3.p.best",
                         "score4.p.median", "score4.p.best")

   new.cols.to.keep <- c("loc", "motifname", "motifscore", "pval")
   merge(tbl[, old.cols.to.keep], tbl.out[, new.cols.to.keep], by="loc")
    
} # refimoForPiq
#------------------------------------------------------------------------------------------------------------------------
test.refimoForPiq <- function(tbl)
{
   printf("--- test.refimoForPiq")
   chrom <- "chr19"
       # start with a two-row piq table, duplicated loc, 2 motifs mapped
       # we want to see 3 different motifs here, each 1 base longer than the 2 oddly reported by piz
   start <- 45423532
   end <- 45423542
   tbl.p <- createPiqTable(chrom, start, end, collapseOnStrand=TRUE)

      # make sure the result is as expected - and NOT what we want to see
   checkEquals(dim(tbl.p), c(2, 12))   
   checkTrue(all(tbl.p$loc == "chr19:45423532-45423542"))
   checkEquals(sort(tbl.p$motif.p),  c("MA0039.2", "MA0599.1"))

   tbl.pf <- refimoForPiq(tbl.p)
   checkEquals(dim(tbl.pf), c(3, 14))
   checkEquals(colnames(tbl.pf), c("loc", "samplecount.p", "length.p", "score1.p.median", "score1.p.best",
                                   "score2.p.median", "score2.p.best", "score3.p.median", "score3.p.best",
                                   "score4.p.median", "score4.p.best", "motifname", "motifscore", "pval"))
   checkEquals(sort(tbl.pf$motifname), c("KLF16_DBD", "MA0741.1", "MA0746.1"))

      # repeat with more locs
   start <- 45423500
   end <- 45423550
   tbl.p <- createPiqTable(chrom, start, end, collapseOnStrand=TRUE)
   checkEquals(dim(tbl.p), c(10, 12))
   checkEquals(sort(unique(tbl.p$motif.p)), 
      c("MA0039.2", "MA0056.1", "MA0599.1", "MA0657.1", "MA0685.1", "MA0740.1", "MA0741.1", "MA0742.1", "MA0746.1", "MA0747.1"))

   tbl.pf <- refimoForPiq(tbl.p)
   checkEquals(dim(tbl.pf), c(6, 14))
   checkEquals(sort(unique(tbl.pf$motifname)), c("KLF16_DBD", "Klf12_DBD", "MA0741.1", "MA0742.1", "MA0746.1", "MA0747.1"))
   
} # test.refimoForPiq
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
ensemble <- function(chrom, start, end, test.motifs=NA, test.locs=NA)
{
   tbl.w <- createWellingtonTable(chrom, start, end, collapseOnStrand=TRUE)

   if(!any(is.na(test.motifs)))
       tbl.w <- subset(tbl.w, motif.w %in% test.motifs)
   if(!any(is.na(test.locs)))
       tbl.w <- subset(tbl.w, loc %in% test.locs)
   
   tbl.h <- createHintTable(chrom, start, end, collapseOnStrand=TRUE)
   if(!any(is.na(test.motifs)))
       tbl.h <- subset(tbl.h, motif.h %in% test.motifs)
   if(!any(is.na(test.locs)))
       tbl.h <- subset(tbl.h, loc %in% test.locs)

   tbl.p <- createPiqTable(chrom, start, end, collapseOnStrand=TRUE)

   if(!any(is.na(test.motifs)))
       tbl.p <- subset(tbl.p, motif.p %in% test.motifs)
   if(!any(is.na(test.locs)))
       tbl.p <- subset(tbl.p, loc %in% test.locs)


      tbl.merged <- merge(tbl.w, tbl.p, by.x=c("loc", "motif.w"), by.y=c("loc", "motif.p"), all=TRUE)
      tbl.merged <- merge(tbl.merged, tbl.h, by.x=c("loc", "motif.w"), by.y=c("loc", "motif.h"), all=TRUE)

   tbl.merged <- merge(tbl.w, tbl.h, by.x=c("loc", "motif.w"), by.y=c("loc", "motif.h"), all=TRUE)
   colname.pos <- grep("motif.w", colnames(tbl.merged))

   colnames(tbl.merged)[colname.pos] <- "motif"   # motif.h dropped in merge
   
   tbl.chipseq <- getHits(db.chipseq, chrom, start, end)

   printf ("found %d chipseq hits in %d bases", nrow(tbl.chipseq),  1 + end - start)
   stopifnot(nrow(tbl.chipseq) > 0)
   tbl.fimo <- getFimoHits(chrom, start, end)
   tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")
   tbl.csWithFimo <- addFimoRegions(tbl.chipseq, tbl.fimo)
   if(!any(is.na(test.locs)))
      tbl.csWithFimo <- subset(tbl.csWithFimo, loc %in% test.locs)

   tbl.csWithTrueFimo <- subset(tbl.csWithFimo, bindingSite==TRUE)
   
   tbl.csft <- tbl.csWithTrueFimo[, c("loc", "motifname", "name", "score1", "motifscore", "pval")]
   colnames(tbl.csft) <- c("loc", "csmotif", "csTF", "csscore", "motifscore", "motifpval")
   tbl.wfc <- merge(tbl.merged, tbl.csft, by.x=c("loc", "motif"), by.y=c("loc","csmotif"), all=TRUE)
   classInfo <- sapply(tbl.wfc, class)
   for(colname in names(classInfo)){
      if(classInfo[[colname]] == "factor"){
          printf("correctoring character class in column '%s'", colname)
          tbl.wfc[[colname]] <- as.character(tbl.wfc[[colname]])
          } # if factor
      } # for colname

   
   motifscore <- apply(tbl.wfc[, c("score2.w.best", "score2.h.best", "motifscore")], 1, function(row) max(row, na.rm=TRUE))
   motifpval <-  apply(tbl.wfc[, c("score3.w.best", "score3.h.best", "motifpval")], 1, function(row) min(row, na.rm=TRUE))

   keepers <- c("loc", "motif", "samplecount.w",   "score1.w.best", "samplecount.h", "score1.h.best", "csTF", "csscore")

   tbl.wfcTrimmed <- tbl.wfc[, keepers]
   colnames(tbl.wfcTrimmed) <- c("loc", "motif", "samplecount.w", "score.w", "samplecount.h", "score.h", "csTF", "csscore")
   tbl.wfcTrimmed$motifscore <- motifscore
   tbl.wfcTrimmed$motifpval <- motifpval
   invisible(unique(tbl.wfcTrimmed))

} # ensemble
#------------------------------------------------------------------------------------------------------------------------
ensemble.flex <- function(chrom, start, end, methods, test.motifs=NA, test.locs=NA)
{
   stopifnot(all(methods %in% c("wellington", "hint", "piq", "chipseq")))

   tbl.out <- data.frame()

   if("wellington" %in% methods){
      tbl.w <- createWellingtonTable(chrom, start, end, collapseOnStrand=TRUE)
      if(!any(is.na(test.motifs)))
         tbl.w <- subset(tbl.w, motif.w %in% test.motifs)
      if(!any(is.na(test.locs)))
         tbl.w <- subset(tbl.w, loc %in% test.locs)
      if(nrow(tbl.out) == 0)
         tbl.out <- tbl.w
      else
         tbl.out <- merge(tbl.out, tbl.w, by.x=c("loc", "motif"), by.y=c("loc", "motif.w"), all=TRUE)
      colnames(tbl.out)[grep("motif.w", colnames(tbl.out))] <- "motif"
      } # wellington
   
   if("hint" %in% methods){
      tbl.h <- createHintTable(chrom, start, end, collapseOnStrand=TRUE)
      if(!any(is.na(test.motifs)))
          tbl.h <- subset(tbl.h, motif.h %in% test.motifs)
      if(!any(is.na(test.locs)))
         tbl.h <- subset(tbl.h, loc %in% test.locs)
      if(nrow(tbl.out) == 0)
         tbl.out <- tbl.h
      else
         tbl.out <- merge(tbl.out, tbl.h, by.x=c("loc", "motif"), by.y=c("loc", "motif.h"), all=TRUE)
      } # wellington

   if("piq" %in% methods){
      tbl.p <- createPiqTable(chrom, start, end, collapseOnStrand=TRUE)
      if(!any(is.na(test.motifs)))
          tbl.p <- subset(tbl.p, motif.p %in% test.motifs)
      if(!any(is.na(test.locs)))
         tbl.p <- subset(tbl.p, loc %in% test.locs)
      tbl.p <- refimoForPiq(tbl.p)
         # piq does not return motifscores, nor pval
         # refimoForPiq adds that
         # but the code which follows, below, does not expect that from piq
         # since it was written before the (temporary) need for refimoForPiq 
         # was understood.  so remove those unexpected columns here
      cols.to.remove <- match(c("motifscore", "pval"), colnames(tbl.p))
      if(length(cols.to.remove) > 0)
         tbl.p <- tbl.p[, -(cols.to.remove)]
      if(nrow(tbl.out) == 0)
         tbl.out <- tbl.p
      else{
         tbl.out <- merge(tbl.out, tbl.p, by.x=c("loc", "motif"), by.y=c("loc", "motifname"), all=TRUE)
         }
      } # piq

    if("chipseq" %in% methods){
       tbl.chipseq <- getHits(db.chipseq, chrom, start, end)
       printf ("found %d chipseq hits in %d bases", nrow(tbl.chipseq),  1 + end - start)
       stopifnot(nrow(tbl.chipseq) > 0)
       tbl.fimo <- getFimoHits(chrom, start, end)
       tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")
       tbl.csWithFimo <- addFimoRegions(tbl.chipseq, tbl.fimo)
       if(!any(is.na(test.locs)))
          tbl.csWithFimo <- subset(tbl.csWithFimo, loc %in% test.locs)
       tbl.csWithTrueFimo <- subset(tbl.csWithFimo, bindingSite==TRUE)
       tbl.csft <- tbl.csWithTrueFimo[, c("loc", "motifname", "name", "score1", "motifscore", "pval")]
       colnames(tbl.csft) <- c("loc", "csmotif", "csTF", "csscore", "motifscore", "motifpval")
       if(nrow(tbl.out) == 0)
          tbl.out <- tbl.csft
       else
          tbl.out <- merge(tbl.out, tbl.csft, by.x=c("loc", "motif"), by.y=c("loc","csmotif"), all=TRUE)
       classInfo <- sapply(tbl.out, class)
       for(colname in names(classInfo)){
          if(classInfo[[colname]] == "factor"){
             printf("correctoring character class in column '%s'", colname)
             tbl.out[[colname]] <- as.character(tbl.out[[colname]])
             } # if factor
          } # for colname
      } # if chipseq
   
      # wellington, hint and chipseq all report motif scores, all obtained from fimo
      # collapse these the best from each of these (max score, min pval) in the unexpected
      # case they differ

   motif.score.columns <- c()
   motif.pval.columns <- c()
   other.columns.to.drop <- c()

   if("wellington" %in% methods){
       motif.score.columns <- c(motif.score.columns, "score2.w.best")
       motif.pval.columns  <- c(motif.pval.columns,  "score3.w.best")
       other.columns.to.drop <- c(other.columns.to.drop, "score2.w.median")
       other.columns.to.drop <- c(other.columns.to.drop, "score3.w.median")
       }

   if("hint" %in% methods){
       motif.score.columns <- c(motif.score.columns, "score2.h.best")
       motif.pval.columns  <- c(motif.pval.columns,  "score3.h.best")
       other.columns.to.drop <- c(other.columns.to.drop, "score2.h.median")
       other.columns.to.drop <- c(other.columns.to.drop, "score3.h.median")
       }

   if("chipseq" %in% methods){
       motif.score.columns <- c(motif.score.columns, "motifscore")
       motif.pval.columns  <- c(motif.pval.columns,  "motifpval")
       }

   if(length(motif.score.columns) > 0) {
       motifscore <- apply(tbl.out[, motif.score.columns, drop=FALSE], 1,
                           function(row) if (all(is.na(row))) return(0) else return(max(row, na.rm=TRUE)))
       columns.to.delete <- match(motif.score.columns, colnames(tbl.out))
       stopifnot(length(columns.to.delete) > 0)
       tbl.out <- tbl.out[, -columns.to.delete]
       tbl.out$motifscore <- motifscore
       }

   if(length(motif.pval.columns) > 0){
       motifpval <- apply(tbl.out[, motif.pval.columns, drop=FALSE], 1,
                           function(row) if (all(is.na(row))) return(1) else return(min(row, na.rm=TRUE)))
       columns.to.delete <- match(motif.pval.columns, colnames(tbl.out))
       stopifnot(length(columns.to.delete) > 0)
       tbl.out <- tbl.out[, -columns.to.delete]
       tbl.out$motifpval <- motifpval
       }

    if(length(other.columns.to.drop) > 0){
       other.columns.to.drop.indices <- match(other.columns.to.drop, colnames(tbl.out))
       stopifnot(length(other.columns.to.drop.indices) > 0)
       tbl.out <- tbl.out[, -other.columns.to.drop.indices]
       }                                       
  
 # keepers <- c("loc", "motif", "samplecount.w",   "score1.w.best", "samplecount.h", "score1.h.best", "csTF", "csscore")
   # 
   # tbl.wfcTrimmed <- tbl.wfc[, keepers]
   # colnames(tbl.wfcTrimmed) <- c("loc", "motif", "samplecount.w", "score.w", "samplecount.h", "score.h", "csTF", "csscore")
   # tbl.wfcTrimmed$motifscore <- motifscore
   # tbl.wfcTrimmed$motifpval <- motifpval
   # invisible(unique(tbl.wfcTrimmed))

   tbl.out
   
} # ensemble.flex
#------------------------------------------------------------------------------------------------------------------------
test.ensemble_hint_wellington_empty_chipseq <- function()
{
   printf("--- test.ensemble_hint_wellington_empty_chipseq")

   chrom <- "chr19"
   start <- 45423000
   end   <- 45424000
     # 229 bases
   test.motifs <- c("MA0512.2", "MA0677.1", "MA0855.1", "MA0856.1")
   test.locs <- c("chr19:45423537-45423550")

   ft <- ensemble(chrom, start, end, test.motifs=test.motifs, test.locs=test.locs)

      # for snooping around fimo scores, repeated for different methods
   coi <- c(grep("score", colnames(ft)), grep("motif", colnames(ft)), grep("csTF", colnames(ft)))
   checkEquals(nrow(ft), 4)
      # in preparation for collapsing all these fimo scores, make sure they are all equal
      # which they should be, since all are for the same motif at the same location

   # now test without the loc constraint, just one motif
   test.motifs <- "MA0677.1"
   test.locs <- c("chr19:45423532-45423542", "chr19:45423532-45423543",
                  "chr19:45423537-45423550", "chr19:45423560-45423567",
                  "chr19:45423560-45423576", "chr19:45423562-45423572",
                  "chr19:45423562-45423575", "chr19:45423562-45423576")

                  
     # with these cherry-picked locs, and no constraint on motifs, we see all combinations
     # of wellington, hint and chipseq
   test.locs <- c("chr19:45428888-45428899", "chr19:45428888-45428899", "chr19:45428888-45428899", "chr19:45428888-45428899",
                  "chr19:45428888-45428899", "chr19:45428905-45428919", "chr19:45428905-45428919", "chr19:45428935-45428949",
                  "chr19:45428937-45428944", "chr19:45428939-45428946", "chr19:45428940-45428952", "chr19:45428941-45428951",
                  "chr19:45428941-45428951", "chr19:45428948-45428960", "chr19:45428948-45428960", "chr19:45428948-45428960",
                  "chr19:45428951-45428961", "chr19:45428951-45428961", "chr19:45428951-45428961", "chr19:45428955-45428968",
                  "chr19:45428956-45428970", "chr19:45428957-45428967", "chr19:45428957-45428971", "chr19:45428958-45428968",
                  "chr19:45428958-45428968", "chr19:45428958-45428972", "chr19:45428959-45428969", "chr19:45428959-45428969",
                  "chr19:45428959-45428969", "chr19:45428959-45428969", "chr19:45428960-45428970", "chr19:45428960-45428970",
                  "chr19:45428960-45428973", "chr19:45428961-45428967", "chr19:45428971-45428984", "chr19:45428982-45428991",
                  "chr19:45428982-45428991", "chr19:45428991-45429000", "chr19:45428991-45429003", "chr19:45428995-45429008",
                  "chr19:45428996-45429005", "chr19:45429072-45429082") 
   chrom <- "chr19"
   start <- 45420000
   end   <- 45430000
   ft2 <- ensemble(chrom, start, end, test.motifs=NA, test.locs=NA)
   checkEquals(dim(ft2), c(540, 10))

} # test.ensemble_hint_wellington_empty_chipseq
#------------------------------------------------------------------------------------------------------------------------
test.ensemble_hint_wellington_chipseq <- function()
{
   printf("--- test.ensemble_hint_wellington_chipseq")

     # now  call ensemble with a region and locs in which a TF (CTCF) is found, and some footprints too

   chrom <- "chr19"
   start <- 45420000
   end   <- 45430000
   test.locs <- c("chr19:45423529-45423545", "chr19:45423530-45423547",
                  "chr19:45423592-45423608", "chr19:45423594-45423607",
                  "chr19:45423910-45423922", "chr19:45423911-45423920")
   test.motifs <- NA
   ft <- ensemble(chrom, start, end, test.motifs=test.motifs, test.locs=test.locs[1])
   checkEquals(nrow(ft), 1)
   checkEquals(colnames(ft),  c("loc", "motif", "samplecount.w", "score.w", "samplecount.h",
                                "score.h", "csTF", "csscore", "motifscore",  "motifpval"))

   checkEquals(ft$motif, "MA0685.1")
   checkEquals(ft$csTF, "SP1")
   checkEquals(ft$samplecount.w, 4)
   checkEquals(ft$samplecount.h, 3)

   chrom <- "chr19"
   start <- 45420000
   end   <- 45430000
   test.motifs <- c("MA0740.1", "MA0742.1")
   test.locs <- c("chr19:45423531-45423544", "chr19:45423531-45423545")
   ft <- ensemble(chrom, start, end, test.motifs=test.motifs, test.locs=test.locs)

} # test.ensemble_hint_wellington_chipseq
#------------------------------------------------------------------------------------------------------------------------
test.ensemble.flex.1.loc.only <- function()
{
   printf("--- test.ensemble.flex.1.loc.only")
    
   chrom <- "chr19"
   start <- 45420000
   end   <- 45430000
   test.locs <- c("chr19:45423532-45423542")
   test.motifs <- NA

     # cherry-piced locs which give results for piq, wellington and hint
     #, "chr19:45423532-45423543", "chr19:45423532-45423544",
     #             "chr19:45423537-45423550", "chr19:45423560-45423577", "chr19:45423562-45423576",
     #             "chr19:45423563-45423573")

   
   ft.w <- ensemble.flex(chrom, start, end, "wellington", test.motifs=test.motifs, test.locs=test.locs[1])
   checkEquals(dim(ft.w), c(3,8))

   ft.h <- ensemble.flex(chrom, start, end, "hint", test.motifs=test.motifs, test.locs=test.locs[1])
   checkEquals(dim(ft.h), c(3,8))

   ft.p <- ensemble.flex(chrom, start, end, "piq", test.motifs=test.motifs, test.locs=test.locs[1])
   checkEquals(dim(ft.p), c(3, 12))

   ft.c <- ensemble.flex(chrom, start, end, "chipseq", test.motifs=test.motifs, test.locs=test.locs[1])
   checkEquals(dim(ft.c), c(3, 6))
   
   methods <- c("wellington", "hint", "piq", "chipseq")
   ft.all <- ensemble.flex(chrom, start, end, methods, test.motifs=test.motifs, test.locs=test.locs[1])
   checkEquals(dim(ft.all), c(3, 24))

} # test.ensemble.flex.1.loc.only
#------------------------------------------------------------------------------------------------------------------------
test.ensemble.flex.all.locs <- function()
{
   printf("--- test.ensemble.flex.all.locs")
    
   chrom <- "chr19"
   start <- 45423800
   end   <- 45424000

   test.locs <- NA
   test.motifs <- NA

   ft.w <- ensemble.flex(chrom, start, end, "wellington", test.motifs=test.motifs, test.locs=test.locs)
   checkEquals(dim(ft.w), c(35,8))

   ft.h <- ensemble.flex(chrom, start, end, "hint", test.motifs=test.motifs, test.locs=test.locs)
   checkEquals(dim(ft.h), c(42,8))

   ft.p <- ensemble.flex(chrom, start, end, "piq", test.motifs=test.motifs, test.locs=test.locs)
   checkEquals(dim(ft.p), c(2, 12))

   ft.c <- ensemble.flex(chrom, start, end, "chipseq", test.motifs=test.motifs, test.locs=test.locs)
   checkEquals(dim(ft.c), c(47, 6))
   
   methods <- c("wellington", "hint", "piq", "chipseq")
   ft.all <- ensemble.flex(chrom, start, end, methods, test.motifs=test.motifs, test.locs=test.locs)
   checkEquals(dim(ft.all), c(69, 24))

} # test.ensemble.flex.all.locs
#------------------------------------------------------------------------------------------------------------------------
test.ensemble <- function()
{
   printf("--- test.ensemble")

     # now  call ensemble with a region and locs in which a TF (CTCF) is found, and some footprints too

   chrom <- "chr19"
   start <- 45423662
   end   <- 45423850
   test.locs <- "chr19:45423814-45423832"
   test.motifs <- NA
   ft <- ensemble(chrom, start, end, test.motifs=test.motifs, test.locs=test.locs)

   checkEquals(dim(ft), c(33, 22))
   tfHits <- ft[which(!is.na(ft$csTF)), "csTF"]
   checkEquals(length(tfHits), 1)
   checkEquals(tfHits, "CTCF")
    
     # test in one 12 base region, very rich in motifs, includes all most combinations of present & missing motifs
   chrom <- "chr19"
   start <- 44903672
   end   <- 44903900

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
cleanFeatureTable <- function(tbl, totalSampleCount)
{
   x <- tbl
   keepers <- c("loc", "motif", "samplecount.w", "length.w", "score1.w.best", "samplecount.h", "length.h", "score1.h.best",
                "samplecount.p", "length.p", "score1.p.best", "score2.p.best", "score3.p.best", "score4.p.best",
                "csTF", "csscore", "motifscore", "motifpval")

   missing <- setdiff(keepers, colnames(tbl))
   if(length(missing) > 0){
       printf("missing columns in cleanFeatureTable: %s", paste(missing, collapse=","))
       stop()
      }
   x <- x[, keepers]
   x$samplecount.w[is.na(x$samplecount.w)] <- 0
   x$score1.w.best[is.na(x$score1.w.best)] <- -99

   x$length.w[is.na(x$length.w)] <- 0
   x$length.h[is.na(x$length.h)] <- 0
   x$length.p[is.na(x$length.p)] <- 0
   
   x$samplecount.h[is.na(x$samplecount.h)] <- 0
   x$score1.h.best[is.na(x$score1.h.best)] <- 0

   x$samplecount.p[is.na(x$samplecount.p)] <- 0
   x$score1.p.best[is.na(x$score1.p.best)] <- 0
   x$score2.p.best[is.na(x$score2.p.best)] <- 0
   x$score3.p.best[is.na(x$score3.p.best)] <- 0
   x$score4.p.best[is.na(x$score4.p.best)] <- 0

   x$csscore[is.na(x$csscore)] <- 0
   x$totalsamplecount <- totalSampleCount

   colnames(x) <- c("loc", "motif", "samplecount.w", "length.w", "score.w",
                    "samplecount.h", "length.h", "score.h",
                    "samplecount.p", "length.p", "score1.p", "score2.p", "score3.p", "score4.p",
                    "csTF", "csscore", "motifscore", "motifpval", "totalsamplecount")

   invisible(unique(x))
   
} # cleanFeatureTable
#------------------------------------------------------------------------------------------------------------------------
test.cleanFeatureTable <- function()
{
   printf("--- test.cleanFeatureTable")

   load ("ft.all.Rdata")
   checkEquals(dim(ft.all), c(69, 24))

   ft.clean <- cleanFeatureTable(ft.all, 19)
   checkEquals(dim(ft.clean), c(69, 19))

   checkTrue("motif" %in% colnames(ft.clean))
   checkTrue(!"motif.h" %in% colnames(ft.clean))
   checkTrue(!"motif.w" %in% colnames(ft.clean))
   checkTrue(!"motif.w" %in% colnames(ft.clean))
   checkEquals(nrow(ft.all), nrow(ft.clean))

   checkTrue(all(ft.clean$length.w >= 0))
   checkTrue(all(ft.clean$length.h >= 0))
   checkTrue(all(ft.clean$length.p >= 0))

   checkEquals(length(ft.clean$score.w), nrow(ft.clean))
   checkTrue(all(ft.clean$score.w <= 0))

   checkEquals(length(ft.clean$score.h), nrow(ft.clean))
   checkTrue(all(ft.clean$score.h >= 0))

   checkTrue(all(ft.clean$motifscore >= 0))
   checkTrue(all(ft.clean$motifpval >= 0))
   checkTrue(all(ft.clean$motifpval <= 1))
   checkTrue(all(ft.clean$csscore >= 0))

} # test.cleanFeatureTable
#------------------------------------------------------------------------------------------------------------------------
test.big10M <- function()
{
   ft.big <- ensemble("chr19", 0, 1000000)
   save(ft.big, file="ft.big.RData")
   
} # test.big10M
#------------------------------------------------------------------------------------------------------------------------
test.run50k <- function()
{
   printf("--- test.run50k")

   chrom <- "chr19"
   start <- 42000000
   end   <- 42050000

   print(system.time(ft <- ensemble(chrom, start, end)))
   totalSampleCount <- nrow(dbGetQuery(db.hint, "select distinct sample_id from hits"))
   ft <- cleanFeatureTable(ft, totalSampleCount)

   tbl.hit.all <- subset(ft, !is.na(csTF) & samplecount.w > 0  & samplecount.h > 0)
   tbl.cs.only <- subset(ft, !is.na(csTF) & samplecount.w == 0 & samplecount.h == 0)
   tbl.w.only  <- subset(ft, is.na(csTF)  & samplecount.w > 0  & samplecount.h == 0)
   tbl.h.only  <- subset(ft, is.na(csTF)  & samplecount.w == 0 & samplecount.h > 0)
   tbl.wh.only <- subset(ft, is.na(csTF)  & samplecount.w > 0  & samplecount.h > 0)
   tbl.wc.only <- subset(ft, !is.na(csTF) & samplecount.w > 0  & samplecount.h == 0)
   tbl.hc.only <- subset(ft, !is.na(csTF) & samplecount.w == 0  & samplecount.h > 0)

   checkEquals(nrow(ft), sum(nrow(tbl.hit.all), nrow(tbl.cs.only), nrow(tbl.w.only),
                             nrow(tbl.h.only), nrow(tbl.wh.only), nrow(tbl.wc.only), nrow(tbl.hc.only)))
   
      # do some positive checks
   for(r in 1:nrow(tbl.hit.all)){
      loc <- locStringToBedTable(tbl.hit.all$loc[r])
      hits.raw <- with(loc, getHits(db.hint, chrom, start, end))
      hits <- subset(hits.raw, name==tbl.hit.all$motif[r])
      checkEquals(nrow(hits), tbl.hit.all$samplecount.h[r])
      hits.raw <- with(loc, getHits(db.wellington, chrom, start, end))
      hits <- subset(hits.raw, name==tbl.hit.all$motif[r])
      checkEquals(nrow(hits), tbl.hit.all$samplecount.w[r])
      } # for r

      # now some negative checks
   for(r in 1:nrow(tbl.w.only)){
      loc <- locStringToBedTable(tbl.w.only$loc[r])
      hits.raw <- with(loc, getHits(db.hint, chrom, start, end))
      if(nrow(hits.raw) > 0) {
         hits <- subset(hits.raw, name==tbl.w.only$motif[r])
         checkEquals(nrow(hits), 0)
         }
      } # for r
    
   for(r in 1:nrow(tbl.h.only)){
      loc <- locStringToBedTable(tbl.h.only$loc[r])
      hits.raw <- with(loc, getHits(db.wellington, chrom, start, end))
      if(nrow(hits.raw) > 0) {
         hits <- subset(hits.raw, name==tbl.h.only$motif[r])
         checkEquals(nrow(hits), 0)
         }
      } # for r
    
} # test.run50k
#------------------------------------------------------------------------------------------------------------------------
run <- function()
{
   chrom <- "chr19"
   start <- 42000000
   end   <- 47000000

   test.locs <- NA
   test.motifs <- NA
   methods <- c("wellington", "hint", "piq", "chipseq")
   ft  <- ensemble.flex(chrom, start, end, methods, test.motifs=test.motifs, test.locs=test.locs)
   tbl <- cleanFeatureTable(ft, totalSampleCount=18)  # 18 lymphoblast samples from encode

   filename <- sprintf("ft.%s.RData", format (Sys.time(), "%a.%b.%d.%Y-%H:%M:%S"))
   save(ft, file=filename)
   ft.clean <- cleanFeatureTable(ft, 18)
   save(ft.clean, file=sprintf("ftClean.%s.RData", format (Sys.time(), "%a.%b.%d.%Y-%H:%M:%S")))

   printf("tbl.hit.all: %d", nrow(subset(ft.clean, !is.na(csTF) & samplecount.w > 0  & samplecount.h > 0)))
   printf("tbl.cs.only: %d", nrow(subset(ft.clean, !is.na(csTF) & samplecount.w == 0 & samplecount.h == 0)))
   printf("tbl.w.only: %d", nrow(subset(ft.clean, is.na(csTF)  & samplecount.w > 0  & samplecount.h == 0)))
   printf("tbl.h.only: %d", nrow(subset(ft.clean, is.na(csTF)  & samplecount.w == 0 & samplecount.h > 0)))
   printf("tbl.wh.only: %d", nrow(subset(ft.clean, is.na(csTF)  & samplecount.w > 0  & samplecount.h > 0)))
   printf("tbl.wc.only: %d", nrow(subset(ft.clean, !is.na(csTF) & samplecount.w > 0  & samplecount.h == 0)))
   printf("tbl.hc.only: %d", nrow(subset(ft.clean, !is.na(csTF) & samplecount.w == 0  & samplecount.h > 0)))

   
} # run
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
    run()
