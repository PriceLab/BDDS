#-------------------------------------------------------------------------------
readWellingtonTable <- function(directory, sampleID, nrows=NA, chromosome=NA)
{
  filename <- grep(sampleID, list.files(directory), v=TRUE)
  full.path <- file.path(directory, filename)
  
  if(!file.exists(full.path))
    return(data.frame)
  
  tbl <- read.table(full.path, sep="\t", as.is=TRUE)
  colnames(tbl) <- c("chrom", "start", "end", "name", "score", "strand")
  #tbl$chrom <- paste("chr", tbl$chrom, sep="")
  if(!is.na(chromosome))
    tbl <- subset(tbl, chrom==chromosome)
  
  if(!is.na(nrows))
    tbl <- tbl[1:nrows,]
  
  invisible(tbl)
  
} # readWellingtonTable
#-------------------------------------------------------------------------------
mergeFimoWithFootprints <- function(tbl.fp, sampleID)
{
  chromosome <- unique(tbl.fp$chrom)
  # enforce treatment of just one chromosome at a time
  stopifnot(length(chromosome) == 1)
  min.pos <- min(tbl.fp$start)
  max.pos <- max(tbl.fp$end)
  
  fimo.chromosome <- sub("chr", "", chromosome)
  query <- sprintf("select * from hg38 where chr='%s' and start >= %d and endpos <= %d",
                   fimo.chromosome, min.pos, max.pos)
  
  tbl.fimo <- dbGetQuery(db.fimo, query)
  colnames(tbl.fimo) <- c("motif", "chrom", "motif.start", "motif.end", "motif.strand", "fimo.score",
                          "fimo.pvalue", "empty", "motif.sequence")
  tbl.fimo <- tbl.fimo[, -grep("empty", colnames(tbl.fimo))]
  tbl.fimo$chrom <- paste("chr", tbl.fimo$chrom, sep="")
  
  gr.fimo <- with(tbl.fimo, GRanges(seqnames=chrom, IRanges(start=motif.start, end=motif.end)))
  
  # --- get some footprints
  
  gr.wellington <- with(tbl.fp,   GRanges(seqnames=chrom, IRanges(start=start, end=end)))
  tbl.overlaps <- as.data.frame(findOverlaps(gr.fimo, gr.wellington, type="within"))
  
  tbl.fimo$loc <- with(tbl.fimo, sprintf("%s:%d-%d", chrom, motif.start, motif.end))
  tbl.fimo$method <- "WELLINGTON"
  tbl.fimo$sample_id <- sampleID
  tbl.regions <- tbl.fimo[tbl.overlaps$queryHits,]
  
  tbl.regions <- cbind(tbl.regions, wellington.score=tbl.fp[tbl.overlaps$subjectHits, "score"])
  invisible(tbl.regions)
  
} # mergeFimoWithFootprints
#-------------------------------------------------------------------------------
splitTableIntoRegionsAndWellingtonHits <- function(tbl, minid)
{
  tbl.regions <- unique(tbl[, c("loc", "chrom", "motif.start", "motif.end")])
  colnames(tbl.regions) <- region.schema() # 29
  # c("loc", "chrom", "motif_start", "motif_end")
  
  new.locs <- setdiff(tbl.regions$loc, names(knownLocs))
  # enter these new.locs into the hash
  lapply(new.locs, function(loc) knownLocs[[loc]] <- 0)   
  printf("novel locs: %d new/%d possible (%d now known, %d dups)",
         length(new.locs), nrow(tbl.regions), length(names(knownLocs)), nrow(tbl.regions) - length(new.locs))
  
  if(length(new.locs) > 0){
    keepers <- match(new.locs, tbl.regions$loc)
    tbl.regions <- tbl.regions[keepers,]
  } else {
    tbl.regions <- data.frame()
  }
  
  tbl.hits <- tbl[, c("loc", "motif", "motif.strand", "sample_id", "method", "wellington.score",
                      "fimo.score", "fimo.pvalue")]
  tbl.hits$length <- with(tbl, 1 + motif.end - motif.start)
  tbl.hits$provenance <- minid
  tbl.hits$score4 <- NA
  tbl.hits$score5 <- NA
  tbl.hits$score6 <- NA
  tbl.hits$type <- "motif.in.footprint"
  coi <- c("loc", "type", "motif", "length", "motif.strand", "sample_id", "method", "provenance", 
           "wellington.score", "fimo.score", "fimo.pvalue", "score4", "score5", "score6")
  tbl.hits <- tbl.hits[, coi]
  colnames(tbl.hits) <- hit.schema()
  invisible(list(regions=tbl.regions, hits=tbl.hits))
  
} # splitTableIntoRegionsAndWellingtonHits    
#-------------------------------------------------------------------------------
