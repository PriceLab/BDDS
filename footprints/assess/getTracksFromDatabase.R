library(RPostgreSQL)
library(igvR)   # on pricelab repo at github

#--------------------------------------------------------------------------------
# preliminary steps, external to this script:
# 1) start igv
# 2) igv: load from file: chr19 lymphoblast bam file, one sample only:
#    ~/github/BDDS/footprints/functionalTests/data/ENCSR000DBY.19.chr.bam
# 3) igv: load wellington footprints
#    ~/github/BDDS/footprints/functionalTests/output/wellington/ENCSR000DBY.19.chr.bam.chr19.peaks.400.bed.WellingtonFootprints.FDR.0.01.bed 
# 4) igv zoom into chr19:44859094-44956080
#--------------------------------------------------------------------------------

githubRepoRoot <- "~/github"
full.path <- file.path(githubRepoRoot, "BDDS/trenadb/src/utils.R")
stopifnot(file.exists(full.path))
source(full.path)  # opens database handles, provides 'getHits'

chrom <- "chr19"
start <- 44859094
end   <- 44956080

igv <- igvR()
stopifnot(connected(igv))
                                                     #  rows   unique locs
tbl.w <- getHits(db.wellington, chrom, start, end)   #  2657           657
tbl.h <- getHits(db.hint, chrom, start, end)         #  4622          1599
tbl.p <- getHits(db.piq, chrom, start, end)          # 75357          3207

displayBedTable(igv, unique(tbl.w[, c("chrom", "start", "endpos")]), "wellington")
displayBedTable(igv, unique(tbl.h[, c("chrom", "start", "endpos")]), "hint")
displayBedTable(igv, unique(tbl.p[, c("chrom", "start", "endpos")]), "piq")

