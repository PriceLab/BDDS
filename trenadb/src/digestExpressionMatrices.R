library(org.Hs.eg.db)
f <- "/Volumes/local/Cory/Alzheimers/repo_counts/rosmap_counts_matrix_normalized.txt"
tbl <- read.table(f, sep="\t", nrows=-1, header=TRUE, as.is=TRUE)    # 642 columns
colnames(tbl) <- sub("^X", "s", colnames(tbl))
ensembl.geneIDs <- sub("\\..*$", "", tbl$ensembl_id)
tbl.map <- select(org.Hs.eg.db, keys=ensembl.geneIDs, keytype="ENSEMBL", columns=c("ENSEMBL", "SYMBOL"))
dups <- which(duplicated(tbl.map$ENSEMBL))  # 799/66658
if(length(dups) > 0)
   tbl.map <- tbl.map[-dups,]
tbl <- cbind(tbl.map, tbl)
dim(tbl) # [1] 65859   641
unmapped <- which(is.na(tbl$SYMBOL))  #  38199
length(unmapped)    # 31848
if(length(unmapped) > 0)
  tbl <- tbl[-unmapped,]  
dups <- which(duplicated(tbl$SYMBOL))
if(length(dups) > 0)
    tbl <- tbl[-dups,]

rownames(tbl) <- tbl$SYMBOL
delete.these.columns <- match(c("ENSEMBL", "SYMBOL", "ensembl_id"), colnames(tbl))
if(length(delete.these.columns) > 0)
    tbl <- tbl[, -delete.these.columns]

dim(tbl) # 24593   638
mtx <- as.matrix(tbl)
mtx[1:5, 1:5]
#         s01_120405 s02_120405 s03_120405 s04_120405 s05_120405
# TSPAN6     4.321084   2.550606   2.626589   3.685091   3.442105
# TNMD       0.076254   0.067835   0.056486   0.057805   0.034769
# DPM1       3.380612   6.227278   4.462377   6.734323   4.207018
# SCYL3      2.897668   4.884139   5.224935   4.899003   3.720255
# C1orf112   0.432108   0.624084   0.536615   0.881532   0.451994
range(mtx) # [1]      0.0 239943.8 
save(mtx, file="rosmap_rnaseq_fpkm_geneSymbols_24593x638.RData")

