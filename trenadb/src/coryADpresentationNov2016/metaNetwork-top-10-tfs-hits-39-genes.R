if(!exists("renderAsNetwork"))  # a function in renderGeneModel.R
    source("~/github/BDDS/trenadb/src/renderGeneModel.R")

table.names <- load("/Volumes/local/Cory/for_Paul/all_targets.RData")

gene.names <- unlist(lapply(strsplit(table.names, split=".", fixed=TRUE), "[", 2))
stopifnot(length(table.names) == length(gene.names))

tbl.all <- data.frame(stringsAsFactors=FALSE)

threshold.IncNodePurity <- 1.0
threshold.beta <- 0.01
threshold.geneCorrelation <- 0.1

tfs.by.targetGene <- list()

for(i in 1:length(gene.names)){  # gene.names are the target genes proposed by all of cory's collabotors
    targetGene <- gene.names[i]
    tbl.gene <- eval(parse(text=table.names[i]))
    tbl.sub <- subset(tbl.gene, abs(beta) < 0.5)
    tbl.sub <- subset(tbl.sub, IncNodePurity >= threshold.IncNodePurity &
                                abs(beta) >= threshold.beta)
    if(nrow(tbl.sub) == 0){ # improvise
        printf ("no IncNodePurity above threshold for %s, max(abs(beta)): %f, max(abs(gene.cor)): %f, max(purity): %f",
                table.names[i], max(abs(tbl.gene$beta)), max(abs(tbl.gene$gene.cor)), max(tbl.gene$IncNodePurity))
       tbl.sub <- subset(tbl.gene, abs(gene.cor) > 0.2 & abs(beta) > 0.01)
       }
    if(nrow(tbl.sub) == 0){ # still nothing, skip to next one
       printf("no good tf candidates in %s", table.names[i])
       next;
       }
    tbl.sub$targetGene <- targetGene
    tfs.by.targetGene[[targetGene]] <- unique(tbl.sub$gene)
    printf("targetGene %10s, %d rows, %d tfs", targetGene, nrow(tbl.sub), length(unique(tbl.sub$gene)))
    tbl.all <- rbind(tbl.all, tbl.sub)
    } # for i

#  find the most common tfs, by footprint incidence
tail(sort(table(tbl.all$gene)), n=10)
# NPAS1  ELK3   SP9 STAT4  RXRA   SP1  EGR4 KLF16 KLF15  EGR3 
#   55    57    58    76    97   122   124   142   146   171 

# find the most common tfs, by target.gene count
tail(sort(table(unlist(tfs.by.targetGene)), use.names=FALSE), n=10)
#  FOXO4  PKNOX2    TBR1    ELF4   RUNX1    TAL1   ESRRG    RXRA NEUROD6   STAT4 
#      6       6       6       7       7       7       8      10      13      15 


top.10.tfs <- names(tail(sort(table(unlist(tfs.by.targetGene)), use.names=FALSE), n=10))
top.10.tbls <- list()

for(target.gene in gene.names){
   tmp <- subset(tbl.all, gene %in% top.10.tfs & targetGene == target.gene)
   printf("target gene %s, rows: %d", target.gene, nrow(tmp))
   if(nrow(tmp) > 0)
      top.10.tbls[[target.gene]] <- tmp
   }

g <- tableToReducedGraph(top.10.tbls)
rcy <- RCyjs(10000:10100, title=paste(targets, collapse=" + "), graph=g)
httpSetStyle(rcy, "style-reducedGraph-footprintCountsIgnored.js")

# restoreLayout(rcy, "stat4-first5targets.RData")
layout(rcy, "cose")
fit(rcy)

