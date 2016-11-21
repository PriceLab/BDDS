library(RPostgreSQL)
library(TReNA)
library(RCyjs)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
genome.db.uri    <- "postgres://whovian/hg38"             # has gtf and motifsgenes tables
footprint.db.uri <- "postgres://whovian/wholeBrain"       # has hits and regions tables
if(!exists("mtx.rosmap")){
   #load("~/s/work/priceLab/cory/module-109/rosmap_rnaseq_fpkm_geneSymbols_24593x638.RData")
     # copy whovian file to your laptop, needed to render network into your locally running web browser
     # here is a tempoary load command for testing on whovian.  todo: this is brittle - fix!
   load("/users/pshannon/tmp/rosmap_rnaseq_fpkm_geneSymbols_24593x638.RData")
   mtx.rosmap <- mtx  # 24593   638
   }
#------------------------------------------------------------------------------------------------------------------------
getTSSTable <- function()
{
   db.gtf <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="gtf", host="whovian")
   query <- "select * from hg38human where moleculetype='gene' and gene_biotype='protein_coding'"
   tbl <- dbGetQuery(db.gtf, query) [, c("chr", "gene_name", "start", "endpos", "strand")]

} # getTSSTable
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.tss"))
    tbl.tss <- getTSSTable()
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test.createModel()
    
} # runTests
#------------------------------------------------------------------------------------------------------------------------
createModel <- function(target.gene, promoter.shoulder, 
                        mtx.expression=NA, mtx.classification=NA,
                        absolute.lasso.beta.min=0.01,
                        randomForest.purity.min=5,
                        absolute.expression.correlation.min=0.1)
{
   stopifnot(target.gene %in% rownames(mtx.expression))
    
   #query <- sprintf("select * from hg38human where gene_name='%s' and moleculetype='gene'", target.gene)
   #tbl.tmp <- dbGetQuery(db.gtf, query)
   #gene.info <- list(chrom=tbl.tmp[1, "chr"], start=tbl.tmp[1, "start"])

   fpf <- FootprintFinder(genome.db.uri, footprint.db.uri, quiet=TRUE)
   tbl.fp <- getFootprintsForGene(fpf, target.gene, size.upstream=promoter.shoulder,
                                  size.downstream=promoter.shoulder)
   candidate.tfs <- sort(unique(tbl.fp$tf_name))
   candidate.tfs <- intersect(rownames(mtx.expression), candidate.tfs)
   goi <- sort(unique(c(target.gene, candidate.tfs)))

   mtx.matched <- mtx.expression[goi,]
   mtx.matched <- asinh(mtx.matched)    

   trena.lasso <- TReNA(mtx.matched, solver="lasso")
   trena.ranfor <- TReNA(mtx.matched, solver="randomForest")
    
   out.lasso <- solve(trena.lasso, target.gene, candidate.tfs)
   out.ranfor <- solve(trena.ranfor, target.gene, candidate.tfs)

   tbl.01 <- out.lasso
   tbl.01 <- tbl.01[, -(grep("^intercept$", colnames(tbl.01)))]
   tbl.01$gene <- rownames(out.lasso)
   rownames(tbl.01) <- NULL
   tbl.01 <- subset(tbl.01, abs(beta) >= absolute.lasso.beta.min &
                            abs(gene.cor) >= absolute.expression.correlation.min)

   tbl.02 <- out.ranfor$edges
   tbl.02$gene <- rownames(tbl.02)
   rownames(tbl.02) <- NULL
   tbl.02 <- subset(tbl.02, IncNodePurity >= randomForest.purity.min  &
                            abs(gene.cor) >= absolute.expression.correlation.min)

   tbl.03 <- merge(tbl.01, tbl.02)
   rownames(tbl.03) <- tbl.03$gene
   tbl.03 <- tbl.03[, -(grep("^gene$", colnames(tbl.03)))]    
   tbl.04 <- tbl.03[order(abs(tbl.03$gene.cor), decreasing=TRUE),]
   
   footprint.start <- unlist(lapply(rownames(tbl.04), function(gene) subset(tbl.fp, tf_name==gene)$mfpstart[1]))
   gene.info <- subset(tbl.tss, gene_name==target.gene)[1,]
   if(gene.info$strand  == "+")
      gene.start <- gene.info$start
   else
      gene.start <- gene.info$end
   distance <- gene.start - footprint.start
   tbl.04$distance <- distance
   tbl.04$beta[which(is.na(tbl.04$beta))] <- 0
   tbl.04$IncNodePurity[which(is.na(tbl.04$IncNodePurity))] <- 0

   tbl.04

} # createModel
#------------------------------------------------------------------------------------------------------------------------
test.createModel <- function()
{
   printf("--- test.createModel")
   tbl <- createModel("TREM2", promoter.shoulder=100,
                      mtx.expression=mtx.rosmap,
                      absolute.lasso.beta.min=0.01,
                      randomForest.purity.min=5,
                      absolute.expression.correlation.min=0.1)

   checkEquals(ncol(tbl), 4)
   checkEquals(colnames(tbl), c("gene.cor", "beta", "IncNodePurity", "distance"))
   checkTrue(nrow(tbl) >= 4)
   checkTrue(all(c("ELF4", "FLI1", "CEBPA", "ELK3") %in% head(rownames(tbl))))

     # eliminate thresholds, ensure that more tfs are returned
   tbl.2 <- createModel("TREM2", promoter.shoulder=100,
                      mtx.expression=mtx.rosmap,
                      absolute.lasso.beta.min=0.0,
                      randomForest.purity.min=0,
                      absolute.expression.correlation.min=0.)
   checkEquals(ncol(tbl.2), 4)
   checkTrue(nrow(tbl.2) > 25)
    
} # test.createModel
#------------------------------------------------------------------------------------------------------------------------
renderAsNetwork <- function(tbl, target.gene)
{
   tfs <- rownames(tbl)
   footprints <- unlist(lapply(tbl$distance, function(x)
       if(x < 0)
         sprintf("fp.downstream.%05d", abs(x))
       else
           sprintf("fp.upstream.%05d", x)))
   
   tbl$footprint <- footprints

   all.nodes <- unique(c(target.gene, tfs, footprints))
   
   g <- graphNEL(nodes=all.nodes, edgemode = "directed")
   nodeDataDefaults(g, attr = "type") <- "undefined"
   nodeDataDefaults(g, attr = "label") <- "default node label"
   nodeDataDefaults(g, attr = "distance") <- 0
   nodeDataDefaults(g, attr = "gene.cor") <- 0
   nodeDataDefaults(g, attr = "beta") <- 0
   nodeDataDefaults(g, attr = "purity") <- 0

   edgeDataDefaults(g, attr = "edgeType") <- "undefined"
   edgeDataDefaults(g, attr = "beta") <- 0
   edgeDataDefaults(g, attr = "purity") <- 0

   nodeData(g, target.gene, "type") <- "targetGene"
   nodeData(g, tfs, "type")         <- "TF"
   nodeData(g, footprints, "type")  <- "footprint"
   nodeData(g, all.nodes, "label")  <- all.nodes
   nodeData(g, footprints, "label") <- tbl$distance
   nodeData(g, footprints, "distance") <- tbl$distance

   nodeData(g, tfs, "gene.cor") <- tbl$gene.cor
   nodeData(g, tfs, "beta") <- tbl$gene.cor
   nodeData(g, tfs, "purity") <- tbl$IncNodePurity

   g <- graph::addEdge(rownames(tbl), tbl$footprint, g)
   edgeData(g, rownames(tbl), tbl$footprint, "edgeType") <- "bindsTo"
   
   g <- graph::addEdge(tbl$footprint, target.gene, g)
   edgeData(g, tbl$footprint, target.gene, "edgeType") <- "regulatorySiteFor"
   
   #edgeData(g, tfs, target.gene, "beta") <- tbl$beta
   #edgeData(g, tfs, target.gene, "purity") <- tbl$IncNodePurity

   rcy <- RCyjs(10000:10100, title="TReNA", graph=g)
   httpSetStyle(rcy, "style.js")
   rcy

} # renderAsNetwork
#------------------------------------------------------------------------------------------------------------------------
test.renderAsNetwork <- function()
{
   printf("--- test.renderAsNetwork")
  
   if(!exists("tbl.trem2"))
      load("tbl.trem2.RData")

   tbl <- subset(tbl.trem2, !is.na(IncNodePurity))  # just 6 TF targets
   target.gene <- "TREM2"
   rcy <- renderAsNetwork(tbl, target.gene)
   layoutByFootprintPosition(rcy)

   rcy

} # test.renderAsNetwork
#------------------------------------------------------------------------------------------------------------------------
layoutByFootprintPosition <- function(rcy)
{
   target.gene.y <- -200
   tfs.y <- 200

   g <- rcy@graph

   fp.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "footprint")]
   xPos <- nodeData(g, fp.nodes, attr="distance")
   setPosition(rcy, data.frame(id=names(xPos), x=as.integer(xPos), y=0, stringsAsFactors=FALSE))
   target.gene <- names(which(nodeData(g, attr="type") == "targetGene"))
   setPosition(rcy, data.frame(id=target.gene, x=0, y=-200, stringsAsFactors=FALSE))

   fpTFedges <- inEdges(g)[fp.nodes]

   for(fp in names(fpTFedges)){
     tfs <- fpTFedges[[fp]]
     xpos.base <- as.integer(nodeData(g, fp, "distance"))
     for(tfi in 1:length(tfs)){
       tf <- tfs[tfi]
       ypos <- tfs.y + ((tfi-1) * 200)
       xpos <- xpos.base + ((tfi-1) * 100)
       setPosition(rcy, data.frame(id=tf, x=xpos, y=ypos, stringsAsFactors=FALSE))
       }
     } # for edge

   fit(rcy, 100)

} # layoutByFootprintPosition
#------------------------------------------------------------------------------------------------------------------------
demo <- function()
{
   target.gene <- "TREM2"
   target.gene <- "APOE"
   target.gene <- "MBP"
   target.gene <- "MOBP"

   tbl <- createModel(target.gene, promoter.shoulder=2000,
                      mtx.expression=mtx.rosmap,
                      absolute.lasso.beta.min=0.1,
                      randomForest.purity.min=3,
                      absolute.expression.correlation.min=0.05)

   rcy <- renderAsNetwork(tbl, target.gene)
   layoutByFootprintPosition(rcy)

} # demo
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()

