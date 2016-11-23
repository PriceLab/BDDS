library(RPostgreSQL)
library(TReNA)
library(RCyjs)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
genome.db.uri    <- "postgres://whovian/hg38"             # has gtf and motifsgenes tables
footprint.db.uri <- "postgres://whovian/wholeBrain"       # has hits and regions tables
if(!exists("fpf"))
   fpf <- FootprintFinder(genome.db.uri, footprint.db.uri, quiet=TRUE)

#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx.rosmap")){
   load("~/s/work/priceLab/cory/module-109/rosmap_rnaseq_fpkm_geneSymbols_24593x638.RData")
     # copy whovian file to your laptop, needed to render network into your locally running web browser
     # here is a tempoary load command for testing on whovian.  todo: this is brittle - fix!
   #load("/users/pshannon/tmp/rosmap_rnaseq_fpkm_geneSymbols_24593x638.RData")
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
   test.tableToFullGraph()
   test.tableToReducedGraph()
    
} # runTests
#------------------------------------------------------------------------------------------------------------------------
createModel <- function(target.gene, promoter.shoulder, 
                        mtx.expression=NA, mtx.classification=NA,
                        absolute.lasso.beta.min=0.0,
                        randomForest.purity.min=1,
                        absolute.expression.correlation.min=0.1)
{
   stopifnot(target.gene %in% rownames(mtx.expression))
    
   #query <- sprintf("select * from hg38human where gene_name='%s' and moleculetype='gene'", target.gene)
   #tbl.tmp <- dbGetQuery(db.gtf, query)
   #gene.info <- list(chrom=tbl.tmp[1, "chr"], start=tbl.tmp[1, "start"])

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

   tbl.03 <- merge(tbl.01, tbl.02, all.y=TRUE)
   #rownames(tbl.03) <- tbl.03$gene
   tbl.03$beta[is.na(tbl.03$beta)] <- 0
   tbl.03$IncNodePurity[is.na(tbl.03$IncNodePurity)] <- 0

   fpStarts.list <- lapply(tbl.02$gene, function(gene) subset(tbl.fp, tf_name==gene)[, c("tf_name", "mfpstart")])
   tbl.fpStarts <-  unique(do.call('rbind', fpStarts.list))

   tbl.04 <- merge(tbl.03, tbl.fpStarts, by.x="gene", by.y="tf_name")
   tbl.04 <- tbl.04[order(abs(tbl.04$gene.cor), decreasing=TRUE),]
   
   # footprint.start <- unlist(lapply(rownames(tbl.04), function(gene) subset(tbl.fp, tf_name==gene)$mfpstart[1]))
   gene.info <- subset(tbl.tss, gene_name==target.gene)[1,]
   if(gene.info$strand  == "+"){
      gene.start <- gene.info$start
      tbl.04$distance <- gene.start - tbl.04$mfpstart
   }else{
      gene.start <- gene.info$end
      tbl.04$distance <-  tbl.04$mfpstart - gene.start
      }

   tbl.04

} # createModel
#------------------------------------------------------------------------------------------------------------------------
test.createModel <- function()
{
   printf("--- test.createModel")
   tbl <- createModel("TREM2", promoter.shoulder=100,
                      mtx.expression=mtx.rosmap,
                      absolute.lasso.beta.min=0.0,
                      randomForest.purity.min=1,
                      absolute.expression.correlation.min=0.1)

   checkEquals(ncol(tbl), 6)
   checkEquals(colnames(tbl), c("gene", "gene.cor", "beta", "IncNodePurity", "mfpstart", "distance"))
   checkTrue(nrow(tbl) >= 40)
   checkTrue(all(c("ELF4", "FLI1", "CEBPA", "ELK3") %in% tbl$gene))

     # eliminate thresholds, ensure that more tfs are returned, should be more than 100
   tbl.2 <- createModel("TREM2", promoter.shoulder=100,
                      mtx.expression=mtx.rosmap,
                      absolute.lasso.beta.min=0.0,
                      randomForest.purity.min=0,
                      absolute.expression.correlation.min=0.)
   checkEquals(ncol(tbl.2), 6)
   checkTrue(nrow(tbl.2) > 100)
    
} # test.createModel
#------------------------------------------------------------------------------------------------------------------------
tableToFullGraph <- function(tbl.list)
{
   g <- graphNEL(edgemode = "directed")
   nodeDataDefaults(g, attr = "type") <- "undefined"
   nodeDataDefaults(g, attr = "label") <- "default node label"
   nodeDataDefaults(g, attr = "distance") <- 0
   nodeDataDefaults(g, attr = "gene.cor") <- 0
   nodeDataDefaults(g, attr = "beta") <- 0
   nodeDataDefaults(g, attr = "purity") <- 0

   edgeDataDefaults(g, attr = "edgeType") <- "undefined"
   edgeDataDefaults(g, attr = "beta") <- 0
   edgeDataDefaults(g, attr = "purity") <- 0

   for(target.gene in names(tbl.list)){
      tbl <- tbl.list[[target.gene]]
      tfs <- tbl$gene
      footprints <- unlist(lapply(tbl$distance, function(x)
         if(x < 0)
            sprintf("%s.fp.downstream.%05d", target.gene, abs(x))
         else
            sprintf("%s.fp.upstream.%05d", target.gene, x)))
   
      tbl$footprint <- footprints
      all.nodes <- unique(c(target.gene, tfs, footprints))
      new.nodes <- setdiff(all.nodes, nodes(g))
      g <- addNode(new.nodes, g)
   
      nodeData(g, target.gene, "type") <- "targetGene"
      nodeData(g, tfs, "type")         <- "TF"
      nodeData(g, footprints, "type")  <- "footprint"
      nodeData(g, all.nodes, "label")  <- all.nodes
      nodeData(g, footprints, "label") <- tbl$distance
      nodeData(g, footprints, "distance") <- tbl$distance
   
      nodeData(g, tfs, "gene.cor") <- tbl$gene.cor
      nodeData(g, tfs, "beta") <- tbl$gene.cor
      nodeData(g, tfs, "purity") <- tbl$IncNodePurity
   
      g <- graph::addEdge(tbl$gene, tbl$footprint, g)
      edgeData(g, tbl$gene, tbl$footprint, "edgeType") <- "bindsTo"
      
      g <- graph::addEdge(tbl$footprint, target.gene, g)
      edgeData(g, tbl$footprint, target.gene, "edgeType") <- "regulatorySiteFor"
      } # for target.gene
   
   g

} # tableToFullGraph
#------------------------------------------------------------------------------------------------------------------------
test.tableToFullGraph <- function()
{
   printf("--- test.tableToFullGraph")

       #---- first, just one target.gene, one tbl
   genes <- c("STAT4", "STAT4", "STAT4", "STAT4", "TBR1", "HLF")
   gene.cors <- c(0.9101553, 0.9101553, 0.9101553, 0.9101553, 0.8947867, 0.8872238)
   betas <- c(0.1255503, 0.1255503, 0.1255503, 0.1255503, 0.1448829, 0.0000000)
   IncNodePurities <- c(35.77897, 35.77897, 35.77897, 35.77897, 23.16562, 19.59660)
   mfpstarts <- c(5627705, 5628679, 5629563, 5629581, 5629563, 5629100)
   distances <- c(-2995, -2021, -1137, -1119, -1137, -1600)

   tbl.1 <- data.frame(gene=genes, gene.cor=gene.cors, beta=betas, IncNodePurity=IncNodePurities, 
                       mfpstart=mfpstarts, distance=distances, stringsAsFactors=FALSE)
   
   target.gene.1 <- "EPB41L3"
   tbl.list <- list(tbl.1)
   names(tbl.list) <- target.gene.1
   g1 <- tableToFullGraph(tbl.list)

   checkTrue(all(c(genes, target.gene.1) %in% nodes(g1)))

   fp.nodes <- sort(grep("^fp", nodes(g1), v=TRUE))
   fp.nodes <- sub("fp.downstream.", "", fp.nodes, fixed=TRUE)
   fp.nodes <- as.integer(sub("fp.upstream.", "", fp.nodes, fixed=TRUE))
   checkTrue(all(fp.nodes %in% abs(distances)))

     # one edge from footprint to TF for every footprint
     # one edge from every unique footprint to the target.gene
   expectedEdgeCount <- length(tbl.1$distance) + length(unique(tbl.1$distance))
   checkEquals(length(edgeNames(g1)), expectedEdgeCount)
   checkEquals(length(nodes(g1)),
               1 + length(unique(tbl.1$gene)) + length(unique(tbl.1$distance)))

       #---- now, try a second target.gene and its tbl
   genes <- c("STAT4", "STAT4", "HLF", "TFEB", "HMG20B", "HMG20B")
   gene.cors <- c(0.9162040, 0.9162040, 0.9028613, -0.8256594, -0.8226802, -0.8226802)
   betas <- c(0.1988552, 0.1988552, 0.1492857, 0.0000000, 0.0000000, 0.0000000)
   IncNodePurities <- c(43.37913, 43.37913, 36.93302, 11.27543, 10.17957, 10.17957)
   mfpstarts <- c(4451334, 4451640, 4452647, 4453748, 4453489, 4453491)
   distances <- c(-4001, -3695, -2688, -1587, -1846, -1844)

   tbl.2 <- data.frame(gene=genes, gene.cor=gene.cors, beta=betas, IncNodePurity=IncNodePurities, 
                       mfpstart=mfpstarts, distance=distances, stringsAsFactors=FALSE)
   target.gene.2 <- "DLGAP1"
   
   tbl.list <- list(tbl.2)
   names(tbl.list) <- target.gene.2

   g2 <- tableToFullGraph(tbl.list)

   checkTrue(all(c(genes, target.gene.2) %in% nodes(g2)))

   fp.nodes <- sort(grep("^fp", nodes(g2), v=TRUE))
   fp.nodes <- sub("fp.downstream.", "", fp.nodes, fixed=TRUE)
   fp.nodes <- as.integer(sub("fp.upstream.", "", fp.nodes, fixed=TRUE))
   checkTrue(all(fp.nodes %in% abs(distances)))

     # one edge from footprint to TF for every footprint
     # one edge from every unique footprint to the target.gene
   expectedEdgeCount <- length(tbl.2$distance) + length(unique(tbl.2$distance))
   checkEquals(length(edgeNames(g2)), expectedEdgeCount)
   checkEquals(length(nodes(g2)),
               1 + length(unique(tbl.2$gene)) + length(unique(tbl.2$distance)))

     # now try both tables and target genes together
   
   tbl.list <- list(tbl.1, tbl.2)
   names(tbl.list) <- c(target.gene.1, target.gene.2)
   g3 <- tableToFullGraph(tbl.list)

} # test.tableToFullGraph
#------------------------------------------------------------------------------------------------------------------------
tableToReducedGraph <- function(tbl.list)
{
   g <- graphNEL(edgemode = "directed")
   nodeDataDefaults(g, attr = "type") <- "undefined"
   nodeDataDefaults(g, attr = "label") <- "default node label"
   nodeDataDefaults(g, attr = "degree") <- 0

   edgeDataDefaults(g, attr = "edgeType") <- "undefined"
   edgeDataDefaults(g, attr = "geneCor") <- 0
   edgeDataDefaults(g, attr = "beta") <- 0
   edgeDataDefaults(g, attr = "purity") <- 0
   edgeDataDefaults(g, attr = "fpCount") <- 0

   for(target.gene in names(tbl.list)){
      tbl <- tbl.list[[target.gene]]
      tbl.fpCounts <- as.data.frame(table(tbl$gene))
      tbl <- merge(tbl, tbl.fpCounts, by.x="gene", by.y="Var1")
      colnames(tbl)[grep("Freq", colnames(tbl))] <- "fpCount"
      dups <- which(duplicated(tbl$gene))
      if(length(dups) > 0)
         tbl <- tbl[-dups,]

      tfs <- unique(tbl$gene)
      all.nodes <- unique(c(target.gene, tfs))
      new.nodes <- setdiff(all.nodes, nodes(g))
      g <- addNode(new.nodes, g)
   
      nodeData(g, target.gene, "type") <- "targetGene"
      nodeData(g, tfs, "type")         <- "TF"
      nodeData(g, all.nodes, "label")  <- all.nodes
   
      g <- graph::addEdge(tbl$gene, target.gene, g)
      edgeData(g, tbl$gene, target.gene, "edgeType") <- "regulatorySiteFor"
      edgeData(g, tbl$gene, target.gene, "fpCount") <- tbl$fpCount
      edgeData(g, tbl$gene, target.gene, "geneCor") <- tbl$gene.cor
      edgeData(g, tbl$gene, target.gene, "purity") <- tbl$IncNodePurity
      edgeData(g, tbl$gene, target.gene, "beta") <- tbl$beta
      } # for target.gene
   
   node.degrees <- degree(g)
   degree <- node.degrees$inDegree + node.degrees$outDegree
   nodeData(g, names(degree), attr="degree") <- as.integer(degree)
   g

} # tableToReducedGraph
#------------------------------------------------------------------------------------------------------------------------
test.tableToReducedGraph <- function()
{
   genes <- c("STAT4", "STAT4", "STAT4", "STAT4", "TBR1", "HLF")
   gene.cors <- c(0.9101553, 0.9101553, 0.9101553, 0.9101553, 0.8947867, 0.8872238)
   betas <- c(0.1255503, 0.1255503, 0.1255503, 0.1255503, 0.1448829, 0.0000000)
   IncNodePurities <- c(35.77897, 35.77897, 35.77897, 35.77897, 23.16562, 19.59660)
   mfpstarts <- c(5627705, 5628679, 5629563, 5629581, 5629563, 5629100)
   distances <- c(-2995, -2021, -1137, -1119, -1137, -1600)

   tbl.1 <- data.frame(gene=genes, gene.cor=gene.cors, beta=betas, IncNodePurity=IncNodePurities, 
                       mfpstart=mfpstarts, distance=distances, stringsAsFactors=FALSE)
   
   target.gene.1 <- "EPB41L3"
   tbl.list <- list(tbl.1)
   names(tbl.list) <- target.gene.1
   g1 <- tableToReducedGraph(tbl.list)

   checkTrue(all(c(genes, target.gene.1) %in% nodes(g1)))
   checkEquals(sort(edgeNames(g1)), c("HLF~EPB41L3", "STAT4~EPB41L3", "TBR1~EPB41L3"))

   checkEquals(as.integer(edgeData(g1, from="STAT4", to="EPB41L3", attr="fpCount")), 4)
   checkEquals(as.integer(edgeData(g1, from="HLF", to="EPB41L3", attr="fpCount")), 1)
   checkEquals(as.integer(edgeData(g1, from="TBR1", to="EPB41L3", attr="fpCount")), 1)

   checkEqualsNumeric(unlist(edgeData(g1, attr="purity"), use.names=FALSE), c(19.59660, 35.77897, 23.16562))
   checkEqualsNumeric(unlist(edgeData(g1, attr="geneCor"), use.names=FALSE), c(0.8872238, 0.9101553, 0.8947867))
   checkEqualsNumeric(unlist(edgeData(g1, attr="beta"), use.names=FALSE), c(0.0000000, 0.1255503, 0.1448829))

} # test.tableToReducedGraph
#------------------------------------------------------------------------------------------------------------------------
renderAsNetwork <- function(tbl, target.gene)
{
   tbl.list <- list(tbl)
   names(tbl.list) <- target.gene
   g <- tableToFullGraph(tbl.list)
   rcy <- RCyjs(10000:10100, title=target.gene, graph=g)
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
   fpXcoords <- -1 * as.integer(xPos)    # negative positions are downstream of TSS, traditionally to the right
   setPosition(rcy, data.frame(id=names(xPos), x=fpXcoords, y=0, stringsAsFactors=FALSE))
   target.gene <- names(which(nodeData(g, attr="type") == "targetGene"))
   setPosition(rcy, data.frame(id=target.gene, x=0, y=-200, stringsAsFactors=FALSE))

   fpTFedges <- inEdges(g)[fp.nodes]

   for(fp in names(fpTFedges)){
     tfs <- fpTFedges[[fp]]
     xpos.base <- -1 * as.integer(nodeData(g, fp, "distance"))
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
hAlign <- function(rcy)
{
  selected.node.ids <- getSelectedNodes(rcy)$id
  if(length(selected.node.ids) == 0){
      printf("no nodes selected, nothing to hAlign")
      return()
      }
  tbl.pos <- getPosition(rcy, selected.node.ids)
  tbl.pos$y <- mean(tbl.pos$y)
  setPosition(rcy, tbl.pos)

} # hAlign
#------------------------------------------------------------------------------------------------------------------------
vAlign <- function(rcy)
{
  selected.node.ids <- getSelectedNodes(rcy)$id
  if(length(selected.node.ids) == 0){
     printf("no nodes selected, nothing to vAlign")
     return();
     }

  tbl.pos <- getPosition(rcy, selected.node.ids)
  tbl.pos$x <- mean(tbl.pos$x)
  setPosition(rcy, tbl.pos)

} # vAlign
#------------------------------------------------------------------------------------------------------------------------
demo <- function()
{
   target.gene <- "APOE"
   target.gene <- "MBP"
   target.gene <- "MOBP"
   target.gene <- "TREM2"
   target.gene <- "SCN2A"
   target.gene <- "TYROBP"

   tbl <- createModel(target.gene, promoter.shoulder=10000,
                      mtx.expression=mtx.rosmap,
                      absolute.lasso.beta.min=0.1,
                      randomForest.purity.min=1,
                      absolute.expression.correlation.min=0.05)

   rcy <- renderAsNetwork(tbl, target.gene)
   layoutByFootprintPosition(rcy)

} # demo
#------------------------------------------------------------------------------------------------------------------------
demo.reducedGraph <- function()
{
   target.genes <- c("DLGAP1", "EPB41L3")

   tbl.1 <- createModel(target.genes[1], promoter.shoulder=5000,
                        mtx.expression=mtx.rosmap,
                        absolute.lasso.beta.min=0.1,
                        randomForest.purity.min=1,
                        absolute.expression.correlation.min=0.05)

   tbl.2 <- createModel(target.genes[2], promoter.shoulder=5000,
                        mtx.expression=mtx.rosmap,
                        absolute.lasso.beta.min=0.1,
                        randomForest.purity.min=1,
                        absolute.expression.correlation.min=0.05)

   tbl.list <- list(tbl.1, tbl.2)
   names(tbl.list) <- target.genes
   g <- tableToReducedGraph(tbl.list)

   rcy <- RCyjs(10000:10100, graph=g)
   httpSetStyle(rcy, "tyle-reducedGraph.js")
   selectNodes(rcy, names(which(nodeData(g, attr="degree") == 1)))
   hideSelectedNodes(rcy)
   layout(rcy, "cose")

   rcy

} # demo.reducedGraph
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()

