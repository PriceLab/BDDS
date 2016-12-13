library(gplots)
stopifnot(packageVersion("TReNA") >= '0.99.57')
source("../renderGeneModel.R")
#------------------------------------------------------------------------------------------------------------------------
probeHeatmapOfSelectedGene <- function()
{
   tbl.selection <- getSelectedNodes(rcy)
   stopifnot(nrow(tbl.selection) == 1)

   geneOfInterest <- tbl.selection$name[1]
   geneNameAdaptedForProbeSearch <- sprintf(":%s.", geneOfInterest)
   probeRows <- grep(geneNameAdaptedForProbeSearch, rownames(mtx.probes))
   stopifnot(length(probeRows) > 1)
   heatmap.2(mtx.probes[probeRows,], trace='none', col=rev(heat.colors(10)), margins=c(15,15), cexCol=0.1, cexRow=1)
       
    
} # probeHeatmapOfSelectedGene
#------------------------------------------------------------------------------------------------------------------------
target.gene <- "COL1A1"
mtx.filename <- "mtxSkin17858x576.RData"
print(load(mtx.filename))
mtx.skin <- mtx
mtx.byProbe.filename <- "mtxSkinByProbe.RData"
print(load(mtx.byProbe.filename))

tbl <- createModel(target.gene, promoter.shoulder=5000,
                   mtx.expression=mtx.skin,
                   absolute.lasso.beta.min=0.01,
                   randomForest.purity.min=1,
                   absolute.expression.correlation.min=0.2)

rcy <- renderAsNetwork(tbl, target.gene)
layoutByFootprintPosition(rcy)
httpSetStyle(rcy, "../coryADpresentationNov2016/style-purityControlsNodeSize.js")

