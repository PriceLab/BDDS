source("../renderGeneModel.R")

target.gene <- "TYROBP"
target.gene <- "TREM2"
target.gene <- "REST"
target.gene <- "SNAP25"
target.gene <- "VGF"
print(load("mayo_TCX_counts_matrix_normalized_geneSymbols_25031x277.RData"))
mtx.mayoTCX <- mtx
print(load("rosmap_counts_matrix_normalized_geneSymbols_25031x638.RData"))
mtx.rosmapHuge <- mtx
tbl <- createModel(target.gene, promoter.shoulder=5000,
                   mtx.expression=mtx.rosmap,
                   #mtx.expression=mtx.mayoTCX,
                   #mtx.expression=mtx.rosmap,
                   absolute.lasso.beta.min=0.1,
                   randomForest.purity.min=1,
                   absolute.expression.correlation.min=0.2)

rcy <- renderAsNetwork(tbl, target.gene)
layoutByFootprintPosition(rcy)
httpSetStyle(rcy, "style-purityControlsNodeSize.js")

