source("../renderGeneModel.R")

target.gene <- "TYROBP"
target.gene <- "TREM2"

tbl <- createModel(target.gene, promoter.shoulder=1000,
                   mtx.expression=mtx.rosmap,
                   absolute.lasso.beta.min=0.1,
                   randomForest.purity.min=1,
                   absolute.expression.correlation.min=0.2)

rcy <- renderAsNetwork(tbl, target.gene)
layoutByFootprintPosition(rcy)
httpSetStyle(rcy, "style-purityControlsNodeSize.js")

