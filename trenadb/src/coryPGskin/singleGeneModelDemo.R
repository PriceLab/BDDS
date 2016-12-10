source("../renderGeneModel.R")

target.gene <- "COL1A1"
mtx.filename <- "mtxSkin17858x576.RData"
print(load(mtx.filename))
mtx.skin <- mtx

tbl <- createModel(target.gene, promoter.shoulder=5000,
                   mtx.expression=mtx.skin,
                   #mtx.expression=mtx.mayoTCX,
                   #mtx.expression=mtx.rosmap,
                   absolute.lasso.beta.min=0.1,
                   randomForest.purity.min=1,
                   absolute.expression.correlation.min=0.2)

rcy <- renderAsNetwork(tbl, target.gene)
layoutByFootprintPosition(rcy)
httpSetStyle(rcy, "../coryADpresentationNov2016/style-purityControlsNodeSize.js")

