# assembleMatrix
#------------------------------------------------------------------------------------------------------------------------
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test.matrixForDirectory()
    
} # runTests
#------------------------------------------------------------------------------------------------------------------------
matrixForDirectory <- function(dir, testing=FALSE)
{
   full.dir <- file.path(dir, "vec")
   files <- list.files(full.dir, pattern=".*vec$")
   file.max <- length(files)
   cols <- list()

   if(testing)
      file.max <- 5

   for(file in files[1:file.max]){
      full.path <- file.path(full.dir, file)
      vector <- scan(full.path, sep="\n", what=numeric(), quiet=TRUE)
      if(testing) vector <- vector[1:8]
      colname <- sub(".vec", "", file, fixed=TRUE)
      cols[[colname]] <- vector
      } # for file

   mtx <- do.call("cbind", cols)
   invisible(mtx)

} # matrixForDirectory
#------------------------------------------------------------------------------------------------------------------------
test.matrixForDirectory <- function()
{
   printf("--- matrixForDirectory")
   mtx <- matrixForDirectory("2046", testing=TRUE)
   checkEquals(dim(mtx), c(8, 5))
   checkTrue(typeof(mtx) == "double")

   mtx <- matrixForDirectory("2046")
   checkEquals(dim(mtx), c(32061, 96))
    
} # test.matrixForDirectory
#------------------------------------------------------------------------------------------------------------------------
assembleMatrices <- function()
{
   experiment.directories <- c("2046", "2051", "2078", "2080", "2159", "2167")
   mtx.all <- list()
   for(dir in experiment.directories){
      mtx.single <- matrixForDirectory(dir)
      mtx.all[[dir]] <- mtx.single
      } # for dir

   mtx <- do.call("cbind", mtx.all)
   invisible(mtx)

   tbl.ids <- read.table("features.tab", sep="\t", header=FALSE, as.is=TRUE)
   colnames(tbl.ids) <- c("gene", "variant", "mystery")
   tokens <- strsplit(tbl.ids$gene, ":")
   genes <-  unlist(lapply(tokens, "[", 2))
   tbl.ids$gene <- genes
   save(tbl.ids, file="tbl.ids.withGeneNamesSomeDuplicated.RData")
   dups <- which(duplicated(tbl.ids$gene))
   stopifnot(length(dups) > 0)
   unique.genes <- tbl.ids$gene[-dups]
   mtx <- mtx[-dups,]
   rownames(mtx) <- unique.genes
   save(mtx, file="mtxSkin17858x576.RData")
   
} # assembleMatrices
#------------------------------------------------------------------------------------------------------------------------
spotCheckFinalMatrix <- function()
{
   load("tbl.ids.withGeneNamesSomeDuplicated.RData")
   load("mtxSkin17858x576.RData")

   checkEquals(dim(mtx), c(17858, 576))

   for(i in 1:50){
      gene <- sample(rownames(mtx),1)
      sampleName <- sample(colnames(mtx), 1)
      directory <- sub("^GSS_*", "", sampleName)
      directory <- sub("_.*$", "", directory)
      filename <- file.path(directory, "vec", sprintf("%s.vec", sampleName))
      checkTrue(file.exists(filename))
      vector <- scan(filename, what=numeric(), sep="\n", quiet=TRUE)
         # which element/s in the vector correspond to our gene?
      search.term <- sprintf("^%s$", gene)
      hits <- grep(search.term, tbl.ids$gene)
      value.in.matrix <- mtx[gene, sampleName]
      printf("validating mtx[%s, %s]: %f", gene, sampleName, value.in.matrix)
      checkTrue(value.in.matrix %in% vector[hits])
      }

} # spotCheckFinalMatrix
#------------------------------------------------------------------------------------------------------------------------
