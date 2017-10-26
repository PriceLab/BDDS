# This script houses some useful shared functions created for analysis purposes
#----------------------------------------------------------------------------------------------------
# Function that reads the FIMO data for one chromosome
getRegionsForChrom <- function(chromosome){

    fimo.conn <- dbConnect(PostgreSQL(),
                           user = "trena",
                           password = "trena",
                           port = "5432",
                           host = "khaleesi",
                           dbname = "fimo"
                           )

    # Print out what chromosome we're using
    to.print <- sprintf("Working on %s", chromosome)
    print(to.print)

    # Make a string for the chromosome
    query <- sprintf("select chrom, start, endpos from fimo_hg38 where chrom = '%s';",
                     chromosome)
    fimo.data <- dbGetQuery(fimo.conn, query)
    dbDisconnect(fimo.conn)

    # Print out the number of records retrieved
    to.print <- sprintf("%d records retrieved for %s", nrow(fimo.data), chromosome)
    return(fimo.data)
    
} #getRegionsForChrom
#----------------------------------------------------------------------------------------------------
# Function for counting coverage given a genomic range
calculateGenomeCoverage <- function(my.range){

    # Use the "reduce" function to calculate the total areas covered
    reduced.range <- GenomicRanges::reduce(my.range)

    # Sum the areas using the "width" operator to get each area
    total.bp <- sum(as.numeric(width(reduced.range)))

    return(total.bp)
    
} #calculateGenomeCoverage
#----------------------------------------------------------------------------------------------------
# Function that takes tissue name and figures out the pathway to its footprints on Khaleesi
compileFootprintPaths <- function(tissue){

    # Create the 4 strings from the tissue
    hint.20 <- paste0("/ssd/mrichard/data/footprints/",tissue,"_hint_20")
    hint.16 <- paste0("/ssd/mrichard/data/footprints/",tissue,"_hint_16")
    well.20 <- paste0("/ssd/mrichard/data/footprints/",tissue,"_wellington_20")
    well.16 <- paste0("/ssd/mrichard/data/footprints/",tissue,"_wellington_16")

    # Create all the paths
    hint.20.paths <- list.files(hint.20)
    hint.20.paths <- file.path(hint.20, hint.20.paths)

    hint.16.paths <- list.files(hint.16)
    hint.16.paths <- file.path(hint.16, hint.16.paths)

    well.20.paths <- list.files(well.20)
    well.20.paths <- file.path(well.20, well.20.paths)

    well.16.paths <- list.files(well.16)
    well.16.paths <- file.path(well.16, well.16.paths)

    # Stick the paths together and return them
    all.paths <- c(hint.20.paths,
                   hint.16.paths,
                   well.20.paths,
                   well.16.paths
                   )

    return(all.paths)
    
} #compileFootprintPaths
#----------------------------------------------------------------------------------------------------
# Function that takes pathways and uses them to read FPs
readFootprintTable <- function(fp.path){
    tbl <- read.table(fp.path, sep = "\t", as.is = TRUE)
    colnames(tbl) <- c("chrom", "start", "end", "name", "score", "strand")
    # Make sure it's just the first 6 columns
    tbl <- tbl[,c("chrom", "start", "end", "name", "score", "strand")]
    return(tbl)
    
} #readFootprintTable
#----------------------------------------------------------------------------------------------------
# Define a function that takes a list of footprint dataframes and makes it into a genomic range
footprintsToGenomicRange <- function(fp.list){

    # Remove the non-start/end/stop columns
    to.keep <- c("chrom","start","end")
    fp.list <- lapply(fp.list,dplyr::select, one_of(to.keep))

    # Collapse into one data frame
    fp.df <- purrr::reduce(fp.list, dplyr::union)

    # Convert it into a genomic range
    my.gr <- with(fp.df, GRanges(seqnames=chrom, IRanges(start=start, end=end)))

    return(my.gr)
} #footprintsToGenomicRange
#----------------------------------------------------------------------------------------------------
# Build a function that takes the tissue and goes all the way through to coverage using
# previously defined functions
coverageFromTissue <- function(tissue){

    # First find the footprints
    fp <- compileFootprintPaths(tissue)

    # Then use a parallel process to read all the footprints
    my.message <- sprintf("Reading %d files for %s", length(fp), tissue)
    print(my.message)
    register(MulticoreParam(workers = 30, stop.on.error = TRUE, log = FALSE), default = TRUE)
    fp.list <- bptry(bplapply(fp, readFootprintTable))

    # Turn the footprint list into a genomic range
    fp.gr <- footprintsToGenomicRange(fp.list)
    rm(fp.list)

    # Find the coverage
    fp.coverage <- calculateGenomeCoverage(fp.gr)
    return(fp.coverage)
    
} #coverageFromTissue
#----------------------------------------------------------------------------------------------------
# Function to turn tissue into reduced GR; used for getting a list of them
createFullCoverage <- function(tissue){

    # First find the footprints
    fp <- compileFootprintPaths(tissue)

    # Then use a parallel process to read all the footprints
    my.message <- sprintf("Reading %d files for %s", length(fp), tissue)
    print(my.message)
    register(MulticoreParam(workers = 30, stop.on.error = TRUE, log = FALSE), default = TRUE)
    fp.list <- bptry(bplapply(fp, readFootprintTable))

    # Turn the footprint list into a genomic range
    fp.gr <- footprintsToGenomicRange(fp.list)
    rm(fp.list)

    # Reduce the genomic range
    reduced.gr <- GenomicRanges::reduce(fp.gr)

    # Return a reduced GR
    return(reduced.gr)
    
} #createFullCoverage
#----------------------------------------------------------------------------------------------------
