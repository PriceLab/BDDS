library(tools)
library(parallel)
cores <- detectCores()

runPiq <- function(motifNumber, bamName) {
    function(motifNumber, bamName) {
        cmd <- sprintf("/home/parsing/piq/parse_output.sh  %d %s", motifNumber, bamName)
        print(cmd)
        system(cmd)
        }
    }

f <- runPiq(motifNumber, bamName)

sprintf("cores = %d", cores)

bam_temp = list.files("/scratch/lympho_bam/mod_bam2")
bam_temp2 = basename(bam_temp)
bam_list = file_path_sans_ext(bam_temp2)

motifNumbers <- 1:519

for (bamName in bam_list) {
        system.time(mclapply(motifNumbers, function(i) f(i, bamName), mc.cores = cores))
}