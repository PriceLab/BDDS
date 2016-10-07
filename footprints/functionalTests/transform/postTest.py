import sys, os
from pandas import *
#--------------------------------------------------------------------------------
if(len(sys.argv) !=2):
    print("usage: postTest.py <file.bed>")
    exit(1)

filename = sys.argv[1]
assert(os.path.isfile(filename))
tbl = read_csv(filename, sep="\t")
#assert(tbl.shape == (2970, 3))
assert(tbl.shape == (2976, 3))
assert(min(tbl.iloc[:, 1].tolist()) == 69929)
assert(max(tbl.iloc[:, 1].tolist()) == 58577629)
chromosomes = unique(tbl.iloc[:5, 0]).tolist()
assert(len(chromosomes) == 1)
assert(chromosomes[0] == "chr19")
