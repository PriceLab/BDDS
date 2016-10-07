import sys, os
from pandas import *
#--------------------------------------------------------------------------------
if(len(sys.argv) !=2):
    print("usage: postTest.py <file.bed>")
    exit(1)

filename = sys.argv[1]
assert(os.path.isfile(filename))
tbl = read_csv(filename, sep="\t")
assert(tbl.shape[1] == 3)
assert((tbl.shape[0] > 2800) and (tbl.shape[0] < 2900))
assert(min(tbl.iloc[:, 1].tolist()) > 69000)
assert(max(tbl.iloc[:, 1].tolist()) < 58578000)
chromosomes = unique(tbl.iloc[:5, 0]).tolist()
assert(len(chromosomes) == 1)
assert(chromosomes[0] == "chr19")
