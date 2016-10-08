import sys, os
from pandas import *
#--------------------------------------------------------------------------------
if(len(sys.argv) !=2):
    print("usage: postTest.py <file.bed>")
    exit(1)

filename = sys.argv[1]
assert(os.path.isfile(filename))
tbl = read_csv(filename, sep="\t")
assert((tbl.shape[0] > 17000) and tbl.shape[0] < 19000 )
#print("table shape: " + str(tbl.shape))
assert(min(tbl.iloc[:, 1].tolist()) > 64000)
assert(max(tbl.iloc[:, 1].tolist()) > 58606000)
chromosomes = unique(tbl.iloc[:5, 0]).tolist()
assert(len(chromosomes) == 1)
assert(chromosomes[0] == "chr19")
