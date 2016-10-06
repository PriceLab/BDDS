import sys, os
from pandas import *
#--------------------------------------------------------------------------------
if(len(sys.argv) !=2):
    print("usage: postTest.py <file.bed>")
    exit(1)

filename = sys.argv[1]
assert(os.path.isfile(filename))
tbl = read_csv(filename, sep="\t")
#assert(tbl.shape == (1279016, 6))
assert(tbl.shape == (1108415,6))
assert(min(tbl.iloc[:, 1].tolist()) == 60015)
assert(max(tbl.iloc[:, 1].tolist()) == 58607530)
chromosomes = unique(tbl.iloc[:5, 0]).tolist()
assert(len(chromosomes) == 1)
assert(chromosomes[0] == 'chr19')
