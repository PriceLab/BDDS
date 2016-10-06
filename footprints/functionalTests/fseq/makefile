# a test for fseq functionality
BEDFILE=../data/ENCSR000DBY.19.chr
OUTFILE=../data/chr19

default: preTest $(OUTFILE).bed postTest

clean:
	- rm $(OUTFILE).bed

preTest:
	python preTest.py

$(OUTFILE).bed: $(BEDFILE).bed
	fseq -f 0 -o ../data -t 2.5 -of bed $(BEDFILE).bed

postTest:
	python postTest.py $(OUTFILE).bed
