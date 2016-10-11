# a test for fseq functionality
BEDFILE=../data/ENCSR000DBY.19.chr
OUTFILE=../data/chr19

default: preTest $(OUTFILE).bed postTest

.PHONY: clean
clean:
	@rm $(OUTFILE).bed

.PHONY: preTest
preTest:
	python preTest.py

$(OUTFILE).bed: $(BEDFILE).bed
	@fseq -f 0 -o ../data -t 2.5 -of bed $(BEDFILE).bed

.PHONY: postTest
postTest:
	python postTest.py $(OUTFILE).bed
