# a test for fseq functionality
BEDFILE=../data/ENCSR000DBY.19

default: preTest $(BEDFILE).bed postTest

clean:
	- rm $(BEDFILE).peaks.bed

preTest:
	python preTest.py

$(BEDFILE).peaks.bed: $(BEDFILE).bed
	fseq -f 0 -o ../data pt 2.5 -of bed $(BEDFILE).peaks.bed

postTest:
	python postTest.py $(BEDFILE).peaks.bed
