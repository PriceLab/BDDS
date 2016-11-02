# a test for hint functionality
HOST=$(shell hostname)
ifeq (${HOST},whovian.systemsbiology.net)
	export PYTHONPATH := $(PYTHONPATH):/local/tools/HINT-1.1.1/lib/python2.7/site-packages/

endif
SAMPLE=ENCSR000DBY
BAMFILE=../data/ENCSR000DBY.19.chr.bam
INDEX=../data/ENCSR000DBY.19.chr.bam.bai
BEDFILE=../data/chr19.peaks.400.bed
INPUTMATRIX=./inputmatrix.txt
OUTDIR=../output/hint/
#OUTFILE=FOO

default: preTest $(INPUTMATRIX) $(OUTFILE) postTest

.PHONY: clean
clean:
	@rm $(INPUTMATRIX)
	@rm $(OUTDIR)*

.PHONY: preTest
preTest:
	python preTest.py

$(INPUTMATRIX): $(BAMFILE)
	@touch $(INPUTMATRIX)
	@echo -e "name\ttype\tfile\tdata\tgroup" >> $(INPUTMATRIX)
	@echo -e "HS2\tregions\t$(BEDFILE)\tHS\t$(SAMPLE)" >> $(INPUTMATRIX)
	@echo -e "DNase\treads\t$(BAMFILE)\tDNASE\t$(SAMPLE)" >> $(INPUTMATRIX)

#$(OUTFILE): $(BEDFILE) $(BAMFILE) $(INDEX) $(INPUTMATRIX)
    
#.PHONY: postTest
#postTest:
#	python postTest.py $(OUTFILE)
