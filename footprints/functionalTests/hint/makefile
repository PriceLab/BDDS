# a test for hint functionality
#
# hint requires python 2, but the tests take python3...
HOST=$(shell hostname)
ifeq (${HOST},whovian.systemsbiology.net)
	export PYTHONPATH := $(PYTHONPATH):/local/tools/HINT-1.1.1/lib/python2.7/site-packages/
	export PATH := /usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/dell/srvadmin/bin:/local/tools/bin
endif
SAMPLE=ENCSR000DBY
BAMFILE=../data/ENCSR000DBY.19.chr.bam
INDEX=../data/ENCSR000DBY.19.chr.bam.bai
BEDFILE=../data/chr19.peaks.400.bed
INPUTMATRIX=./inputmatrix.txt
OUTDIR=../output/hint/
OUTFILE=$(SAMPLE).bed

#default: preTest $(INPUTMATRIX) $(OUTFILE) postTest
default: preTest $(INPUTMATRIX) $(OUTFILE)

.PHONY: clean
clean:
	@rm $(INPUTMATRIX)
	@rm $(OUTDIR)*

.PHONY: preTest
preTest:
ifeq (${HOST},whovian.systemsbiology.net)
	python3 preTest.py
else
	python preTest.py
endif

$(INPUTMATRIX): $(BAMFILE)
	@touch $(INPUTMATRIX)
	@echo -e "name\ttype\tfile\tdata\tgroup" >> $(INPUTMATRIX)
	@echo -e "HS2\tregions\t$(BEDFILE)\tHS\t$(SAMPLE)" >> $(INPUTMATRIX)
	@echo -e "DNase\treads\t$(BAMFILE)\tDNASE\t$(SAMPLE)" >> $(INPUTMATRIX)

$(OUTFILE): $(BEDFILE) $(BAMFILE) $(INDEX) $(INPUTMATRIX)
	rgt-hint --output-location $(OUTDIR) $(INPUTMATRIX)

#.PHONY: postTest
#postTest:
#ifeq (${HOST},whovian.systemsbiology.net)
#	python3 postTest.py $(OUTDIR)$(OUTFILE)
#else
#	python postTest.py $(OUTDIR)$(OUTFILE)
#endif
