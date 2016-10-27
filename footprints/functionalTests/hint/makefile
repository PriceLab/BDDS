# a test for hint functionality
HOST=$(shell hostname)
ifeq(${HOST},whovian.systemsbiology.net)
	export PYTHONPATH := $(PYTHONPATH):/local/tools/HINT-1.1.1/lib/python2.7/site-packages/

endif

BAMFILE=../data/ENCSR000DBY.19.chr.bam
BEDFILE=../data/chr19.peaks.400.bed
OUTFILE=../output/wellington/ENCSR000DBY.19.chr.bam.chr19.peaks.400.bed.WellingtonFootprints.FDR.0.01.bed
OUTDIR=../output/wellington/

default: preTest $(OUTFILE) postTest

.PHONY: clean
clean:
	@find $(OUTDIR)*.wig -delete 
	@find $(OUTDIR)p* -delete

.PHONY: preTest
preTest:
	python preTest.py

$(OUTFILE): $(BEDFILE) $(BAMFILE)
	@find $(OUTDIR)* -delete
ifeq ($(PYTHON),3)
	wellington_footprints3.py -fdrlimit=0 $(BEDFILE) $(BAMFILE) $(OUTDIR)
else
	wellington_footprints.py -fdrlimit=0 $(BEDFILE) $(BAMFILE) $(OUTDIR)
endif

.PHONY: postTest
postTest:
	python postTest.py $(OUTFILE)
