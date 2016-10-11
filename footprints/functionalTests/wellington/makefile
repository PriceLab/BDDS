# a test for fseq functionality
BAMFILE=../data/ENCSR000DBY.19.chr.bam
BEDFILE=../data/chr19.peaks.400.bed
OUTFILE=../output/wellington/ENCSR000DBY.19.chr.bam.chr19.peaks.400.bed.WellingtonFootprints.FDR.0.01.bed
OUTDIR=../output/wellington/

default: preTest $(OUTFILE) postTest

.PHONY: clean
clean:
	@ rm $(OUTDIR)p\ value\ cutoffs/*
	@ rmdir $(OUTDIR)p\ value\ cutoffs/
	@ rm $(OUTDIR)*

.PHONY: preTest
preTest:
	python preTest.py

$(OUTFILE): $(BEDFILE) $(BAMFILE)
ifeq ($(PYTHON),3)
	wellington_footprints3.py $(BEDFILE) $(BAMFILE) $(OUTDIR)
else
	wellington_footprints.py $(BEDFILE) $(BAMFILE) $(OUTDIR)
endif

.PHONY: postTest
postTest:
	python postTest.py $(OUTFILE)
