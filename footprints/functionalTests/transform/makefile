# a test for fseq functionality
BEDFILE=../data/chr19.bed
OUTFILE=../data/chr19.peaks

default: preTest $(OUTFILE).bed postTest

.PHONY: clean
clean:
	@rm $(OUTFILE).std.bed
	@rm $(OUTFILE).dif.bed
	@rm $(OUTFILE).400.bed

.PHONY: preTest
preTest:
	python preTest.py

$(OUTFILE).bed: $(BEDFILE)
	@cut -f 1,2,3 $(BEDFILE) > $(OUTFILE).std.bed
	@awk '{a=$$2-$$3;print $$0 "\t" a;}' $(OUTFILE).std.bed \
	> $(OUTFILE).dif.bed
	@cat $(OUTFILE).dif.bed | awk '($$4 < -400)' | cut -f 1,2,3 \
	> $(OUTFILE).400.bed

.PHONY: postTest
postTest:
	python postTest.py $(OUTFILE).400.bed
