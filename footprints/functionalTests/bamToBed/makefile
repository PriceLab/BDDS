BAMFILE=../data/ENCSR000DBY.19.chr

default: preTest $(BAMFILE).bed postTest

.PHONY: clean
clean:
	- rm $(BAMFILE).bam.bai
	- rm $(BAMFILE).bed

.PHONY: preTest
preTest:
	python3 preTest.py

$(BAMFILE).bed: $(BAMFILE).bam $(BAMFILE).bam.bai
	bamToBed -i $(BAMFILE).bam > $(BAMFILE).bed

$(BAMFILE).bam.bai: $(BAMFILE).bam
	samtools index $(BAMFILE).bam

.PHONY: postTest
postTest:
	python3 postTest.py $(BAMFILE).bed
