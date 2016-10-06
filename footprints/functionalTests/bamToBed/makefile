BAMFILE=../data/ENCSR000DBY.19.chr

default: preTest $(BAMFILE).bed postTest

clean:
	- rm $(BAMFILE).bam.bai
	- rm $(BAMFILE).bed

preTest:
	python3 preTest.py

$(BAMFILE).bed: $(BAMFILE).bam $(BAMFILE).bam.bai
	bamToBed -i $(BAMFILE).bam > $(BAMFILE).bed

$(BAMFILE).bam.bai: $(BAMFILE).bam
	samtools index $(BAMFILE).bam

postTest:
	python3 postTest.py $(BAMFILE).bed
