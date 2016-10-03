BAMFILE=../data/ENCSR000DBY.19

default:  $(BAMFILE).bed
clean:
	- rm $(BAMFILE).bam.bai
	- rm $(BAMFILE).bed

$(BAMFILE).bed: $(BAMFILE).bam $(BAMFILE).bam.bai
	python bamToBed.py $(BAMFILE).bam

$(BAMFILE).bam.bai: $(BAMFILE).bam
	samtools index $(BAMFILE).bam
