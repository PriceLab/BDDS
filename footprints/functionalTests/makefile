# makefile to run all tests for wellington functionality
BAMFILE = data/ENCSR000DBY.19.chr
SUBDIRS = bamToBed fseq transform
SUBCLEAN = $(addsuffix .clean,$(SUBDIRS))

## all        : run all tests
.PHONY: all $(SUBDIRS)
all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

bamToBed: $(BAMFILE).bam

fseq: bamToBed

transform: fseq

## clean      : clean up intermediate files from testing
.PHONY: clean $(SUBCLEAN)
clean: $(SUBCLEAN)

$(SUBCLEAN): %.clean:
	$(MAKE) -C $* -f makefile clean
	
.PHONY: help
help: makefile
	@sed -n 's/^##//p' $<
