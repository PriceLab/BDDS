# makefile to run all tests for wellington functionality
BAMFILE = data/ENCSR000DBY.19.chr
SUBDIRS = bamToBed fseq transform wellington hint
SUBCLEAN = $(addsuffix .clean,$(SUBDIRS))

## all        : run all tests - if PYTHON=3 argument given, runs python3 version of wellington
.PHONY: all $(SUBDIRS)
all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

## bamToBed   : bamToBed test
bamToBed: $(BAMFILE).bam

## fseq       : fseq test
fseq: bamToBed

## transform  : transform test
transform: fseq

## wellington : wellington test - if PYTHON=3 argument given, runs python3 version
wellington: transform

## hint       : hint test
hint: transform

## clean      : clean up intermediate files from testing
.PHONY: clean $(SUBCLEAN)
clean: $(SUBCLEAN)

$(SUBCLEAN): %.clean:
	$(MAKE) -C $* -f makefile clean
	
.PHONY: help
help: makefile
	@sed -n 's/^##//p' $<
