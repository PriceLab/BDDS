The *functionalTests* directory contains a recursive, self documenting makefile-based set of test and example functions that define our DNase hypersensitivity data processing footprint identification workflow(s).

Documentation is available with `make help`. Specific steps can be run, such as `make hint`. After running, use `make clean` to remove intermediate files that are produced in a test run.

The test processes data for chromosome 19 for a single sample of lymphoblast data.

For example, all tests can be run with `make all PYTHON=3`. In this case, the source .bam file from the data directory is converted to a bed file (`make bamToBed`), fseq is used to identify peaks (`make fseq`), the fseq output is rearranged (`make transform`, and footprints are called with wellington (`make wellington`) and with hint (`make hint`). Intermediate files are written to the data directory, and the output is written to the output directory. If all's well, `make clean` deletes these files.
