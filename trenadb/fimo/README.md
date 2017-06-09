# Instructions for creating/updating fimo

When new motifs become available, they should be added to the existing fimo database. For the purposes of this short tutorial, we'll assume said motifs have already been added to the PriceLab/MotifDb repository. If not, please refer to that repository, located [here](https://github.com/PriceLab/MotifDb). At this point, the workflow for doing so is as follows:

1. Install the latest MotifDb version from the forked PriceLab repo
2. Using MotifDb in R, export .meme files of your desired motifs
3. Using the fimo program, intersect the motif .meme files with all chromosomes
4. Copy the existing fimo database dump from Amazon S3
5. Restore the database (I do this on an Amazon EC2 instance)
6. Using the template from [`create_fimo_table.sh`](https://github.com/PriceLab/BDDS/blob/master/trenadb/fimo/create_fimo_table.sh), create a shell script that copies the fimo output into the existing fimo database
7. Create indices using the commands in [index.sql](https://github.com/PriceLab/BDDS/blob/master/trenadb/fimo/index.sql)
8. Dump the new fimo database locally
9. Copy the database dump to Amazon S3

What follows is a more detailed description of how to carry out all these steps, with examples using the fimo version created on June 8, 2017. 

## 1. Install the latest MotifDb version from the forked PriceLab repo

First, clone the repo onto your machine:

`git clone https://github.com/PriceLab/MotifDb.git`

Next, open up R and install the necessary dependency packages from BioConductor using the following:

`source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocGenerics", "S4Vectors", "IRanges", "Biostrings", "rtracklayer"))
`
Once the dependency packages are installed (which will take a few minutes), exit R and navigate to the MotifDb root directory (e.g. /scratch/github/MotifDb) and install the MotifDb package locally:

`R CMD INSTALL .`

This should install the MotifDb library on your machine. 

## 2. Using MotifDb in R, export .meme files of your desired motifs

Open up R again and load the MotifDb package:

`library(MotifDb)`

As detailed in the [MotifDb vignette](http://bioconductor.org/packages/release/bioc/vignettes/MotifDb/inst/doc/MotifDb.pdf), you can filter out certain types of motifs using the `query` function. For example, we can select the human (*H.sapiens*) motifs from JASPAR 2016 as follows:

`jaspar_human <- query(query(MotifDb, "hsapiens"),"jaspar2016")`

In order to run these motifs through the fimo program, we want a .meme file; this can be directly done using MotifDb as follows:

`export(jaspar_human, con = "./jaspar_human.meme", format = "meme")`

We've now created the file "jaspar_human.meme" of human motifs found in the JASPAR 2016 database. 

## 3. Using the fimo program, intersect the motif .meme files with all chromosomes

Assuming you've installed the fimo tools from the MEME suite (if not, find instructions [here](http://meme-suite.org/doc/install.html?man_type=web), you'll need to run fimo on your new motif .meme file and do so on all chromosomes. Your .meme file is already in your workspace, but you will need the hg38 chromosomes. We have them broken up into sensibly-sized groups in an Amazon S3 bucket, so copy all of them to your workspace with the following:

`aws s3 cp s3://marichards/GRCh38 . --recursive`

There are 17 groups of chromosomes; you can run these individually, but if you are on an AWS instance then it makes sense to run them all in parallel as these will take a while. As a guide, please see the example script created for running all HOMER motifs, [`run_all_homer.sh`]((https://github.com/PriceLab/BDDS/blob/master/trenadb/fimo/run_all_homer.sh)). 

If you insist on doing chromosomes individually, the notation is as follows:

`fimo --text --oc . --no-qvalue ../meme/homer_all.meme ../chromosomes/1.fa > ./01_homer_all_fimo.txt`

In this example, we're using a .meme file of all HOMER motifs and chromosome 1, then dumping the output into an appropriately-labeled text file.

**Please note: this step will take a long time, so plan accordingly**

## 4. Copy the and restore the existing fimo database dump from Amazon S3

Ultimately, you'll want to add your text files created in step 3 to the fimo database. The most recent version of fimo was created on June 08, 2017 and can be copied to your machine as follows:

`aws s3 cp s3://cory-dbtest/2017_06_08_fimo.dump .`

The database is around 34 GB and will take some time to download. Once it does, you'll need to restore it in PostgreSQL using the following command:

`sudo pg_restore --verbose --clean --no-acl --no-owner --dbname=fimo --create 2017_06_08_fimo.dump`

**This command will also take quite a while to run (probably a couple of hours), so plan accordingly**

## 5. Using the template from [`create_fimo_table.sh`](https://github.com/PriceLab/BDDS/blob/master/trenadb/fimo/create_fimo_table.sh), create a shell script that copies the fimo output into the existing fimo database

## 6. Create indices using the commands in [index.sql](https://github.com/PriceLab/BDDS/blob/master/trenadb/fimo/index.sql)

## 7. Dump the new fimo database locally

## 8. Copy the database dump to Amazon S3

