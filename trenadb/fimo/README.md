## Instructions for creating/updating fimo

When new motifs become available, they should be added to the existing fimo database. For the purposes of this short tutorial, we'll assume said motifs have already been added to the PriceLab/MotifDb repository. If not, please refer to that repository, located [here](https://github.com/PriceLab/MotifDb). At this point, the workflow for doing so is as follows:


1. Clone the latest MotifDb version from the forked PriceLab repo
2. Using MotifDb in R, export .meme files of your desired motifs
3. Using the fimo program, intersect the motif .meme files with all chromosomes
4. Copy the existing fimo database dump from Amazon S3
5. Restore the database (I do this on an Amazon EC2 instance)
6. Using the template from [`create_fimo_table.sh`](https://github.com/PriceLab/BDDS/blob/master/trenadb/fimo/create_fimo_table.sh), create a shell script that copies the fimo output into the existing fimo database
7. Create indices using the commands in [index.sql](https://github.com/PriceLab/BDDS/blob/master/trenadb/fimo/index.sql)
8. Dump the new fimo database locally
9. Copy the database dump to Amazon S3

What follows is a more detailed description of how to carry out all these steps, with examples using the fimo version created on June 8, 2017. 

## 1

## 2

## 3

## 4

## 5

## 6

## 7

## 8
