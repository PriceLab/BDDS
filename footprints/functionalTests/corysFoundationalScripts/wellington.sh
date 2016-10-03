# Line commands for creating the TFBS database
# I think I originally did this on CentOS, though I've since moved to Ubuntu. There may be some aspects of installing the programs that is different because of this

# F-seq required a special install package, ant, and my OS didn't have the right java version. This is how I fixed that
sudo -s

# Get the correct version of jdk (not sure if this is needed for Ubuntu)
aws s3 cp s3://cory-trena/jdk-7u79-linux-x64.gz .
tar xvfz jdk-7u79-linux-x64.gz

## modify bashrc to point to jdk
vi ~/.bashrc

export JAVA_HOME=/mnt/jdk1.7.0_79

## install apt

wget https://archive.apache.org/dist/ant/source/apache-ant-1.9.6-src.tar.bz2
tar xvfj apache-ant-1.9.6-src.tar.bz2
rm apache-ant-1.9.6-src.tar.bz2

wget http://anduin.linuxfromscratch.org/BLFS/other/junit-4.11.jar
wget http://hamcrest.googlecode.com/files/hamcrest-1.3.tgz

tar -xvf hamcrest-1.3.tgz
rm hamcrest-1.3.tgz
cd ..
mkdir /lib/optional
cd /mnt
cp -v junit-4.11.jar hamcrest-1.3/hamcrest-core-1.3.jar /lib/optional

cd apache-ant-1.9.6/
./build.sh -Ddist.dir=/opt/ant-1.9.6 dist &&
ln -v -sfn ant-1.9.6 /opt/ant

## install fseq
cd /home
git clone https://github.com/aboyle/F-seq.git
cd F-seq/
# ant is the weird install package needed for F-seq
ant
cd dist~/
tar -xvf fseq.tgz
echo "export PATH=\$PATH:/mnt/F-seq/dist~/fseq/bin/" >> ~/.bashrc

# increase memory for fseq
# should be a sed command, I did it manually. 
cd /home/F-seq/dist~/fseq/bin
JAVAOPTS="-Xmx32g"

### install samtools

cd /home
sudo wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2
sudo tar xvfj samtools-1.3.tar.bz2
cd samtools-1.3
sudo make
echo "export PATH=\$PATH:/home/samtools-1.3/" >> ~/.bashrc
source ~/.bashrc

## install htslib (includes tabix
wget https://github.com/samtools/htslib/releases/download/1.3/htslib-1.3.tar.bz2
tar xvfj htslib-1.3.tar.bz2
cd htslib-1.3
make
make prefix=/user/ install 

## if bedtools not installed...
git clone https://github.com/arq5x/bedtools2.git
cd bedtools2
make
make prefix=/usr/ install

## install pyBedTools and pyDNase
pip install pybedtools
pip install pyDNase

## install SNAP
sudo git clone https://github.com/amplab/snap.git
cd snap
sudo make
echo "export PATH=\$PATH:/home/snap/" >> ~/.bashrc
source ~/.bashrc

## install SRA-toolkit
### install sra-toolkit
cd /home
sudo wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
sudo tar -vxzf sratoolkit.tar.gz

echo "export PATH=\$PATH:/home/sratoolkit.2.5.7-ubuntu64/bin/" >> ~/.bashrc
source ~/.bashrc


############################

# Create indices for SNAP

mkdir /scratch/resources
aws s3 cp s3://snapr-ref-assemblies/human/Homo_sapiens.GRCh38.dna.SORTED.fa /scratch/resources/
aws s3 cp s3://snapr-ref-assemblies/human/Homo_sapiens.GRCh38.83.gtf.gz /scratch/resources/
gunzip /scratch/resources/Homo_sapiens.GRCh38.83.gtf.gz

mkdir /scratch/genome
snap-aligner index /scratch/resources/Homo_sapiens.GRCh38.dna.SORTED.fa /scratch/genome/ -s 16 -locationSize 5 -bSpace

# run SNAP aligner (this is for single reads)
snap-aligner single /scratch/genome/ /scratch/sra_data/.fastq.gz -o /scratch/snap_out/SRR412463.bam -mrl 26

#######

# Once I had the aligned bam files, I converted them to bed files and merged them prior to running F-seq and Wellington

#!/bin/bash

ls -lt /mnt3/tissue/*bam | rev | cut -f 1 -d "/" | cut -f 2 -d "." | rev > bam_list

for file in `cat bam_list`
do
                NAME=`echo $file | cut -f 1 -d "."`
                echo "$NAME"
                samtools index /mnt3/tissue/$NAME.bam
                bedtools bamtobed -i /mnt3/tissue/$NAME.bam > /mnt/tissue/$NAME.bed # notice the switch from mnt3 to mnt
                mkdir /mnt/tissue
                mkdir /mnt/tissue/$NAME-fseq
                cd /mnt
                fseq -f 0 -o /mnt/tissue/$NAME-fseq -t 2.5 -of bed /mnt/tissue/$NAME.bed
                cd /mnt/tissue/$NAME-fseq
				cat *.bed > $NAME.peaks.bed
				cut -f 1,2,3 $NAME.peaks.bed > $NAME.peaks.std.bed
                awk '{a=$2-$3;print $0 "\t" a;}' $NAME.peaks.std.bed > $NAME.peaks.dif.bed
                cat $NAME.peaks.dif.bed | awk '($4 < -400)' | cut -f 1,2,3 > $NAME.peaks.400.bed
                cd ../
                mkdir /mnt/tissue/$NAME-well_out
                wellington_footprints.py /mnt/tissue/$NAME-fseq/$NAME.peaks.400.bed $file /mnt/tissue/$NAME-well_out
done


###################

# FIMO search
# Not sure if we want to have this as part of the pipeline as we already have run this, but we are thinking of parameter tuning this and changing how it's done.
# I broke up the genome into groups of chromosomes with the groups being close in size. Here I ran it on hg19, but have since moved to GRCh38. I only ran it on the primary assembly (no additional haplotypes)
# I need to get you the JASPAR file with that also has some of the MEME TFs that is used for this.

#!/bin/bash

# attempt to run FIMO on each chromosome

# create directory for the work
cd /mnt
mkdir fimo
cd fimo

# download genome file from s3
aws s3 cp s3://snapr-ref-assemblies/human/Homo_sapiens.GRCh38.dna.SORTED.fa .

# split the genome file into chromosomes. each chromosome will be named like this: 10.fa
cat Homo_sapiens.GRCh38.dna.SORTED.fa | awk 'BEGIN { CHROM="" } { if ($1~"^>") CHROM=substr($1,2); print $0 > CHROM".fa" }'

# get the size and name of each chromosome
ls -lt | grep fa | grep -v Homo | rev | cut -f 1,6 -d " " | rev | sort -n > c_list

# combine chromosomes so to make files of comparable size
cat 13.fa 14.fa > 13-14.fa
cat 15.fa 16.fa > 15-16.fa
cat 17.fa 18.fa 20.fa > 17_18_20.fa
cat 19.fa Y.fa 21.fa 22.fa MT.fa > 19_21_22_Y_MT.fa
rm 13.fa
rm 14.fa
rm 15.fa
rm 16.fa
rm 17.fa
rm 18.fa
rm 19.fa
rm 20.fa
rm 21.fa
rm 22.fa
rm Y.fa
rm MT.fa

# for the non-redundant additional motifs
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/fimo_out_nonred/ --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/Homo_sapiens.GRCh37.68.genome.fa & 

nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_X --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/X.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_1 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/1.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_2 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/2.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_3 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/3.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_4 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/4.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_5 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/5.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_6 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/6.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_7 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/7.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_8 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/8.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_9 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/9.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_10 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/10.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_11 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/11.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_12 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/12.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_13-14 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/13-14.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_15-16 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/15-16.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_17-18-20 --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/17_18_20.fa &
nohup fimo --oc /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_19-21-22-Y-MT --thresh 1e-5 /local/Cory/fimo/nonredundant_motifs.meme /local/Cory/fimo/GRCh37/19_21_22_Y_MT.fa &

# concatinate the rest of the fimo outputs
cat /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_1/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_2/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_3/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_4/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_5/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_6/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_7/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_8/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_9/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_10/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_11/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_12/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_13-14/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_15-16/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_17-18-20/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_19-21-22-Y-MT/fimo.gff /local/Cory/fimo/alt_pvalue/new_jaspar/fimo_out_X/fimo.gff > new_nonred_main.gff

##############################

### getting ready to do a bedtools intersect of the fimo output and wellington fp

# concatinate the nonredundant motifs Seth put together
cat /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_11/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_12/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_10/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_9/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_X/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_8/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_7/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_5/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_6/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_4/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_15-16/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_19-21-22-Y-MT/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_3/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_17-18-20/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_2/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_1/fimo.gff /local/Cory/fimo/GRCh37/fimo_dir/nonred/fimo_out_13-14/fimo.gff > /local/Cory/fimo/GRCh37/fimo_dir/nonred/combined/nonred.gff

# concatinate the rest of the fimo outputs
cat /local/Cory/fimo/alt_pvalue/fimo_out_1/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_2/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_3/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_4/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_5/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_6/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_7/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_8/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_9/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_10/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_11/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_12/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_13-14/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_15-16/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_17-18-20/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_19-21-22-Y-MT/fimo.gff /local/Cory/fimo/alt_pvalue/fimo_out_X/fimo.gff > /local/Cory/fimo/GRCh37/fimo_dir/nonred/combined/nonred_main.gff

# combine the two files
cat nonred_main.gff nonred.gff > fimo_all_temp.gff

# get rid of a couple of extraneous columns
cat fimo_all_temp.gff | cut -f 1,4,5,6,7,8,9 > fimo_all_temp2.gff 

# needed for sort to properly work
export LC_ALL=C

# sort the concatinated bed file
sort -k1,1 -k2,2n fimo_all_temp2.gff > fimo_all.gff

#### script to modify the wellington files

#!/bin/bash

# download the footprint files from s3 (getting the most liberal call)
aws s3 ls s3://cory-tfdb/GRCh38-results/ --recursive | grep WellingtonFootprints.-10.bed | rev | cut -f 1-3 -d " " | rev > file_list
FILES=`cat file_list`
# because the directory name includes p value cutoff (including spaces) I need to specify to loop over every /n
IFS=$'\n'

for i in $FILES
do
        aws s3 cp s3://cory-tfdb/$i .
done

files=`ls *.bed`
IFS=$'\n'

for i in $files
do
                prefix=`echo $i | cut -f 1 -d "."`
                echo "$prefix"i
                # because the bed file uses chr at the begining, I need to get rid of it so it will match the downstream files that don't include chr
                cat $i | cut -f 2- -d "r" > $prefix.400.WellingtonFP.temp.-10.bed
                export LC_ALL=C
                sort -k1,1 -k2,2n $prefix.400.WellingtonFP.temp.-10.bed > $prefix.400.WellingtonFP.-10.bed
                rm $i
                rm $prefix.400.WellingtonFP.temp.-10.bed
done

#####

# do intersect with bed tools
bedtools intersect -wo -sorted -a fimo_sort.gff -b well_out/*.bed > fimo_fp.bed

# replace one of the columns with just the TF in preparation to sorting on uniqueness. The sed part fixes the spacing (converts space to tab) which awk messes up.
cat fimo_fp.bed | cut -f 7 | cut -f 1 -d ";" | cut -f 2 -d "=" > tf_list
awk 'FNR==NR{a[NR]=$1;next}{$6=a[FNR]}1' tf_list fimo_fp.bed | sed -e 's/ [ ]*/\t/g' > fimo_fp2.bed


# for the purpose of comparing to the old tfdb, get a unique list of footprints on motifs
sort -u -k1,1 -k2,2 -k6,6 fimo_fp2.bed > uniq_fimo_fp.bed


#### pulling out the TSS for each gene from the gtf file

# starting with Homo_sapiens.GRCh37.75.gtf because it's the last update for GRCh37
# now doing it with Homo_sapiens.GRCh38.83.gtf
# pull out just the "gene" and "protein_coding" fields
#cat Homo_sapiens.GRCh38.83.gtf | awk '$3 ~ /gene/' | awk /protein_coding/ > GRCh38.83.gene.gtf
cat Homo_sapiens.GRCh38.83.gtf | awk '$3 ~ /transcript/' > GRCh38.83.transcript.gtf

# create a list of unique gene_ids (not transcripts) by which to do the TSS grouping
cat GRCh38.83.gene.gtf | cut -f 2 -d "\"" | sort -u > gene_list



# separate + and - strands
cat GRCh38.83.transcript.gtf | awk '$7 ~ /\+/' > GRCh38.83.transcript.plus.gtf
cat GRCh38.83.transcript.gtf | awk '$7 ~ /\-/' > GRCh38.83.transcript.minus.gtf

# create a bed file that has +/- 1mb TSS for each gene
# sed removes unwanted " and ; transcript is field 14 and gene is field 10
cat GRCh38.83.transcript.plus.gtf | awk '{ print $1,"\t", $4 - 1000000,"\t", $4 + 1000000,"\t",$14,"\t",$12 }' | sed 's|[";]||g' > GRCh38.83.plus.1mb.bed
cat GRCh38.83.transcript.minus.gtf | awk '{ print $1,"\t", $5 - 1000000,"\t", $5 + 1000000,"\t",$14,"\t",$12 }' | sed 's|[";]||g' > GRCh38.83.minus.1mb.bed

#### for just getting the TSS of every transcript:
cat GRCh38.83.transcript.plus.gtf | awk '{ print $1,"\t", $4,"\t",$14,"\t",$10 }' | sed 's|[";]||g' > GRCh38.83.plus.1mb.bed
cat GRCh38.83.transcript.minus.gtf | awk '{ print $1,"\t", $5,"\t",$14,"\t",$10 }' | sed 's|[";]||g' > GRCh38.83.minus.1mb.bed
#combine and sort
cat GRCh38.83.minus.1mb.bed GRCh38.83.plus.1mb.bed > GRCh38.83.full.1mb.bed
sort -k1,1 -k2,2n GRCh38.83.full.1mb.bed > GRCh38.83.full.sort.1mb.bed
cat GRCh38.83.full.sort.1mb.bed | awk '{printf ("%s\t%s\t%s\t%s\n", $1,$2,$3,$4,$5)}' > GRCh38.83.tss.1mb.bed
####

#combine and sort
cat GRCh38.83.minus.1mb.bed GRCh38.83.plus.1mb.bed > GRCh38.83.full.1mb.bed
sort -k1,1 -k2,2n GRCh38.83.full.1mb.bed > GRCh38.83.full.sort.1mb.bed
# get the format right for bed
cat GRCh38.83.full.sort.1mb.bed | awk '{printf ("%s\t%s\t%s\t%s\n", $1,$2,$3,$4,$5)}' > GRCh38.83.full.sort2.1mb.bed

# get rid of negative numbers (which introduces spaces which need to be replaced with tabs)
cat GRCh38.83.full.sort2.1mb.bed | awk '$2<0 {$2=0} 1' | sed -e 's/ [ ]*/\t/g' > GRCh38.83.1mb.bed

# money shot
bedtools intersect -wo -a GRCh38.83.1mb.bed -b uniq_fimo_fp.bed > tfdb_1mb_GRCh38.83.bed

#calc number of unique genes
sort -u -k1,1 -k2,2n tfdb_1mb_GRCh38.83bed > uniq_tfdb_1mb_GRCh38.83.-20.bed

