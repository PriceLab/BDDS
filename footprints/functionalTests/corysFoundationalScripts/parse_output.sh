#/bin/bash/

MOTIF_NUMBER=$1
BAM_FILE_NAME=$2
OUTPUT_DIR=/scratch/piq/results-$MOTIF_NUMBER-$BAM_FILE_NAME
#BAM_DIR=/scratch/forpaul/
BAM_DIR=/scratch/lympho_bam/mod_bam
mkdir -p $OUTPUT_DIR
mkdir -p /scratch/piq/agg
mkdir -p /scratch/piq/tmp
cd /home/piq-single
Rscript pwmmatch.exact.r common.r /home/parsing/piq/jaspar_2016_forPIQ.txt $MOTIF_NUMBER $OUTPUT_DIR
wait
Rscript bam2rdata.r common.r $OUTPUT_DIR/$BAM_FILE_NAME.RData $BAM_DIR/$BAM_FILE_NAME.bam
wait
Rscript pertf.r common.r $OUTPUT_DIR/ /scratch/piq/tmp $OUTPUT_DIR $OUTPUT_DIR/$BAM_FILE_NAME.RData $MOTIF_NUMBER
wait

FILE_NAME=`ls -1 $OUTPUT_DIR/*.RC-calls.csv | rev | cut -f 1 -d "/" | rev | cut -f 2 -d "-" | rev | cut -c 4- | rev`
echo "$FILE_NAME"
sed 's/,/\t/g' $OUTPUT_DIR/$MOTIF_NUMBER-$FILE_NAME-calls.all.csv | sed 's/""/num/g' | sed 's/"//g' > $OUTPUT_DIR/csv1
sed -i '1d' $OUTPUT_DIR/csv1
sed 's/,/\t/g' $OUTPUT_DIR/$MOTIF_NUMBER-$FILE_NAME.RC-calls.all.csv | sed 's/""/num/g' | sed 's/"//g' > $OUTPUT_DIR/csv2
sed -i '1d' $OUTPUT_DIR/csv2
cat $OUTPUT_DIR/csv1 $OUTPUT_DIR/csv2 > $OUTPUT_DIR/csv3
sed '1d' $OUTPUT_DIR/$MOTIF_NUMBER-$FILE_NAME-calls.all.bed > $OUTPUT_DIR/temp1
sed '1d' $OUTPUT_DIR/$MOTIF_NUMBER-$FILE_NAME.RC-calls.all.bed > $OUTPUT_DIR/temp2
cat $OUTPUT_DIR/temp1 $OUTPUT_DIR/temp2 > $OUTPUT_DIR/bed1
paste $OUTPUT_DIR/bed1 $OUTPUT_DIR/csv3 > $OUTPUT_DIR/all_temp
awk -F "\t" '{print $1,$2,$3,$4,$6,$10,$11,$12,$13}' $OUTPUT_DIR/all_temp | tr ' ' '\t' | sed 's/.RC//g' | sort -k1,1V -k2,2n > /scratch/piq/agg/$FILE_NAME.$BAM_FILE_NAME.bed
aws s3 cp /scratch/piq/agg/ s3://cory-temp/piq_out/ --exclude "*" --include "$FILE_NAME.$BAM_FILE_NAME.bed" --recursive

