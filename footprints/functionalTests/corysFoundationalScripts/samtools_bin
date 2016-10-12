#/bin/bash/

# define bin
bin=1000
# length of chr19 is 58617616
bin_length=$((58617616 / $bin))
#echo $bin_length

START=0
END=$bin
 
for (( c=$START; c<=$END; c++ ))
do
	temp_start=$((($bin_length * $c) + 1 ))
	temp_end=$(($bin_length * ( 1 + $c)))
	samtools view -h ENCSR000DBY.19.bam "chr19:$temp_start-$temp_end" > chr19_"$temp_start"_"$temp_end".bam
	samtools view -F 0x4 chr19_"$temp_start"_"$temp_end".bam | wc -l
	rm chr19_"$temp_start"_"$temp_end".bam
done
