#!/bin/bash

# Veronika Bernhauerova   October 28, 2018
# Last update:     October 28, 2018


usage=$(cat <<'END_USAGE'

 Veronika Bernhauerova October 28, 2018

  The script will run only if the reference genome file and sequence file are in the same path
  or directory.

  Usage: call_consensus.sh [-h] <reference> <bwt2index> <stock>
 
   -h         print this help 

   <reference>  indexed reference file only in FASTA format; input name ID without .fasta extension

   <bwt2index> bowtie2 reference index

   <stock>      sequence file only in FASTQ.GZ format; input name ID without .FASTQ.GZ extension
  
   Examples:
   /path/to/call_consensus.sh old_reference_name sequences_name 

END_USAGE
)

if [ $# -eq 0 ]
then
   echo -e "$usage\n"
   exit
fi

while getopts "hswa:e:" optionName; do
case "$optionName" in
h) echo "$usage\n"; exit;;
[?]) echo "error: invalid argument"; echo "$usage\n"; exit;;
esac
done

### input file name
old_reference=$1
index=$2
stock=$3

inRef=${old_reference}.fasta
indexRef=${index}
inFile=${stock}.fastq.gz

### BOWTIE2 alignment 
bowtie2 -x $indexRef -U $inFile -p 4 | samtools sort -o ${stock}.bowtie2.bam -O bam -T bowtie2.deleteme

### index sorted bam 
samtools index ${stock}.bowtie2.bam;

### generate a variant calling file from the sorted bam
bcftools mpileup -f $inRef ${stock}.bowtie2.bam | bcftools call -mv -Oz -o ${stock}.vcf.gz

### index the variant calling file
bcftools index ${stock}.vcf.gz

### normalize indels
bcftools norm -f $inRef ${stock}.vcf.gz -Ob -o ${stock}.norm.bcf

### filter adjacent indels within 5bp
bcftools filter --IndelGap 5 ${stock}.norm.bcf -Ob -o ${stock}.norm.flt-indels.bcf

### get the genome consensus
cat $inRef | bcftools consensus ${stock}.vcf.gz > ${old_reference}_consensus.fasta


### align the old reference to new reference and index the bam
bowtie2 -x $indexRef -f ${old_reference}_consensus.fasta -p 4 | samtools sort -o ${old_reference}_consensus.bowtie2.bam -O bam -T bowtie2.deleteme
samtools index ${old_reference}_consensus.bowtie2.bam

### remove all useless files
rm ${stock}.bowtie2.bam
rm ${stock}.bowtie2.bam.bai
rm ${stock}.vcf.gz
rm ${stock}.vcf.gz.csi
rm ${stock}.norm.bcf
rm ${stock}.norm.flt-indels.bcf

exit 0