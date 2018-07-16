#!/bin/bash
set -x

#If you have a single sample, run this script pointed towards the directory containing that sample's data. If you have more than one fastq file in that directory, use the full script in the parent directory.
#Usage: ./Example_Single_Sample.sh [pe/se] [full_path_to_reference_fasta] [path_to_directory_containing_reads] [path_to_directory_for_alignments] [YYYYMMDD date] [individual_identifier]

Pairing=$1
Reference=$2
Path_To_Input=$3
Alignment_Directory=$4
Date_Of_SNP_Calling=$5
Sample_ID=$6

mkdir ${Alignment_Directory}

if ls ${Reference}* | grep -q ${Reference}.bwt; then
	echo "There is already a bwa reference for this, continuing"
else
	echo "There is not already a bwa reference for this, creating one"
	programs/bwa index ${Reference}
fi

bwa mem -R "@RG\tID:${Sample_ID}\tSM:${Sample_ID}" ${Reference} ${Path_To_Input}/*R1* ${Path_To_Input}/*R2* > ${Alignment_Directory}/${Sample_ID}.sam
samtools view -F 4 -q 10 -bS ${Alignment_Directory}/${Sample_ID}.sam | samtools sort - ${Alignment_Directory}/${Sample_ID}
samtools index ${Alignment_Directory}/${Sample_ID}.bam
rm ${Alignment_Directory}/${Sample_ID}.sam
freebayes -f ${Reference} ${Alignment_Directory}/*bam > ${Date_Of_SNP_Calling}_${Sample_ID}.vcf
