#!/bin/bash
set -x

Pairing=$1
Reference=$2
Path_To_Input=$3
Alignment_Directory=$4
Date_Of_SNP_Calling=$5
Project_ID=$6

mkdir ${Alignment_Directory}

##Check if a BWA reference is already present. If not, make one.
if ls ${Reference}* | grep -q ${Reference}.bwt; then
	echo "There is already a bwa reference for this, continuing"
else
	echo "There is not already a bwa reference for this, creating one"
	programs/bwa index ${Reference}
fi

##Get a list of your files. This is the bit that you may need to change depending on the naming structure of your data
##As long as your files are in pretty standard (e.g., contain R1 (and R2, if paired) in file names) and do not contain ".f" outside the file extension this will work fine. gzipped fastq are fine too.
if [[ ${Pairing} == "pe" ]]; then
	#Create a 3-column file "files.txt" containing your R1 files, R2 files, and Readgroup ID's that will serve as the input for the main loop
	paste <(ls ${Path_To_Input}/*R1*) <(ls ${Path_To_Input}/*R2*) <(ls ${Path_To_Input}/*R1* | sed 's/^.*[\/]//' | sed 's/\.f.*.$//') > files.txt
else
	 if [[ ${Pairing} == "se" ]]; then
		#Create a 2-column file "files.txt" containing your R1 files and Readgroup ID's that will serve as the input for the main loop
		paste <(ls ${Path_To_Input}/*R1*) <(ls ${Path_To_Input}/*R1* | sed 's/^.*[\/]//' | sed 's/\.f.*.$//') > files.txt
	else
		#Echo polite error and exit
		echo "Sorry for the confusion, but you need to specify if your input fastq files are paired end 'pe' or single end 'se' (see: usage)"
		exit 0
	fi
fi

##Count the number of files you have, this loop will run on the file you created (files.txt)
##If your data are paired, it will run this loop:
if [[ ${Pairing} == "pe" ]]; then
	Lines=$(wc -l < files.txt) #Count the number of fastq files you specified in input directory and set it as a variables '${Lines}'
	x=1 #set a variable x to 1
	while [ $x -le $Lines ] #initiate count from 1 to number of fastq files
	do
		string="sed -n ${x}p files.txt" #create a string with the x'th line of your file
	        str=$($string)
	        var=$(echo $str | awk -F"\t" '{print $1, $2, $3}') #split the line such that each column is a variable
	        set -- $var #set a series of variables
	        R1=$1 #the first variable (column 1 in files.txt)
		R2=$2 #second variable (column 2 in files.txt) This script reads one full line at a time and executes the loop on each column.
		ID=$3 #you guessed it! The third column in files.txt

		#bwa with default parameters, but adding a readgroup to each sample so the resulting vcf will have a column for each sample
		programs/bwa mem -R "@RG\tID:${ID}\tSM:${ID}" ${Reference} ${R1} ${R2} > ${Alignment_Directory}/${ID}.sam

		#samtools: -F 4 flag tells it to remove unmapped reads; -q 10 tells it to remove reads that have a 0.1 likelihood to be mapped incorrectly (based on MAPQ score)
		programs/samtools view -F 4 -q 10 -bS ${Alignment_Directory}/${ID}.sam | programs/samtools sort - ${Alignment_Directory}/${ID}

		##index BAM file
		programs/samtools index ${Alignment_Directory}/${ID}.bam

		#clean up temp file
		rm ${Alignment_Directory}/${ID}.sam

	x=$(( $x + 1 )) #loop this until x = the number of fastq files and terminate the loop at the end
	done

##if data are specified as single-end (see section above for comments if you are curious about specific lines):
else
	Lines=$(wc -l < files.txt)
        x=1
        while [ $x -le $Lines ]
        do
		string="sed -n ${x}p files.txt"
		str=$($string)
		var=$(echo $str | awk -F"\t" '{print $1, $2}')
		set -- $var
		R1=$1
		ID=$2

	        programs/bwa mem -R "@RG\tID:${ID}\tSM:${ID}" ${Reference} ${R1} > ${Alignment_Directory}/${ID}.sam
	        programs/samtools view -F 4 -q 10 -bS ${Alignment_Directory}/${ID}.sam | programs/samtools sort - ${Alignment_Directory}/${ID}
	        programs/samtools index ${Alignment_Directory}/${ID}.bam
	        rm ${Alignment_Directory}/${ID}.sam

	x=$(( $x + 1 ))
	done
fi

echo -e "\n\nIf the next line gives you a 'generating faidx' warning, that's okay! It is using samtools faidx to index your reference for you.\n\n"
##Run freebayes using default options
##Alternatively, use some options! They are too project-specific for me to generically include, but check out https://github.com/ekg/freebayes and also see what papers doing similar things are doing.
##It should throw an erro
programs/freebayes -f ${Reference} ${Alignment_Directory}/*bam > ${Date_Of_SNP_Calling}_${Project_ID}.vcf

echo -e "\n\nAll finished! Consider parsing the resulting file using vcflib (https://github.com/vcflib/vcflib) or vcftools (http://vcftools.sourceforge.net/)\n"
