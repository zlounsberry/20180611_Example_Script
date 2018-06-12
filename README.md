# Example workflow for reference-guided variant calling 

See below for a list of programs that need installing (included here in the programs directory for convenience, but how to install for Ubuntu included below as well)
For programs in the `programs` directory, you may need to use `chmod +x [file]` each one to make them executable before running.
Once everything is installed, execute the line listed as "Example" below

## General usage: 
`./ExampleWorkflow.sh [paired=pe/single=se] [reference.fasta] [path to input file(s)] [path to output file(s)] [date you called SNPs as YEARMMDD] [name of project]`
## Example: 
`./ExampleWorkflow.sh se ./reference/canFam3_chr26.fa ./data ./alignments 20180611 Dog_SNP_Calling`

This example script takes several fastq files that represent R1's (Read 1) for 2 dogs from the DoGSD project (http://dogsd.big.ac.cn/dogsd/pages/download/fastq.js), aligns them to the canFam3 chr26 reference, and calls SNPs.

Dog subset files were created as follows (from full gzip downloads):

`for file in $(ls ~/Downloads/dog_00* | sed 's/^.*\///g' | sed 's/_1.fq.gz//g'); do zcat ~/Downloads/${file}_1.fq.gz | head -n800000 > ./${file}.R1.fastq; done`

`for file in $(ls dog* | sed 's/.R1.fastq//g'); do split -a1 -d -l200000 --additional-suffix=.fastq ${file}.R1.fastq data/${file}_R1_; done`

`rm *R1.fastq`

The following software is required to run this script 
NOTE: May require sudo privileges to install. Contact administrator for help if needed.
#### bwa (https://sourceforge.net/projects/bio-bwa/files/)
`apt-get install bwa`
#### samtools (https://sourceforge.net/projects/samtools/)
`apt-get install samtools`
#### freebayes (https://github.com/ekg/freebayes)
`git clone --recursive https://github.com/ekg/freebayes
cd freebayes
make`



### When you are finished, your resulting file should match `spoiler_output/20180611_Dog_SNP_Calling.vcf`
