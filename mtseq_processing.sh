####
## example pipeline for mtDNA sequencing analysis
## author: Miles Benton
## created: 170901
## modified: 181008
## 
## WARNING: this is a work in progress script, as such there may still be 'hard coded' links to software/data 
## sources that will need modifying on a per use basis. 

# NOTE: ensure you have a directory called 'bam', this contains bam files for all samples
# also create a directory called 'vcf', this will house resultant vcf files
# and a directory called 'fasta', which will contain fa/fasta files

## the below will check for these directories and make them if they don't exist
mkdir -p bam fasta vcf

# NOTE: for this script to run there must be at least one bam file in the 'bam' directory

## only need to do the reference formatting/indexing once 
# make sure reference sequence is formattd correctly and indexed
bgzip -c rCRS.fa > rCRS.fa.gz
samtools faidx rCRS.fa.gz
# create a GATK dictionary for rCRS
/usr/lib/jvm/java-8-openjdk-amd64/bin/java -jar ~/Downloads/software/GATK/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar CreateSequenceDictionary -R rCRS.fa -O rCRS.dict
##

## create a file containing a list of all samples to be processed
ls bam/*.bam | sed -e 's/.bam//g' | sed -e 's/bam\///g' > bam/bam_list.txt

## NOTE: run the below if you're sequence isn't mapped to the rCRS
# convert bam files back to fq
while read x; do
  echo -e "... running samtools bam2fq on ${x} ..."
  samtools bam2fq bam/"${x}".bam > fastq/"${x}".fq
done <bam/bam_list.txt
# align using bwa mem to rCRS and output as sorted bam
while read x; do
  bwa mem -t 4 rCRS.fa fastq/"${x}".fq | samtools sort -O BAM -o bam/"${x}".bam -
done <bam_list.txt
# index bam files
while read x; do
  samtools index bam/"${x}".bam
done <bam_list.txt
##

## after you have rCRS aligned data move onto the next step

## normalise calls from mpileup and index
while read x; do
  # samtools mpileup
  echo -e "... running samtools mpileup, variant calling and normalisation on ${x} ..."
  samtools mpileup -E -u -v -p -f rCRS.fa bam/"${x}".bam | \
  # normalise and call variants against ref genome
  bcftools norm -f rCRS.fa -m +both -d both | \
  bcftools call -v -m --ploidy 1 | \
  bcftools norm -f rCRS.fa -m -both | \
  bcftools norm -d both -O z -f rCRS.fa > vcf/"${x}".vcf.gz
  # create an index for the above
  echo -e "... creating index for ${x}.vcf.gz ..."
  tabix vcf/"${x}".vcf.gz
done <bam/bam_list.txt

## create fasta file of each sample
while read x; do
  echo -e "... creating fasta file for ${x} ..."
  /usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/java -jar ~/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R rCRS.fa -o fa/"${x}".fa --variant vcf/"${x}".vcf.gz
done <bam/bam_list.txt

## create a merged vcf
vcf-merge --ref-for-missing 0 vcf/*.vcf.gz > all_samples_merged.vcf

## create a merged fasta file
cat rCRS.fa fa/*.fa > fa/merged_fasta.fa
##
##/END
