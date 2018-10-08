# mitochondrial sequencing and analysis pipeline

A series of scripts and tools to perform mitochondrial sequence, alignment and analysis

The script `mtseq_processing.sh` details the process of creating appropriate index files, re-aligning to rCRS (if required), variant calling (generating vcf files) and creatation of fasta files. It also creates merged vcf and fasta files for batched experiments.

## required software

The pipeline requires a Linux environment (should probably run on MacOS, and 'maybe' on Linux-subsytem for Windows 10...). The following software tools are required:

  - samtools: http://www.htslib.org/
  - bcftools: http://www.htslib.org/
  - tabix: http://www.htslib.org/
  - bwa: https://github.com/lh3/bwa
  - GATK
    - latest version (try running version 3.8 first before installing 4+): https://software.broadinstitute.org/gatk/download/
    - version 3.8 (required to use `FastaAlternateReferenceMaker`): https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef
    
**Note:** if you are comfortable installing Anaconda (https://www.anaconda.com/) I recommend using Bioconda (https://bioconda.github.io/) to install and maintain the above tools, it gives a lot of flexibility and easeof-use (i.e. it's fairly simple to create and maintain multiple environments based off different Python versions).

## example of running a single sample

(note: further below is an edited version of a previous email which goes in depth on some issues and annotation)

If you want to check you have the correct tools installed and that the pipeline will run you can attempt processing a single sample. Here is some example code:

```sh
## example 'pipeline' for a single sample

# convert bam back to fastq
samtools bam2fq NI_110050-run15.1-ix20.bam > NI_110050_mt.fq

# align fastq to rCRS
bwa mem -t 4 rCRS.fa NI_110050_mt.fq | samtools sort -O BAM -o NI_110050_mt.bam 

# variant calling with left normalisation
samtools mpileup -E -u -v -p -f rCRS.fa NI_110050_mt.bam | \
  # normalise and call variants against ref genome
  bcftools norm -f rCRS.fa -m +both -d both | \
  bcftools call -v -m --ploidy 1 | \
  bcftools norm -f rCRS.fa -m -both | \
  bcftools norm -d both -O z -f rCRS.fa > NI_110050_mt.vcf.gz

# index vcf file
tabix NI_110050_mt.vcf.gz

# create fasta file from vcf file
java -jar ~/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R rCRS.fa -o NI_110050_mt.fa --variant NI_110050_mt.vcf.gz
```

**NOTE:** the above code assumes that you have the correct references, indexes and directories (GATK), and that all files are in the current working directory. This is different to the pipeline, which requires specific directories to be present. 

## 'Issues' with the Ion Torrent reference sequence

First, the ‘problem’ - on the torrent we use the RSRS consensus sequence to align all our mtDNA. Turns out a lot of tools only allow rCRS, making annotations a bit tricky when using tools like mitomaster.

Solution: remap the samples you’ve run (this is actually the easiest and fastest step of what is to follow).

### Process:

  - convert bam files back to fq (using samtools bam2fq)
  - remap the fq files using bwa mem (against rCRS this time)
  - run all new bam files through our existing pipeline to normalise and call variants (using rCRS)
    - individual fasta files are also generated
  - generated a merged vcf file for all samples

At this point we are back to where we want to be before going into annotation.

### QC and annotation

In terms of the QC and annotation I’ve added a few steps which I think provide some interesting outputs.

#### haplogrep 2 (https://haplogrep.uibk.ac.at/)
        
  - if you load the merged vcf file it will determine the haplotypes for each sample and also give you an idea of the accuracy and quality of the estimation.
    - you’ll see there are a few samples flagged in yellow and red - I suggest following these up (more on this later).
  -  you can output/export some nice information from here, including haplotype files and a phylogenetic tree with all your samples placed according to their mt haplotypes (very cool figure for your talks etc.).
  - generate/export a merged fasta file from here and you can import it into the next step.

#### mitomaster (https://mitomap.org/foswiki/bin/view/MITOMASTER/WebHome)
        
  - feeding the merged fasta file in should correctly annotate all variants.
  - you can export/download a cvs or excel file of the completed annotations which is really nice.

#### mitosuite (http://mitosuite.com/)

  - I’ve mentioned this tool before, it generates really nice QC results as well as haplotype calling and annotation.
  - I’ve spent too long trying to get this running on taurus:
    - the issue is that it only runs on Linux/Mac, but it requires a GUI which is problematic on a headless server like taurus… long story short maybe we could install it on the ‘new’ computer running Linux in Larisa’s office.
    - the other problem this software has is that because it has no command-line element it requires a user to manually run each sample through it. I’ve hacked together a method of making this easier and have run through your samples, but have reached ou to the authors to see if we can work on a better solution - it makes no sense in the age of high-throughput where you have 100s-1000s of samples that you have to do them one at a time!
    - in running it through your samples it’s obvious which ones have QC issues (we’ve discussed this before) so it’s worth going through the results files for each sample and noting which are a bit ‘dicey’. Some things I observed:
    - there are samples that only have sequence for ‘half’ the genome, obviously a set of primers didn’t amplify or something.
    - there is enrichment/excess sequence around primer binding sites, we see this routinely for torrent data - might be worth considering if we reduce the coverage to the average around these region - at least take care in calling/reporting on variants inside these regions.
    - there is a constant lack of sequence observed in every sample at position 3106 - I’m wondering if this is one of the ‘placeholder’ bases in the reference genome.
    - when looking at heteroplasmy there are over-representations in similar positions (309,310), I think these may be due to deletions/insertions and we might need to follow this up by testing a few different aligners - this shouldn’t cause you issues just looking at single point variation but something to keep in mind.
    - there are a non-trivial amount of samples that have magnitudes lower sequence amount than the others, be careful when trying to call variants in these.
