# mitochondrial sequencing and analysis pipeline

A series of scripts and tools to perform mitochondrial sequence, alignment and analysis

The script `mtseq_processing.sh` details the process of creating appropriate index files, re-aligning to rCRS (if required), variant calling (generating vcf files) and creatation of fasta files. It also creates merged vcf and fasta files for batched experiments.

(note: the below is an edited version of a previous email)

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

#### haplogrep 2
        
  - if you load the merged vcf file it will determine the haplotypes for each sample and also give you an idea of the accuracy and quality of the estimation.
    - you’ll see there are a few samples flagged in yellow and red - I suggest following these up (more on this later).
  -  you can output/export some nice information from here, including haplotype files and a phylogenetic tree with all your samples placed according to their mt haplotypes (very cool figure for your talks etc.).
  - generate/export a merged fasta file from here and you can import it into the next step.

#### mitomaster
        
  - feeding the merged fasta file in should correctly annotate all variants.
  - you can export/download a cvs or excel file of the completed annotations which is really nice.

#### mitosuite

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
