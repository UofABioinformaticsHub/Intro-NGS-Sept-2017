
# Sequence Alignment

Once we have cleaned our data of any contaminating sequences, and removed the bases which are more likely to contain errors we can more confidently align our reads to a reference.  Different experiments may have different reference sequences depending on the context.  For example, if we have a sub-sample of the genome associated with restriction sites like RAD-Seq, we would probably align to a reference genome, or if we have RNA-Seq we might choose to align to the transcriptome instead of the whole genome.  Alternatively, we might be interested in de novo genome assembly where we have no reference genome to compare our data to.

## How Aligning Works

Most fast aligners in widespread public use are based on a technique called the Burrows-Wheeler Transform, which is essentially a way of restructuring, or indexing, the genome to allow very rapid searching.  This technique comes from computer science & is really beyond the scope of what most of us need to know.  The essence of it is that we have a very fast searching method, and most aligners use a seed sequence within each read to begin the searching.  These seeds are then expanded outwards to give the best mapping to a sequence.  There are many different alignment tools available today & each one will have a particular strength. For example, bowtie is very good for mapping short reads, whilst bowtie2 or bwa is more suited to reads longer than 50bp.

#### What’s the difference

Some key differences between aligners is in the way they index the genome, and in the a way they are equipped to handle mismatches & indels.  Choosing an aligner can be a difficult decision with the differences often being quite subtle.  Sometimes there is a best choice, other times there really isn’t.  Make sure you’ve researched relatively thoroughly before deciding which to use.

## Aligning our WGS reads

### Downloading A Reference Genome

To align any reads, we first need to download the appropriate (i.e.  latest) genome \& then we can build the index to enable fast searching via the Burrows-Wheeler Transform. Like we’ve seen in the previous sections, our reads today come from the nematode or Roundworm (*Caenorhabditis elegans*).  We have a copy of the genome read for you in the rawData directory,  however you can always redownload the genome (like you can do with all model genomes) by opening Firefox & head to [ftp://ftp.ensembl.org/pub/release-88/fasta/caenorhabditis_elegans/](ftp://ftp.ensembl.org/pub/release-88/fasta/caenorhabditis_elegans/).  (Here you can find the transcriptome (cdna), genome (dna), proteome (pep) & non-coding RNA (ncrna) in separate folders). With our C. elegans genome fasta file, we should have a quick look at the file.  It will contain all chromosomes & the mitochondrial sequences. We can print out the first 10 lines using the head command.

```
cd ~/rawData/
gunzip cel1.fa.gz
head cel1.fa
```

Note that the first line describes the following sequence & begins with a \> symbol.  We can use this to search within the file using regular expressions \& print all of these description lines.

```grep "^>" cel1.fa```

Alternatively could simply count them using the -c option.

```grep -c ">" cel1.fa```

## Building an Index

We will align using the tool bwa which is one of the original Burrows-Wheeler transformation mappers. bwa was developed in 2009 by Heng Li (Harvard/Broad Institute, USA). From the bwa manual, it details the three specific bwa algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate.” Today we will be using bwa-mem to align our C. elegans WGS reads. Once again, we need to check the help pages. Fortunately the bwa page is actually pretty friendly on the screen and appears without the usual -h option.

```
bwa
```

We should also inspect the help page for bwa index which we will use to build the index.

```
bwa index
```

Using this particular process you can usually just run the command on the fasta file and the index will be called by the same filename.  However in this case, we will name the index ”cel1” (thats cel for *C. elegans*) and ”1” for the version, by using the -p flag/parameter Now that we’ve had a look, type to following command which will take a few minutes to run.

```
bwa index -p Cel1 cel1.fa
```

Let’s look at what files have been created.

`ls`

You should be able to open a few of the files with the ”less” command, however the main files (\*.sa, \*.bwt and \*.pac) are the BWT transformed files that are in binary, so we can’t really see what they look like, but these are required by the aligner bwa.

## Aligning the reads

Because we only have a small subset of the actual sequencing run, we should be able to run this alignment in a reasonable period of time

```
bwa mem Cel1 WGS_SRR2003569_1.fastq.gz WGS_SRR2003569_2.fastq.gz | samtools view -bhS -> WGS_SRR2003569.bam
```

Let’s break down this command a little.  The first part of the command:

```
bwa mem Cel1 WGS/WGS_SRR2003569_1.fastq.gz WGS/WGS_SRR2003569_2.fastq.gz
```

will align our compressed sequenced reads to the Cel1 bwa index that we made. Usually you can create a SAM file (see next section) to store all the alignment data.  SAM files however are text files which can take up a significant amount of disk space, so its much more efficient to pipe it to the samtools command and create a compressed binary SAM file (called BAM). To do this, we run the program samtools:

```
samtools view -bhS - > WGS_SRR2003569.bam
```

In this context, samtools view is the general command that allows the conversion of the SAM to BAM. There is another more compressed version of the SAM file, called CRAM, which you can also create using samtools view.  However, we will not use that today.
