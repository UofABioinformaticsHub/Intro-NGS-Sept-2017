* TOC
{:toc}

# Sequence Alignment

Once we have cleaned our data of any contaminating sequences, and removed the bases which are more likely to contain errors we can more confidently align our reads to a reference.  Different experiments may have different reference sequences depending on the context.  For example, if we have a sub-sample of the genome associated with restriction sites like RAD-Seq, we would probably align to a reference genome, or if we have RNA-Seq we might choose to align to the transcriptome instead of the whole genome.  Alternatively, we might be interested in de novo genome assembly where we have no reference genome to compare our data to.

## How Aligning Works

Most fast aligners in widespread public use are based on a technique called the Burrows-Wheeler Transform, which is essentially a way of restructuring, or indexing, the genome to allow very rapid searching.  This technique comes from computer science & is really beyond the scope of what most of us need to know.  The essence of it is that we have a very fast searching method, and most aligners use a seed sequence within each read to begin the searching.  These seeds are then expanded outwards to give the best mapping to a sequence.  There are many different alignment tools available today & each one will have a particular strength. For example, bowtie is very good for mapping short reads, whilst bowtie2 or bwa is more suited to reads longer than 50bp.

#### What’s the difference

Some key differences between aligners is in the way they index the genome, and in the a way they are equipped to handle mismatches & InDels.  Choosing an aligner can be a difficult decision with the differences often being quite subtle.  Sometimes there is a best choice, other times there really isn’t.  Make sure you’ve researched relatively thoroughly before deciding which to use.

## Aligning our WGS reads

### Downloading A Reference Genome

To align any reads, we first need to download the appropriate (i.e.  latest) genome \& then we can build the index to enable fast searching via the Burrows-Wheeler Transform. Like we’ve seen in the previous sections, our reads today come from the nematode or Roundworm (*Caenorhabditis elegans*).  

**Note**: If you want the full genome sequence you can use the command-line program `wget` to download the *C. elegans* genome sequence. If `wget` doesnt work for you, you can always you can always redownload the genome (like you can do with all model genomes) by opening Firefox & head to [ftp://ftp.ensembl.org/pub/release-90/fasta/caenorhabditis_elegans/](ftp://ftp.ensembl.org/pub/release-90/fasta/caenorhabditis_elegans/).  

For todays tutorial, we've given you just the sequence of chrI.
It may have been accidentally saved as the file `WGS` so if you have a file called `WGS` and can't see this file call an instructor over.

```
# Have a look at the first few lines
cd ~/WGS
head chrI.fa
```

Note that the first line describes the following sequence & begins with a \> symbol.  We can use this to search within the file using regular expressions \& print all of these description lines.


## Building an Index

We will align using the tool bwa which is one of the original Burrows-Wheeler transformation mappers. bwa was developed in 2009 by Heng Li (Harvard/Broad Institute, USA). From the bwa manual, it details the three specific bwa algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate.” Today we will be using bwa-mem to align our C. elegans WGS reads. Once again, we need to check the help pages. Fortunately the bwa page is actually pretty friendly on the screen and appears without the usual -h option.

```
bwa
```

We should also inspect the help page for bwa index which we will use to build the index.

```
bwa index
```

Using this particular process you can usually just run the command on the fasta file and the index will be called by the same filename.  However in this case, we will name the index "Celegans_chrI" by using the `-p` flag/parameter Now that we’ve had a look, type to following command which will take a few minutes to run.

```
bwa index ~/WGS/chrI.fa -p Celegans_chrI
```

Let’s look at what files have been created.

```
ls
```

You should be able to open a few of the files with the ”less” command, however the main files (\*.sa, \*.bwt and \*.pac) are the BWT transformed files that are in binary, so we can’t really see what they look like, but these are required by the aligner bwa.

## Aligning the reads

Because we only have a small subset of the actual sequencing run, we should be able to run this alignment in a reasonable period of time

```
cd ~/WGS/trimmedData/fastq
bwa mem -t 4 ~/WGS/Celegans_chrI SRR2003569_sub_1.fastq.gz SRR2003569_sub_2.fastq.gz | samtools view -bhS -F4 -> SRR2003569_chI.bam
```

Let’s break down this command a little.  The first part of the command:

```
bwa mem -t 4 ~/WGS/Celegans_chrI SRR2003569_sub_1.fastq.gz SRR2003569_sub_2.fastq.gz
```

will align our compressed sequenced reads to the Celegans_chrI `bwa` index that we made. Usually you can create a SAM file (see next section) to store all the alignment data.  SAM files however are text files which can take up a significant amount of disk space, so its much more efficient to pipe it to the `samtools` command and create a compressed binary SAM file (called BAM). To do this, we run the program `samtools`:

```
samtools view -bhS - > SRR2003569_chI.bam
```

In this context, `samtools` view is the general command that allows the conversion of the SAM to BAM. There is another more compressed version of the SAM file, called CRAM, which you can also create using `samtools` view.  However, we will not use that today.

**Note:** By using the `-t 4` parameter, we can take advantage of modern computers that allow multi-threading or parallelisation. This just means that the command can be broken up into 4 chunks and run in parallel, speeding up the process. Check you computers system settings, but you should be able to use at least 2 or 4 threads to run this alignment!

To find out information on your resulting alignment you can `samtools`:

```
samtools stats SRR2003569_chI.bam
```

This is basically the same as another command `samtools flagstat`, but it gives additional information.

### Questions

1. How many reads aligned to our genome?

2. How many reads aligned as a pair?

3. What information does `samtools stats` provide that `samtools flagstat` does not?

4. How many aligned as a "proper" pair? ..what the hell is a proper pair anyway??

# Viewing the alignments

A common tool used for viewing alignments is IGV browser.
We can open this just by entering `igv` in the terminal.

```
igv
```

Once you've open IGV, go to the `Genomes` menu & select `Load genome from file`.
Navigate to where you have `chrI.fa` and load this file.
Although this isn't the full genome, it will have everything we've aligned.

Now go to the `File` menu and select `Load from File` and navigate to your alignments.
Unfortunately you won't see anything until you zoom in.
This is so IGV doesn't hold the entire set of alignments in memory which would slow your computer to a stand-still.
Keep zooming in until some alignments appear then have a look around.

*What does all of the information mean when you hover over an alignment?*

We'll come back and have a look again after we've finished calling variants.
