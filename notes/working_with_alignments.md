* TOC
{:toc}

# SAM/BAM files

## The SAM/BAM file format

Reads that have been to a reference are no longer stored in fastq format but are stored in either SAM or BAM format.
These two formats are virtually identical, however the SAM format is a text file which is easily readable to human eyes, whilst a BAM file is the same information converted to binary.
This conversion means that file sizes are smaller, and that computational processes can be performed more efficiently.
Typically, we work with BAM files as these provide gains in storage space & analytic speed.
The tools we use to inspect these files are provided in the package samtools, which has been installed on your VM.

The reads from the previous dataset which mapped to *chrI* of *C. elegans* are in the folder `~/WGS/trimmedData/fastq`

```
cd ~/WGS/trimmedData/fastq
ls
```

## Conversion to BAM format

The BAM format is much more convenient computationally, so we have converted our alignments into BAM format using `samtools view` during the alignment process.
SAM files are plain text, whilst BAM files are compressed and much easier for the computer to read/write.
As BAM files are in binary format they will look like gibberish if we try to read them directly.
Instead we can inspect them by using `samtools view` as in the line above.

There is also a header seection to each file which details the fasta files used in the alignments.
To view this header we use

```
samtools view SRR2003569_chI.bam | head
```

As this data can easily spill across lines, it might be helpful to maximise your terminal to see the complete line structure.

## The SAM/BAM data structure

If we understand what information is contained within a file, we can know what decisions to make as we progress with our analysis, so let's have a look at what the data structure is for a SAM/BAM file.
A SAM/BAM file is `tab-delimited`, which means that each field is separated by a tab, giving a data structure effectively consisting of columns (or fields).
In order, these are:

| 1 | QNAME | Query template/pair NAME |
| 2 | FLAG | bitwise FLAG |
| 3 | RNAME | Reference sequence NAME |
| 4 | POS | 1-based leftmost POSition/coordinate of clipped sequence |
| 5 | MAPQ | MAPping Quality (Phred-scaled) |
| 6 | CIGAR | extended CIGAR string |
| 7 | MRNM | Mate Reference sequence NaMe (`=' if same as RNAME) |
| 8 | MPOS | 1-based Mate POSistion |
| 9 | TLEN | inferred Template LENgth (insert size) |
| 10 | SEQ | query SEQuence on the same strand as the reference |
| 11 | QUAL | query QUALity (ASCII-33 gives the Phred base quality) |
| 12 | OPT | variable OPTional fields in the format TAG:VTYPE:VALUE |



Several of these fields contain useful information, so looking the the first few lines which we displayed above, you can see that these reads are mapped in pairs as consecutive entries in the QNAME field are often (but not always) identical.
Most of these fields are self-explanatory, but some require exploration in more detail.


## SAM Flags

These are quite useful pieces of information, but can be difficult at first look.
Head to http://broadinstitute.github.io/picard/explain-flags.html to see a helpful description.
The simplest way to understand these is that it is a bitwise system so that each description heading down the page increases ina binary fashion.
The first has value 1, the second has value 2, the third has value 4 & so on until you reach the final value of 2048.
The integer value contained in this file is the unique sum of whichever attributes the mapping has.
For example, if the read is paired \& mapped in a proper pair, but no other attributes are set, the flag field would contain the value 3.

#### Questions
{:.no_toc}


1. *What value could a flag take if the read was 1 - paired; 2 - mapped in a proper pair; 3 - it was the first in the pair \& 4 - the alignment was a supplementary alignment.*
2. *Some common values in the bam file are 99, 147 & 145. Look up the meanings of these values.*



Things can easily begin to confuse people once you start searching for specific flags, but if you remember that each attribute is like an individual flag that is either on or off (i.e. it is binary).
If you searched for flags with the value 1, you wouldn't obtain the alignments with the exact value 1, rather you would obtain the alignments for which the first flag is set & these can take a range of values.


Let's try this using the command `samtools view` with the option `-f N` to include reads with a flag set and the option `-F N` to exclude reads with a specific flag set.
Let's get the first few reads which are mapped in a proper pair, so the flag `2` will be set.

```
samtools view -f 2 SRR2003569_chI.bam | head
```

Note that none of the flags actually have the value 2, but if you typed the values 99, 147 or 163 into the webpage, you'll see that this flag is set for all of these values.
Similarly if we wanted to extract only the reads which are NOT mapped in a proper pair we would change the option to a upper-case F.

```
samtools view -F 2 SRR2003569_chI.bam | head
```

Again, try entering a few of these sample values into the webpage and you will see that this flag is not set for any of these values.

This can be a very helpful tool for extract subsets of your aligned reads.
For example, we can create a new BAM file with only the reads which were aligned in a proper pair by entering the following command.

```
samtools view -f 2 -bo SRR2003569_chI.bam
ls -lh
```

You can pull out highly specific combinations of alignments should you so choose


## CIGAR strings

These give useful information about the type of alignment that has been performed on the read.
In the first few reads we called up earlier, most had the value `..M` where `..` is some number.
These are the perfect Matches, where the sequence has aligned exactly.
The other abbreviations in common use are I (insertion), D (deletion) & S(substitution).

What is the interpretation of the first `CIGAR` string in your set of alignments.
