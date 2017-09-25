# Introduction and Setup
{:.no_toc}

* TOC
{:toc}

# Introduction

## General Information

Thank you for your attendance & welcome to the *Introduction to NGS Data* Workshop.
This is an offering by the University of Adelaide, Bioinformatics Hub which is a centrally funded initiative from the Department of Vice-Chancellor (Research), with the aim of assisting & enabling researchers in their work.
Training workshops & seminars such as this one are an important part of this initiative.

Some additional resources run by the Bioinformatics Hub which may be of interest beyond today are:

- A University web-page at http://www.adelaide.edu.au/bioinformatics-hub/
- To be kept up to date on upcoming events and workshops, please join the internal Bioinformatics mailing list on http://list.adelaide.edu.au/mailman/listinfo/ bioinfo
- A Twitter account https://twitter.com/UofABioinfoHub/
- An active Slack team for discussing Bioinformatics questions with the local community. Slack teams do require an invitation to join, so please email the Hub on bioinf_hub@adelaide.edu.au to join the community. All are welcome.

Today’s workshop has been put together by Jimmy Breen and Steve Pederson.
We hope it will be useful in enabling you to continue and to advance your research.

## Course Summary

Next-generation sequencing (NGS) has become an important tool in assessing biological signal within an organism or population. Stemming from previous technologies that were costly and time-consuming to run, NGS platforms are relatively cheap and enable the investigation of the genome, transcriptome, methylome etc at extremely high resolution. The high-throughput of these machines also has unique challenges, and it is important that scientists are aware of the potential limitations of the platforms and the issues involved with the production of good quality data.

In this course, we will illustrate the fundamentals of Illumina NGS technology (the current market leader in the production of sequencing data), describe the affect that library preparation can have on downstream data, as well as running through a basic workflow used to analyse NGS data.

## Today's Instructors

Today your helpers will be

- Steve Pederson, Hien To, Alastair Ludington (Bioinformatics Hub)
- Jimmy Breen (Robinson Research Institute & Bioinformatics Hub)
- Terry Bertozzi (SA Museum)
- Rick Tearle, Lloyd Low (JS Davies Research Centre)

## Recommendations For Working

In the following pages, we strongly encourage you to manually type all commands.
The mistakes you will inevitably make will actually be important learning steps.
Additionally, in your work beyond today, you will probably not have any instructions to follow.
The experience of typing these commands will equip you for future work far better than if you simply copy & paste.

For today’s session, you will also be provided with red post-it notes.
**Please use these to signal whether you need help or not by placing them on your monitors**.
These are easy for instructors to spot so we can make our way over, although do be aware that there will be times when all instructors are busy.
This will be important as we all set our computers up in the following section.

## Computer Setup

For those running OSX or Ubuntu we will be able to run today's session on your own machines.
Many bioinformticians use these operating systems so the suite of tools is relatively mature.
In order to correctly configure your computer for today's session, please follow the links below.

- [Mac/OSX](../install/osxInstall)
- [Windows](../install/windowsInstall)
- [Ubuntu](../install/ubuntuInstall)

## Data For Today's Workshop

For those using a VM, the data will already be on your machines.
The rest of us will need to download then extract the data.
For convenience we recommend placing this in a directory called `~/WGS/rawData/fastq`

```
cd
mkdir -p WGS/rawData/fastq
```

**VM users will not need to execute the following command.**
All others will need to perform this todownload today's data.
```
wget -c "https://universityofadelaide.box.com/shared/static/cc0sgo2kya68zs2qu6r4ql2o38kits76.gz" -O "Intro-NGS-Sept-2017-files.tar.gz"
```

**All users will need to execute the following**.

```
tar -xzvf Intro-NGS-Sept-2017-files.tar.gz
mv *.fastq.gz WGS/rawData/fastq/
mv chrI.fa WGS
```


[Home](../)
