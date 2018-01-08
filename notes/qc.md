* TOC
{:toc}

# Quality Control

## Using fastqc

A common tool for checking the quality of a fastq file is the program fastqc.
As with all programs on the command line, we need to see how it works before we use it.
The following command will open the help file in the less pager which we used earlier.
To navigate through the file, use the `<spacebar>` to move forward a page, `<b>` to move back a page & `<q>` to exit the manual.

```
fastqc -h | less
```

Fastqc will create an html report on each file you submit which can be opened from any web browser, such as firefox.
As seen in the help page, fastqc can be run from the command line or from a graphic user interface (GUI).
Using a GUI is generally intuitive so today we will look at the command line usage, as that will give you more flexibility & options going forward.
Some important options for the command can be seen in the manual.
As you will see in the manual, setting the `-h` option as above will call the help page.
Look up the following options to find what they mean.

| Option | Usage |
|:------ |:------|
| -o     |       |
| -t     |       |


As we have two files, we will first need to create the output directory, then we can run fastqc using 2 threads which will ensure the files are processed in parallel.
This can be much quicker when dealing with large experiments.

```
cd ~/WGS/rawData/
mkdir FastQC
cd fastq
fastqc -o ../FastQC -t 2 *gz
```

It’s probably a good idea to scribble a note next to each line if you didn’t understand what you did.
If you haven’t seen the command `mkdir` before, check the help page `man mkdir`.

The above command gave both files to `fastqc`, told it where to write the output (`-o  ̃FastQC`) & requested two threads (`-t 2`).
The reports are in the html files, with all of the plots & additional data stored in the zip files.

To look at the QC report for each file, we can use `firefox`.

```
cd ~/WGS/rawData/FastQC
ls -lh
```

The reports are in `html` files, which may be in the `FastQC` directory, or may be in the directories for the individual files, (depending on your version of FastQC).
Find the `html` files and open using your favourite browser.
If you have FileZilla installed, you can download these to your local machine and open them from there.
Otherwise open them in Ubuntu using the following command.
(The best browser for those on the VMs is firefox).

```
firefox *html
```

The left hand menu contains a series of click-able links to navigate through the report, with a quick guideline about each section given as a tick, cross of exclamation mark.
Two hints which may make your inspection of these files easier are:

1. To zoom out in firefox use the shortcut `Ctrl-`. Reset using `Ctrl0` and zoom in using `Ctrl+`
2. You can open these directly from a traditional directory view by double clicking on the .html file.

If your terminal seems busy after you close firefox, use the `Ctrl C` shortcut to stop whatever is keeping it busy.

#### Questions
{:.no_toc}

1. *How many sequences are there in both files?*
2. *How long are the sequences in these files?*

## Interpreting the FASTQC Report

As we work through the QC reports we will develop a series of criteria for filtering and cleaning up our files.
There is usually no perfect solution, we just have to make the best decisions we can based on the information we have.
Some sections will prove more informative than others, and some will only be helpful if we are drilling deeply into our data.
Firstly we’ll just look at a selection of the plots.
We’ll investigate some of the others with some ‘bad’ data later.

### Per Base Sequence Quality
{:.no_toc}

Both of the files should be open in firefox in separate tabs.
Perform the following steps on both files.
Click on the `Per base sequence quality` hyperlink on the left of the page & you will see a boxplot of the QC score distributions for every position in the read.
This is the main plot that bioinformaticians will look at for making informed decisions about later stages of the analysis.

*What do you notice about the QC scores as you progress through the read?*

We will deal with trimming the reads in a later section, but start to think about what you should do to the reads to ensure the highest quality in your final alignment & analysis.

**Per Tile Sequence Quality**
This section just gives a quick visualisation about any physical effects on sequence quality due to the tile within the each flowcell.
For the first file, you will notice an even breakdown in the quality of sequences near the end of the reads across all tiles.
In our second QC report, you will notice a poor quality around the 25th base in the 2nd (or 3rd) tile.
Generally, this would only be of note if drilling deeply to remove data from tiles with notable problems.
Most of the time we don’t factor in spatial effects, unless alternative approaches fail to address the issues we are dealing with.

**Per Sequence Quality Scores** This is just the distribution of average quality scores for each sequence.
There’s not much of note for us to see here.

**Per Base Sequence Content** This will often show artefacts from barcode sequences or adapters early in the reads, before stabilising to show a relatively even distribution of the bases.

**Sequence Length Distribution** This shows the distributions of sequence lengths in our data. Here we have sequences that are all the same lengths, however if the length of your reads is vital (e.g. smallRNA data), then this can also be an informative plot.

**Sequence Duplication Levels** This plot shows about what you’d expect from a typical NGS experiment.
There are a few duplicated sequences (rRNA, highly expressed genes etc.) and lots of unique sequences represented the diverse transcriptome.
This is only calculated on a small sample of the library for computational efficiency and is just to give a rough guide if anything unusual stands out.

**Overrepresented Sequences** Here we can see any sequence which are more abundant than would be expected. Sometimes you'll see sequences here that match the adapters used, or you may see highly expressed genes here.

**Adapter Content** This can give a good guide as to our true fragment lengths. If we have read lengths which are longer than our original DNA/RNA fragments (i.e. inserts) then the sequencing will run into the adapters.
If you have used custom adapters, you may need to supply them to `FastQC` as this only searches for common adapter squences.

**Kmer Content** Statistically over-represented `k`-mers can be seen here & often they will overlap.
In our first plot, the black & blue `k`-mers are the same motif, just shifted along one base.
No information is given as to the source of these sequences, and you would expect to see barcode sequences or motifs that correspond to any digestion protocols here.

## Some More Example Reports

Let’s head to another sample plot at the [FastQC homepage](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)

**Per Base Sequence Quality** Looking at the first plot, we can clearly see this data is
not as high quality as the one we have been exploring ourselves.

**Per Tile Sequence Quality** Some physical artefacts are visible & some tiles seem to
be consistently lower quality. Whichever approach we take to cleaning the data will more
than likely account for any of these artefacts. Sometimes it’s just helpful to know where a
problem has arisen.

**Overrepresented Sequences** Head to this section of the report & scan down the
list. Unlike our sample data, there seem to be a lot of enriched sequences of unknown
origin. There is one hit to an Illumina adaptor sequence, so we know at least one of the
contaminants in the data. Note that some of these sequences are the same as others on
the list, just shifted one or two base pairs. A possible source of this may have been non
random fragmentation.

**Kmer Content**
*Do you notice anything unusual about this plot?*


Interpreting the various sections of the report can take time & experience.
A description of each of the sections [is available from the fastqc authors](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/) which can be very helpful as you're finding your way.

Another interesting report is available at http://www.bioinformatics.babraham.ac.uk/projects/fastqc/RNA-Seq_fastqc.html.
Whilst the quality scores generally look pretty good for this one, see if you can find a point of interest in this data.
This is a good example, of why just skimming the first plot may not be such a good idea.

## Working With Larger Datasets

In our dataset of two samples it is quite easy to think about the whole experiment &
assess the overall quality.
*What about if we had 100 samples?*
Each .zip archive contains text files with the information which can easily be parsed into an overall summary.
We could write a script to extract this information if we had the time.
However, the Bioinformatics Hub has been writing an `R` package to help with this, which is available from https://github.com/UofABioinformaticsHub/ngsReports.

We're hoping to publish this soon and using the package is beyond the scope of today.
However, we've included a [sample report](../data/ngsReports_Fastqc) of a dataset summarised using heatmaps.
This is simply the default report produced and the package is capable of exploring large datasets relatively easily.
Have a look at [this report](../data/ngsReports_Fastqc) and see if you can understand any of the plots.
Call an instructor over if you have any questions.
