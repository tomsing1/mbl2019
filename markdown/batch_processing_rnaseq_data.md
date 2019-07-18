# Setting the stage

This document uses 
[snakemake]()
to coordinate an analysis workflow with the following steps:

1. Quality control of raw FASTQ files with `FastQC`
2. Adapter trimming with `skewer`
3. Mapping of the trimmed reads to the transcriptome with `salmon`
4. Collating a summary QC report with `multiqc`

## Prerequisites

Execution of this workflow requires that the following directory structure
is available in the user's home directory.

**Note: replace the word `SPECIES` with either `fly`, `mouse`, `planaria`,
`fish` or `worm`, as required.

```
~/analysis/
  reads/
  references/SPECIES/
```

The workflow will create additional directories `quantitation`, `qc` and 
`trimmed`.

The workflow will automatically process all FASTQ files in the `reads`
directory.

## Before you start:

1. If you have any older `quantitation`, `qc`, or `trimmed`, please move their
contents elsewhere (or delete them)!
2. Make sure your `reads` directory *only* contains FASTQ files that you want
to process right now!

## Getting the Snakefile

Please download the `Snakefile` from AWS S3 into your `analysis` directory.

```
aws s3 cp s3://mbl.data/scripts/Snakefile ~/analysis/
```

Open the `Snakefile` using nano and RStudio and check **change** `species`,
`FRAGMENT_LENGTH` and `FRAGMENT_SD` parameters at the top of the file. Save
the file if you made any changes.

## Getting the reference files

We have prepared `salmon` indices for the relevant species. Please download
the *right* index for your species of interest. 

**Note: replace the word `SPECIES` with either `fly`, `mouse`, `planaria`,
`fish` or `worm`, as required.

```
aws s3 sync s3://mbl.data/references/SPECIES ~/analysis/references/SPECIES
aws s3 cp s3://mbl.data/references/adapters.fa ~/analysis/references/
```

## Getting the FASTQ files

Please run the code **for the group you are assigned to** to place the
relevant FASTQ files into your `~/analysis/reads` folder.

**Note:** To make sure we only process the FASTQ files we are interested in,
we clean out the `~/analysis/reads` folder first.

### Fly, group 1

```
rm ~/analysis/reads/*.fastq.gz
aws s3 sync s3://mbl.data/reads/pre_mbl/fly_1 ~/analysis/reads
```

### Fly, group 2
```
rm ~/analysis/reads/*.fastq.gz
aws s3 sync s3://mbl.data/reads/pre_mbl/fly_2 ~/analysis/reads
```

### Fly, group 3
```
rm ~/analysis/reads/*.fastq.gz
aws s3 sync s3://mbl.data/reads/pre_mbl/fly_3 ~/analysis/reads
```

### Mouse
```
rm ~/analysis/reads/*.fastq.gz
aws s3 sync s3://mbl.data/reads/pre_mbl/mouse ~/analysis/reads
```

### Planaria
```
rm ~/analysis/reads/*.fastq.gz
aws s3 sync s3://mbl.data/reads/pre_mbl/planaria ~/analysis/reads
```

### Worm
```
rm ~/analysis/reads/*.fastq.gz
aws s3 sync s3://mbl.data/reads/pre_mbl/worm ~/analysis/reads
```

## Starting the workflow

To start the workflow, first start a new `screen` session on your server:

```
screen -S rnaseq
```

Within the screen, please move to the `analysis` directory, which should contain
the `Snakefile` you downloaded.

Start a dryrun of your workflow to see if things look ok:

```
cd ~/analysis
snakemake --dryrun --cores 4
```

Check the output - and then start the real run!

```
cd ~/analysis
snakemake --cores 4 2>&1 | tee snake.log
```
