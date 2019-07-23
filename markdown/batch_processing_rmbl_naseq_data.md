# Processing data from RNA-seq libraries generated at MBL

This document walks you through the steps required to 

- preprocess raw RNA-seq data obtained from the libraries *you* generated at MBL
from fish, fly and mouse samples
- create QC reports with FastQC and multiqc
- load the data into R for downstream analysis, e.g. differential expression
analysis.

## Before you start:

QC reports for the first (pre-MBL) datasets are now available in our AWS S3
bucket
[here](https://s3.amazonaws.com/mbl.data/index.html)
and
you can always load the annotated counts as the DGEList into R via the
[rmbl2019 R package](https://github.com/tomsing1/mbl2019/blob/master/rmbl2019/README.md).

Therefore, you can now clean the intermediate results from the previous analysis
from the`analysis` folder in your home directory on your server. 

**Note:** Please routinely save files and outputs that you would like to keep
long-term by exporting them to your own computer. The cloud servers will be
switched off at the end of the Genomics section of the course - and **ALL** of
the files and data on the servers will be deleted forever.


1. Move any files you would like to keep out of the `analysis` folder in your
home directory to another location. (Even better: save them to your laptop via
RStudio's `More -> Export` functionality.)

2. Delete  *all* files and folders within your `analysis` directory. For
example, delete 
- the `Snakefile`
- the `references` directory
- the `reads` directory
- the `quantitation` directory
- the `qc` directory
- the `trimmed` directory
- the `multiqc_report.html` file
- the `multiqc_data` directory

3. Execute the following code to retrieve all necessary files from the new
batch, alongside a fresh copy of the `references` indices and the `Snakefile`.
(These files haven't changed - but let's start with a clean slate.)

```
mkdir -p ~/analysis
aws s3 cp s3://mbl.data/scripts/Snakefile ~/analysis/
aws s3 sync s3://mbl.data/references ~/analysis/references
```

### To process the *26* FASTQ files from the `fly` samples:

```
aws s3 sync s3://mbl.data/reads/mbl/fly ~/analysis/reads
```

### To process the *120* FASTQ files from the `mouse` samples:

```
aws s3 sync s3://mbl.data/reads/mbl/mouse ~/analysis/reads
```

### To process the *67* FASTQ files from the `fish` samples:

```
aws s3 sync s3://mbl.data/reads/mbl/fish ~/analysis/reads
```

4. Edit the `SPECIES` argument in the `Snakefile`, e.g. if you downloaded the

- `fly` samples: `SPECIES = "fly"`
- `mouse` samples: `SPECIES = "mouse"`
- `fish` samples: `SPECIES = "fish"`

**Note:** You can edit the `Snakefile` either using `nano` or with `RStudio`.

That's it - you are ready to process the new FASTQ files. Just follow
the instructions in the next section.

## Processing the RNA-seq data

This document uses 
[snakemake]()
to coordinate an analysis workflow with the following steps:

1. Quality control of raw FASTQ files with `FastQC`
2. Adapter trimming with `skewer`
3. Mapping of the trimmed reads to the transcriptome with `salmon`
4. Collating a summary QC report with `multiqc`

## Prerequisites

Execution of this workflow requires that the following directory structure
is available in your `home` directory on your Amazon cloud server:


```
/home/ubuntu/analysis/
├── Snakefile
├── reads/
└── references/
```

The workflow will create additional directories `quantitation`, `qc` and 
`trimmed`.

The workflow will automatically process *all* FASTQ files in the `reads`
directory.

### Starting the workflow

To start the workflow, first start a new `screen` session on your server:

```
screen -S rnaseq
```

Within the screen, please move to the `analysis` directory, which should contain
the `Snakefile` you downloaded.

Start a dryrun of your workflow to see if things look ok:

```
cd ~/analysis
snakemake --dryrun
```

Snakemake will list the commands it intends to run, as well as their total
number. Do you see what you expected? Great - then start the real processing
run. Because your servers have two cores, you can use both of them by
specifying the `--cores 2` argument:

```
cd ~/analysis
snakemake --cores 2
```

Processing the files will take quite a while, so leave the screen by pressing
the following key combinations, one after the other, to detach the screen:

```
Ctrl-a
Ctrl-d
```

You can always re-attach the screen by typing:

```
screen -rd rnaseq
```

## Completing the pre-processing

When the workflow is finished, the `multiqc_report.html` file is generated
in your `analysis` directory. In other words, when you see the 
`multiqc_report.html` file, then you are done and can examine the quality of
the sequencing data. (Right-click the `multiqc_report.html` file in RStudio
to open it in your web browser.)

