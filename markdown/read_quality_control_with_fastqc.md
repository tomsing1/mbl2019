# Introduction

Before starting the analysis, it is a good idea to check the quality of the raw data
that we received from the sequencing center.

The
[FastQC tool](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
offers a one-stop-shop to create an html file with many different QC metrics.

# Running FastQC

To create a QC report for a single sample, simply point the `fastqc` command at the
location of the file on your server.

```
mkdir -p ~/analysis
cd ~/analysis
fastqc reads/F1_SC_S107_L002_R1_001.fastq.gz
```

Great! Now we have two new files in the `~/analysis/reads` folder:

- F1_SC_S107_L002_R1_001_fastqc.html
- F1_SC_S107_L002_R1_001_fastqc.zip

Let's navigate to the HTML file and take a look.

# Creating reports in a separate directory

The first FastQC report we generated was created in the `reads` subdirectory. That's a
bit messy. 

To instruct FastQC to create its output in a different directory, we can
add the `-o` argument to specify the output location.

**Note:** The output directory has to exist!

```
mkdir -p ~/analysis/qc
cd ~/analysis
fastqc -o qc reads/F1_SC_S107_L002_R1_001.fastq.gz
```

Now that we have another copy of the QC report, please delete the 

- F1_SC_S107_L002_R1_001_fastqc.html
- F1_SC_S107_L002_R1_001_fastqc.zip

files from the `~/analysis/reads` directory.

**Warning:** Take care not to delete any of the FASTQ files in the same location!

