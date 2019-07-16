# Introduction

Each sequencing library contains sequences of different lengths. Each of the DNA elements
in the library is flanked by two adapters, e.g. priming sites for the sequencing
reaction.

Sometimes, especially when the DNA fragments are very short, the sequencer reads through
the DNA sequence of interest and continues into the 3' adapter. It is good practice to
remove such spurious sequences.

Many different tools for adapter trimming exist, including e.g.

- [trim_galore_](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [skewer](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-182)

Here, we will use the `skewer` tool to trim known adapter sequences from our reads.

# Adapter sequences

What do adapters actually look like? Depending on the protocol used to create the
sequencing library, different adapter sequences may be present. 

Claire has prepared a small FASTA file with the adapters that were used to generate our
data. Let's download it from AWS S3:

```
mkdir -p ~/analysis/references/
aws s3 cp s3://mbl.data/references/adapters.fa ~/analysis/references/
```

Let's take a look at the file:

```
cd ~/analysis/references/
cat adapters.fa
```

# Running skewer

The following command trims adapters from the reads in the 
`F1_SC_S107_L002_R1_001.fastq.gz` FASTQ file.

```
mkdir -p trimmed
cd ~/analysis/trimmed
skewer -l 35 -t 4 -z \
    -o ~/analysis/trimmed/F1_SC \
    -x ~/analysis/references/adapters.fa \
    ~/analysis/reads/F1_SC_S107_L002_R1_001.fastq.gz
```

**Note:** 
- The command is very long, that's why we split it into multiple lines for better
readability. The `\` at the end of a line indicates that the command continues in the next
line.
- `skewer` automatically adds the `-trimmed.fastq.gz` suffix to the specified
output prefix.

What do all of those arguments mean (e.g. `-t 4`, `-l 35`, `-z` or `-x`)?
Check out the skewer help! Can you figure it out?

```
skewer --help
```

