# Introduction

[Michael Love]()
has written the 
[tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)
R package to load `salmon` quantitation results into R.

All of the quantification result are found in the `quant.sf` output file.
We need to identify the path to all relevant files.

The following code snipped identifies all `quant.sf` files in the `analysis`
directory and extracts the sample name from the filename.

# Locating quantitation results

```
library(tximport)
library(stringr)

files <- dir("~/analysis", pattern = "quant.sf", recursive = TRUE, 
             full.names = TRUE)

names(files) <- str_match(files, pattern = ".*/(.*)_S\\d*_.*/quant.sf")[, 2]

files
```

# Mapping transcripts to genes

`salmon` quantifies gene expression at the transcript level. To obtain the expression
level at the gene level, e.g. the total abundance of all isoforms, we need to track which
transcripts are from the same gene locus.

For example, the quantitation file `quant.sf` for mouse samples contains rows for
transcripts with the following identifiers: `ENSMUST00000196221.1`,
`ENSMUST00000179664.1`, etc.

The `MUST` part of the identifier indicates that these are mouse *transcripts*. But are
they from the same gene? Or from different ones? 

This information is available e.g. via 
[Ensembl's biomart](http://uswest.ensembl.org/biomart)

```
library(edgeR)
library(stringr)
library(tximport)

# identify the paths to salmon output files in your `analysis` directory
files <- dir("~/analysis", pattern = "quant.sf", recursive = TRUE, 
             full.names = TRUE)
names(files) <- str_match(files, pattern = ".*/(.*)_S\\d*_.*/quant.sf")[, 2]

# read the gene annotation file
gene_anno <- read.delim("~/analysis/references/mouse/gene_annotations.txt.gz",
                        stringsAsFactors = FALSE)

# import the quantitation results into R 
txi <- tximport(files, type = "salmon", tx2gene = gene_anno[, c(2, 1)],
                countsFromAbundance = "lengthScaledTPM")
class(txi)
names(txi)
head(txi$counts)
head(row.names(txi$counts))

# create a gene annotation table
genes <- gene_anno[match(row.names(txi$counts), gene_anno$Gene.stable.ID), c(1, 3, 4)]
colnames(genes) <- c("gene_id", "gene_type", "gene_symbol")

# combine gene annotations & counts into a DGEList object
y <- DGEList(txi$counts, genes = genes)

```