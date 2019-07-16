# Introduction

Often the main goal of an RNA-seq experiment is quantitation of gene expression, either of
individual transcript isoforms or at the gene level.

Traditionally, reads were first aligned to the *genome* and then assigned to overlapping
annotated exons. The total counts for each gene (or isoform) were then obtained by
summing counts of all constituting exons.

For species with a well-annotated transcriptome, 
[Alignment-free methods](https://bioconductor.org/help/course-materials/2017/CSAMA/lectures/2-tuesday/lec07-alignmentfree-rnaseq.pdf)
of quantitation are both faster and more accurate.

Here, we will use the
[salmon]()
tool to map reads to the transcriptome and quantify the abundance of each transcript 
isoform.

# Creating a salmon transcriptome index

To map to the transcriptome we first need - well, a transcriptome. For most model
organisms, the
[Ensembl website](http://uswest.ensembl.org/index.html)
is a great source for biological sequences. Make sure to check their
[EnsemblMetazoa](https://metazoa.ensembl.org/species.html)
page for your species of interest as well.

For example, check out the
[Reference sequences](https://docs.google.com/document/d/10qZeS0D_hGm3o1J2T-KS83Q_6MkuBo3Xl_j2mlTsu7E/edit#bookmark=id.czvmpscnrvmw)
section on the course website to find the URL for the *Transcript sequences* for your
species of interest.

**Note:** For some species, we use sources other than Ensembl:

- For mouse (or human) data, we recommend using the
[Gencode](https://www.gencodegenes.org/)
annotations instead.
- The planarian transcriptome is not hosted by Ensembl, so we get the sequences from
[PlanMine](http://planmine.mpi-cbg.de/planmine/dataCategories.do)
instead.

Let's download the transcript sequences for the mouse:

```
mkdir -p ~/analysis/references/mouse
cd ~/analysis/references/mouse
wget ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
```

To make the transcriptome suitable for analysis with `salmon`, we first have to index it:

```
cd ~/analysis/references/mouse
salmon index \
    -p 4 \
    -k 21 \
    --gencode \
    -t Mus_musculus.GRCm38.cdna.all.fa.gz \
    -i salmon_index
```

**Note:**

- We added `--gencode` argument to inform `salmon index` that the sequence identifiers are
in a Gencode-specific format. This argument must **NOT** be used with transcriptomes from
other sources (e.g. Ensembl).

# Mapping and quantifying with Salmon

Now that we have the transcriptome index, we are ready to quantify our data:

```
mkdir ~/analysis/quantitation
cd ~/analysis/quantitation

salmon quant \
    -i ~/analysis/references/mouse/salmon_index \
    -l A \
    --fldMean 200 \
    --fldSD 100 \
    -r ~/analysis/trimmed/F1_SC-trimmed.fastq.gz \
    -p 8 \
    --validateMappings \
    -o F1_SC
```

The `salmon quant` command has a lot of options, you can see them either
[in the salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#description-of-important-options)
or with the following command: `salmon quant --help-reads`.

