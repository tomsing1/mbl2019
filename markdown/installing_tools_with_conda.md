# Introduction

Conda is a software repository that facilitate the installation of software.
It is not specific to biology or even science.

# Conda channels

Biology-focussed tools are available from a special `channel`, called `bioconda`.
In addition, the `conda-forge` channel also contains useful tools.

To make your conda installation aware of these channels, please issue the following
three commands in your terminal:

```
conda config --add channels bioconda
conda config --add channels conda-forge
```

Now your `conda` installation will also check the `bioconda` and 
`conda-forge` whenever you want to install a new tool.

# Installing tools

To install a single tool, e.g. the adapter trimming tool
[skewer](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-182)
you can use the `conda install` command:

```
conda install skewer
```

This command can also install multiple tools at the same time, e.g.

```
conda install fastqc multiqc salmon skewer snakemake 
```

