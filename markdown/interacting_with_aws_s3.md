# Introduction

AWS S3 is a cloud storage solution for files. Unlike harddrives that are physically
connected to (or mounted on) a computer system, files need to be copied to / from
AWS S3 before they are accessible.

# The AWS command line tool

Amazon provides the `aws` command line too, which is already installed on your server.
The tool has lots of functionalities, but we will only be using the `aws s3` subcommand.

To check whether the `aws` command is available, try the following command:

```
aws s3
```

and

```
aws s3 help
```

# Listing the contents of folders on AWS S3

This is the top level of the `mbl.data` bucket:
```
aws s3 ls s3://mbl.data/
```

Within the bucket, you can navigate the folder structure, e.g. you can find all of the
data generated before the course in this location:

```
aws s3 ls s3://mbl.data/reads/pre_mbl/
```

**Challenge:** Can you get a list of all files available from *mouse* samples?


# Copying files to your server

To work with any of the files, we first have to copy them to your server's harddrive.
The `sync` (synchronize) command can copy complete directories for you.

For example, the following command will copy all FASTQ files for the mouse samples into
a the `reads` folder in your home directory.

```
mkdir -p ~/analysis
aws s3 sync --dryrun s3://mbl.data/reads/pre_mbl/mouse ~/reads
```

Wow, that was fast! Wait a minute - what does the `--dryrun` argument do?
Before you execute any `aws s3` command, it is good practice to preview what will happen 
by including the `--dryrun` argument. Things look good? Great, then rerun the command and
omit the `--dryrun` argument:

```
mkdir -p ~/analysis
aws s3 sync s3://mbl.data/reads/pre_mbl/mouse ~/analysis/reads
```

This command will now copy ~ 3 Gb of data from AWS S3 to your server.
