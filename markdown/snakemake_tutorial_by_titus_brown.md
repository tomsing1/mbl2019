# Running snakemake!

This material is a slightly modified version of the tutorial published
by Titus Brown [here](https://hackmd.io/7k6JKE07Q4aCgyNmKQJ8Iw?view#Running-snakemake).

## Getting started - your first Snakefile

Create a new text file (File, New File, Text file) and write:

```
rule fastqc_a_file:
  shell:
    "fastqc reads/F1_SC_S107_L002_R1_001.fastq.gz"
```

(I suggest copy/pasting this into RStudio.)

Then save it as a file named Snakefile.

Now, run snakemake:

```
snakemake
```

and you should see:

```
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 1
Rules claiming more threads will be scaled down.
...
```

and there will be two new files,

Points to make:

- the snakemake configuration file is by default called Snakefile

## Updating the Snakefile to track inputs and outputs

At the moment this is basically just a shell script with extra syntax… what’s the point?

Well, shell scripts - and this snakefile, too - will rerun the command every time you run the file, even if there’s no reason to do so because the file hasn’t changed.

Digression: This is particularly important for large or long workflows, where you’re dealing with 10s to 100s of files that may take hours to days to process! It can be hard to figure out which files to rerun, but (spoiler alert) snakemake can really help you do this!

It’s hard to track this kind of thing in a shell script - I usually just comment out the lines I don’t want run, or break my commands up into multiple shell scripts so they don’t take so long - but with snakemake, you can annotate the rule with input and output files!

Change your snakefile to look like this:

```
rule fastqc_a_file:
  input:
    "reads/F1_SC_S107_L002_R1_001.fastq.gz"
  output:
    "reads/F1_SC_S107_L002_R1_001_fastqc.html",
    "reads/F1_SC_S107_L002_R1_001_fastqc.zip"
  shell:
    "fastqc reads/F1_SC_S107_L002_R1_001.fastq.gz"
```

here, we’ve annotated the rule with the required
input file, as well as the expected output files.

Question: how do we know what the output files are?

Now run:

```
snakemake
```

and you should see:

```
Building DAG of jobs...
Nothing to be done.
```

What happened??

snakemake looked at the file, saw that the output files existed, and figured out that it didn’t need to do anything!

## Forcibly re-running things
You can tell snakemake to run the rule no matter what with -f:

```
snakemake -f
```

You can also remove an output file and it will automatically re-run:

```
rm reads/html
snakemake
```

note that you don’t need to remove all the output files to rerun a command - just remove one of them.

You can also update the timestamp on an input file, and snakemake will figure out that the output file is older than the input file, and rerun things.

```
touch reads/*.fastq.gz
snakemake
```

This will become important later :)

## Multiple rules

We have many FASTQ files. Can Snakemake help us process multiple files at once?
Let's copy an additional FASTQ file into our `reads` directory:

```
aws s3 cp s3://mbl.data/reads/pre_mbl/mouse/M1_TGR_S113_L002_R1_001.fastq.gz ~/analysis/reads/
```

Let’s add a rule to run fastqc on a second file:

```
rule fastqc_a_file:
  input:
    "reads/F1_SC_S107_L002_R1_001.fastq.gz"
  output:
    "reads/F1_SC_S107_L002_R1_001_fastqc.html",
    "reads/F1_SC_S107_L002_R1_001_fastqc.zip"
  shell:
    "fastqc reads/F1_SC_S107_L002_R1_001.fastq.gz"

rule fastqc_a_file2:
  input:
    "reads/M1_TGR_S113_L002_R1_001.fastq.gz"
  output:
    "reads/M1_TGR_S113_L002_R1_001_fastqc.html",
    "reads/M1_TGR_S113_L002_R1_001_fastqc.zip"
  shell:
    "fastqc reads/M1_TGR_S113_L002_R1_001.fastq.gz"
```

Now, if you run this, the Right Thing won’t happen: snakemake will do nothing. Why?

Well, snakemake only runs the first rule in a Snakefile, by default. You can give a rule name on the command line, if you like, or you can tell snakemake what output file(s) you want. Let’s do the latter:

```
snakemake reads/F1_SC_S107_L002_R1_001_fastqc.html reads/M1_TGR_S113_L002_R1_001_fastqc.html
```

and now you should see the second fastqc command run, with the appropriate output files!

Note that snakemake only runs the second rule, because it looks at the output files and sees that the first file you wanted, 0Hour_001_1_fastqc.html already exists!

Points to make:

- this is pretty long compared to the same shell script…
- specifying which file or rule you want is kind of annoying…

## A first refactoring: adding a better default rule
Let’s start refactoring (cleaning up) this Snakefile.

First, let’s add a rule at the top:

```
rule all:
  input:
    "reads/F1_SC_S107_L002_R1_001_fastqc.html",
    "reads/M1_TGR_S113_L002_R1_001_fastqc.html"

rule fastqc_a_file:
  input:
    "reads/F1_SC_S107_L002_R1_001.fastq.gz"
  output:
    "reads/F1_SC_S107_L002_R1_001_fastqc.html",
    "reads/F1_SC_S107_L002_R1_001_fastqc.zip"
  shell:
    "fastqc reads/F1_SC_S107_L002_R1_001.fastq.gz"

rule fastqc_a_file2:
  input:
    "reads/M1_TGR_S113_L002_R1_001.fastq.gz"
  output:
    "reads/M1_TGR_S113_L002_R1_001_fastqc.html",
    "reads/M1_TGR_S113_L002_R1_001_fastqc.zip"
  shell:
    "fastqc reads/M1_TGR_S113_L002_R1_001.fastq.gz"
```

this rule, by convention called all, is a default rule that produces all the output files. But it’s a bit weird! It’s all input, and no output!

This is a blank rule that gathers together all of the various files you want produced, and says “hey, snakemake, I depend on all of these files for my input - make them for me!” And then, once those files are all there, it …does nothing.

Yep, this is perfectly legal in snakemake, and it’s one way to make your life easier.

Note that `snakemake -f` no longer works properly, because `-f` only forces rerunning a single rule. 

## A second refactoring: doing a bit of templating
There’s a lot of repetition in each of these rules. Let’s collapse it down a little bit by replacing the filename in the fastqc command with a magic variable, {input}.

```
rule all:
  input:
    "reads/F1_SC_S107_L002_R1_001_fastqc.html",
    "reads/M1_TGR_S113_L002_R1_001_fastqc.html"

rule fastqc_a_file:
  input:
    "reads/F1_SC_S107_L002_R1_001.fastq.gz"
  output:
    "reads/F1_SC_S107_L002_R1_001_fastqc.html",
    "reads/F1_SC_S107_L002_R1_001_fastqc.zip"
  shell:
    "fastqc {input}"

rule fastqc_a_file2:
  input:
    "reads/M1_TGR_S113_L002_R1_001.fastq.gz"
  output:
    "reads/M1_TGR_S113_L002_R1_001_fastqc.html",
    "reads/M1_TGR_S113_L002_R1_001_fastqc.zip"
  shell:
    "fastqc {input}"
```

This all works as before, but now the rule is a bit more generic and will work with any input file. Sort of.

## Refactoring 3: templating output files, too
What do I mean, sort of?

Well, the output filenames ALSO depend on the input file names in some way - specifically, fastqc replace part of the filename with _fastqc.html and _fastqc.zip to make its two output files.

Let’s rewrite the rule using some snakemake pattern matching:

```
rule all:
  input:
    "reads/F1_SC_S107_L002_R1_001_fastqc.html",
    "reads/M1_TGR_S113_L002_R1_001_fastqc.html"

rule fastqc_a_file:
  input:
    "{filename}.fastq.gz"
  output:
    "{filename}_fastqc.html",
    "{filename}_fastqc.zip"
  shell:
    "fastqc {input}"

rule fastqc_a_file2:
  input:
    "{filename}.fastq.gz"
  output:
    "{filename}_fastqc.html",
    "{filename}_fastqc.zip"
  shell:
    "fastqc {input}"
```

What we’ve done here is tell snakemake that anytime we say we want a file that ends with _fastqc.html, it should look for a file that ends in .fastq.gz and then run fastqc on it.

Try running this:

```
snakemake
```

Oh no! We get a `AmbiguousRuleException`! What’s going on?

Well, if you look at the rule above, we’ve given snakemake two different rules to produce the same file(s)! fastqc_a_file and fastqc_a_file2 are now identical rules! snakemake doesn’t like that.

Let’s remove one, to get a trimmer, leaner, and above all functional snakefile:

```
rule all:
  input:
    "reads/F1_SC_S107_L002_R1_001_fastqc.html",
    "reads/M1_TGR_S113_L002_R1_001_fastqc.html"

rule fastqc_a_file:
  input:
    "{filename}.fastq.gz"
  output:
    "{filename}_fastqc.html",
    "{filename}_fastqc.zip"
  shell:
    "fastqc {input}"
```

and THAT should now work just fine!

## Rerunning snakemake
Note you can just run snakemake whenever you want. It won’t do anything unless something’s changed.

```
snakemake
```
It’s actually kind of soothing…

```
snakemake
snakemake
snakemake
```

## Building out the workflow

So, we’ve gotten `fastqc` sorted out. What’s next?

Let’s add in a new rule - multiqc, to summarize our fastqc results.

multiqc takes a directory name under which there are one or more fastqc reports, and then summarizes them.

Running it on the command line,

```
multiqc .
```

you can see that it creates two outputs, multiqc_report.html and the directory multiqc_data/ which contains a bunch of files. Let’s create a snakemake rule for this; add:

```
rule run_multiqc:
  output:
    "multiqc_report.html",
    directory("multiqc_data")
  shell:
    "multiqc ."
```

to the bottom of the file. (Note, you need to tell snakemake if an output is a directory.)

Now run it:

```
snakemake run_multiqc
```

This …doesn’t really do what we want, for a few reasons.

- First of all, you have to specify the rule name or else snakemake doesn’t run anything. How do we fix this??

- Second of all, `multiqc_report.html` already exists, so snakemake doesn’t run the rule. How do we actually test the rule??

- Third of all, the multiqc rule has no input dependencies. How do we specify them??

Let’s fix the first two things first:

-add multiqc_report.html to the inputs for the first all.
- then remove multiqc_report.html and re-run snakemake.

Your snakefile should look like:

```
rule all:
  input:
    "reads/F1_SC_S107_L002_R1_001_fastqc.html",
    "reads/M1_TGR_S113_L002_R1_001_fastqc.html"
    "multiqc_report.html"

rule fastqc_a_file:
  input:
    "{filename}.fastq.gz"
  output:
    "{filename}_fastqc.html",
    "{filename}_fastqc.zip"
  shell:
    "fastqc {input}"

rule run_multiqc:
  output:
    "multiqc_report.html",
    directory("multiqc_data")
  shell:
    "multiqc ."
```

Yay, that seems to work!

Points to make:

- other than the first rule, rules can be in any order
- the rule name doesn’t really matter, it’s mostly for debugging. It just needs to be “boring” (text, underscores, etc. only)

## Providing input files explicitly the multiqc rule
The third problem, that multiqc doesn’t have any input dependencies, is a bit harder to fix.

(Why do we want to fix this? Well, this is how snakemake tracks “out of date” files - if we don’t specify input dependencies, then we may update one of the fastqc results that multiqc uses, but snakemake won’t re-run multiqc on it, and our multiqc results will be out of date.)

To use variables, let’s make a Python list at the very top, containing all of our expected output files from fastqc:

```
    "reads/F1_SC_S107_L002_R1_001_fastqc.html"
    "reads/M1_TGR_S113_L002_R1_001_fastqc.html"

fastqc_output = ["reads/F1_SC_S107_L002_R1_001_fastqc.html", "reads/F1_SC_S107_L002_R1_001_fastqc.html"]
```

and modify the all and multiqc rules to contain this list. The final snakefile looks like this:

```
fastqc_output = ["reads/F1_SC_S107_L002_R1_001_fastqc.html", "reads/F1_SC_S107_L002_R1_001_fastqc.html"]

rule all:
  input:
    fastqc_output,
    "multiqc_report.html"

rule fastqc_a_file:
  input:
    "{filename}.fastq.gz"
  output:
    "{filename}_fastqc.html",
    "{filename}_fastqc.zip"
  shell:
    "fastqc {input}"

rule run_multiqc:
  input:
    fastqc_output
  output:
    "multiqc_report.html",
    directory("multiqc_data")
  shell:
    "multiqc ."
```

Points to make:

- quoted strings are …strings (filename, usually)
- you can also use Python to create, manipulate, etc. filenames or lists of them, and it will work fine!
- this includes loading in filenames from spreadsheets, etc. - more on that later.

## Refactoring this to make it slightly more concise –
We’ve got one more redundancy in this file - the fastqc_output is listed in the all rule, but you don’t need it there! Why?

Well, multiqc_report.html is already in the all rule, and the multiqc rule depends on fastqc_output, so fastqc_output already needs to be created to satisfy the all rule, so… specifying it in the all rule is redundant! And you can remove it!

(It’s not a big deal and I usually leave it in. But I wanted to talk about dependencies!)

The Snakefile now looks like this:

```
fastqc_output = ["reads/F1_SC_S107_L002_R1_001_fastqc.html", "reads/F1_SC_S107_L002_R1_001_fastqc.html"]

rule all:
  input:
    "multiqc_report.html"

rule fastqc_a_file:
  input:
    "{filename}.fastq.gz"
  output:
    "{filename}_fastqc.html",
    "{filename}_fastqc.zip"
  shell:
    "fastqc {input}"

rule run_multiqc:
  input:
    fastqc_output
  output:
    "multiqc_report.html",
    directory("multiqc_data")
  shell:
    "multiqc ."
```

and we can rerun it from scratch by doing:

```
rm data/*.html
snakemake
```

## Making a clean rule
It’s kind of annoying to have to delete things explicitly. Snakemake should take care of that for us. Let’s add a new rule, clean, that forces rerunning of things –

```
rule clean:
  shell:
    "rm -f {fastqc_output} multiqc_report.html"
```

and now try rerunning things:

```
snakemake -p clean
snakemake
```

A few things to point out –

- Here we see the use of variables inside a shell command, again - {fastqc_output} means "replace the thing in the curly quotes with the Python values in fastqc_output.
- We’re using snakemake -p to get a printout of the commands that are run.
- What’s particularly nice about the clean rule (the name is a convention, not a requirement) is that you only need to keep track of the expected output files in one or two places - the all rule, and the clean rule.

## Digression: running things in parallel
So, we’ve put all this work into making this snakefile with its input rules and its output rules… and there are a lot of advantages to our current approach already! Let’s list a few of them –

- we’ve completely automated our analysis!
- we can easily add new data files into fastqc and multiqc!
- we can rerun things easily, and (even better) by default only things that need to be run will be run.
- the snakefile is actually pretty reusable - we could drop this into a new project, and, with little effort, run all of these things on new data!
but there’s even more practical value in this, too – because we’ve given snakemake a recipe, rather than a script, we can now run these commands in parallel, too!

Try:

```
snakemake clean
snakemake -j 2
```

this will run up to four things in parallel!


More advanced snakemake
## Dry run
You can use snakemake -n to see what would be run, without actually running it! This is called a “dry run”. This is useful when you have really big compute.

```
snakemake clean
snakemake -n
```

## Outputting the entire workflow diagram
You can visualize your workflow like so!

```
snakemake --dag | dot -Tpng > dag.png
```

# Additional material

- [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html)
- [Snakemake tutorial by Titus Brown](https://hackmd.io/7k6JKE07Q4aCgyNmKQJ8Iw?view#Running-snakemake)
- [Genome Analysis Workshop](https://molb7621.github.io/workshop/Classes/snakemake-tutorial.html)
