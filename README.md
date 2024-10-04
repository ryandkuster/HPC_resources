## Bioinformatics exercises for ISAAC-NG

### Set up practice directory

To begin, we'll make a simplified project directory. There are many ways to organize bioinformatics projects that are better, but this is an okay starting point.
```bash
cd $SCRATCHDIR/
mkdir bioinfo_practice
cd bioinfo_practice

mkdir data
mkdir data/reference
mkdir data/derived
mkdir data/raw
mkdir logs
mkdir results
```

### Get sequence data

This directory is empty. Let's use `wget` to directly download data from a known URL.

Because we'll be transferring data, it's best practice to use a data transfer node (DTN). We can just ssh into one.
```bash
ssh dtn1
```

When we ssh into another node, it is common to end up in a home directory. Let's change all the way back into our practice dir and download some files.

We'll be downloading an _E. coli_ reference genome and _S. aureus_ protein sequences from NCBI ftp.
```bash
cd $SCRATCHDIR/bioinfo_practice/data/reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_protein.faa.gz
```

Let's exit the dtn.
```bash
exit
```

Let's move back to our original project directory, then uncompress the gzipped fasta reference sequence.
```bash
cd $SCRATCHDIR/bioinfo_practice
gunzip data/reference/GCF_000005845.2_ASM584v2_genomic.fna.gz
```

### Analyzing sequence data with BLAST

Let's see how many of the protein sequences from _S. aureus_ map closely to the _E. coli_ reference genome.

Previously, we accessed the dtn1 using `ssh` to transfer data on a node meant for data transfer. Now that we're interested in doing computing that may use more resources than typical login node activity, we need to use a compute node.

Let's request a compute node using a `SLURM` command called `srun`. This will let us interact directly with our job instead of sending it into the background to run.
```bash
srun --account ACF-UTK0011 --partition=short --qos=short --nodes=1 --ntasks=2 --mem=8G --time=0-03:00:00 --pty /bin/bash
```
> [!NOTE]
> _The node this command outputs will be unique to each user_

When we ssh into a node, it typically sets us in our home directory. Let's move back to our practice directory so we can do actual work:
```bash
cd $SCRATCHDIR/bioinfo_practice
```

We're interested in building a database from the _E. coli_ reference genome (GCF_000005845). Using `module load` we can make blast-plus available in our environment.

The command `makeblastdb` should be available. Before running our command we can use the `-version` flag to confirm we have blast-plus loaded.
```bash
module load blast-plus/2.12.0
makeblastdb -version
```

Now make the database. Any output from the tool that would typically print to the screen (stdout) will be redirected with `>` to a log file.

By using `2>&1` after the redirect isn't always needed, but here it implies that we also want to write error messages (stderr) to that same log file.
```bash
makeblastdb -in data/reference/GCF_000005845.2_ASM584v2_genomic.fna -dbtype nucl -out data/derived/e_coli_db > logs/makeblastdb.log 2>&1
```

```bash
gunzip data/reference/GCF_000013425.1_ASM1342v1_protein.faa.gz
tblastn -query data/reference/GCF_000013425.1_ASM1342v1_protein.faa -evalue .1 -db data/derived/e_coli_db -out results/e_coli_s_aureus_tblastn_results.txt -outfmt 6 > logs/tblastn.log 2>&1
```

### Downloading and using pre-compiled binary software

At times, software may not be available using `module`. There are many options that require special configuration such as `conda` or `singularity`, but commonly a simple solution for well-supported software is to download a pre-compiled binary.

Again, we'll use the dtn since we'll be downloading a potentially large amount of data.

```bash
ssh dtn1
mkdir -p ~/bin
wget -P ~/bin/ https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xvzf ~/bin/sratoolkit.current-ubuntu64.tar.gz -C ~/bin
```

Let's download some _S. aureus_ short read RNASeq data.
```bash
~/bin/sratoolkit.3.1.1-ubuntu64/bin/fasterq-dump --split-files SRR29481819 -O $SCRATCHDIR/bioinfo_practice/data/raw
```

Assuming you've previously connected to the `salloc` node, then the dtn1 node, we're getting somewhat nested now!

Let's exit the dtn1 node and see what node we're in:
```bash
exit
```

If we're still in our `srun` allocated node, we can keep working!

Just as we previously loaded the BLAST-PLUS module, we can unload it.
```bash
module unload blast-plus/2.12.0
```

### Creating fastq QC statistics and transferring locally with scp

Let's load fastqc.
```bash
module load fastqc/0.11.9
fastqc -o $SCRATCHDIR/bioinfo_practice/results $SCRATCHDIR/bioinfo_practice/data/raw/SRR29481819_1.fastq $SCRATCHDIR/bioinfo_practice/data/raw/SRR29481819_2.fastq
```
Finally, to transfer the html files for visualization on your local device, let's use secure copy (scp).

> [!NOTE]
> _You'll need to have a working terminal set up on your local device. If you're having trouble don't worry, these files can also be found using OnDemand._

From *your* device, go to a directory where you'd like to download the fastqc output (e.g., Downloads).
```bash
scp -r <your usename>@dtn2.isaac.utk.edu:/lustre/isaac/scratch/<your usename>/bioinfo_practice/results/\*html .
```

You'll be prompted again for your password.

> [!NOTE]
> _When using scp to transfer files to zsh (macOS), you need to escape all asterisks with backslashes_

That completes the exercises, don't forget to exit the `salloc` node by using `exit`. Exiting twice is typically necessary to free up the node.

