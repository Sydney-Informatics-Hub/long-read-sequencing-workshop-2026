# Workshop development notes

## Data

### Downloading FASTQ files

#### Prussing et al., 2020

The following command was used to download the FASTQs from the [Prussing et al. 2020 publication](https://doi.org/10.3389/fmicb.2020.02007):

```bash
for ACC in SRR11909877 SRR11909879 SRR11909881 SRR11909883 SRR11909885; do
    fasterq-dump $ACC
done
```

This downloads the following FASTQs:

```console
SRR11909877.fastq
SRR11909879.fastq
SRR11909881.fastq
SRR11909883.fastq
SRR11909885.fastq
```

**Note** that this command relies on the `fasterq-dump` tool, which is part of the `sra-tools` package. We used the `ncbi/sra-tools:3.4.1` Docker image to run the above command.

#### Roberts et al., 2023

The following commands were used to download the FASTQs from the [Roberts et al., 2023 publication](https://doi.org/10.1099/mgen.0.001048):

```bash
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR828/003/ERR8282753/ERR8282753.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR828/002/ERR8282752/ERR8282752.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR828/001/ERR8282751/ERR8282751.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR828/002/ERR8282742/ERR8282742.fastq.gz
```

#### Ryan Wick's pre- and post-trimmed bacterial ONT FASTQs

Ryan Wick has made several bacterial ONT samples available from his paper [Completing bacterial genome assemblies with multiplex MinION sequencing](https://doi.org/10.1099/mgen.0.000132). The [GitHub repository for the paper](https://github.com/rrwick/Bacterial-genome-assemblies-with-multiplex-MinION-sequencing) contains links to Figshare for both pre- and post-trimmed FASTQs, which are used in this workshop for demonstrating adapter trimming and quality control analysis, respectively.

The Figshare links are:

- Pre-trimmed: [https://figshare.com/articles/dataset/Basecalled_ONT_reads/5170843](https://figshare.com/articles/dataset/Basecalled_ONT_reads/5170843)
- Post-trimmed: [https://figshare.com/articles/dataset/Trimmed_ONT_reads/5170852](https://figshare.com/articles/dataset/Trimmed_ONT_reads/5170852)

The links to the files are in [data/rwick.figshare.tsv](../data/rwick.figshare.tsv):

```
https://ndownloader.figshare.com/files/8811226	barcode01.fastq.gz
https://ndownloader.figshare.com/files/8811229	barcode02.fastq.gz
https://ndownloader.figshare.com/files/8811232	barcode03.fastq.gz
https://ndownloader.figshare.com/files/8811235	barcode04.fastq.gz
https://ndownloader.figshare.com/files/8811238	barcode05.fastq.gz
https://ndownloader.figshare.com/files/8811241	barcode06.fastq.gz
https://ndownloader.figshare.com/files/8811247	barcode07.fastq.gz
https://ndownloader.figshare.com/files/8811250	barcode08.fastq.gz
https://ndownloader.figshare.com/files/8811253	barcode09.fastq.gz
https://ndownloader.figshare.com/files/8811256	barcode10.fastq.gz
https://ndownloader.figshare.com/files/8811259	barcode11.fastq.gz
https://ndownloader.figshare.com/files/8811262	barcode12.fastq.gz
https://ndownloader.figshare.com/files/8811274	barcode01.trimmed.fastq.gz
https://ndownloader.figshare.com/files/8811280	barcode02.trimmed.fastq.gz
https://ndownloader.figshare.com/files/8811283	barcode03.trimmed.fastq.gz
https://ndownloader.figshare.com/files/8811286	barcode04.trimmed.fastq.gz
https://ndownloader.figshare.com/files/8811289	barcode05.trimmed.fastq.gz
https://ndownloader.figshare.com/files/8811292	barcode06.trimmed.fastq.gz
https://ndownloader.figshare.com/files/8811298	barcode07.trimmed.fastq.gz
https://ndownloader.figshare.com/files/8811301	barcode08.trimmed.fastq.gz
https://ndownloader.figshare.com/files/8811304	barcode09.trimmed.fastq.gz
https://ndownloader.figshare.com/files/8811307	barcode10.trimmed.fastq.gz
https://ndownloader.figshare.com/files/8811310	barcode11.trimmed.fastq.gz
https://ndownloader.figshare.com/files/8811313	barcode12.trimmed.fastq.gz
```

To download the files:

```bash
while read URL NAME; do
    curl -L -o ${NAME} ${URL}
done < rwick.figshare.tsv
```

#### Human long read data

To demonstrate contamination of bacterial reads with human sequences, a small test dataset of human long read sequencing data was downloaded from the [`nanoseq` branch of the `nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/nanoseq/fastq/demultiplexed) repository. The dataset is a small FASTQ containing reads from the NA12878 Genome In A Bottle sample.

```bash
wget https://github.com/nf-core/test-datasets/raw/refs/heads/nanoseq/fastq/demultiplexed/NA12878_DNA.fastq.gz
```

The human data is considerably longer in read length than the bacterial reads, so they were trimmed to 1/5 of their length:

```bash
awk 'NR % 4 == 2 || NR % 4 == 0 { l = int(length($0) / 5); print(substr($0, 1, l)); next } { print $0 }' <(zcat NA1
2878_DNA.fastq.gz) > NA12878_DNA.trimmed.fastq
```

##### Alternate source

The above NA12878 data failed to be classififed by Kraken2 and the Kalamari database. An alternate dataset was downloaded to obtain reads that would be successfully classified. This was obtained from the [Oxford Nanopore Human Reference Datasets GitHub repository](https://github.com/nanopore-wgs-consortium/NA12878/tree/master).

```bash
# Download data
wget http://s3.amazonaws.com/nanopore-human-wgs/rel6/FASTQTars/FAB39088-288418386_Multi.tar

# Extract files
tar -xvf FAB39088-288418386_Multi.tar

# Take a subset of the FASTQs
for f in Notts/FAB39088-288418386_Multi/fastq/fastq_runid_067595119a22414e19bcc3d686d28710d59b35b9_*.fastq.gz; do
    zcat $f
done > NA12878.subset.fastq

# Run Kraken2 to find all reads that successfully classify as human
kraken2 \
    NA12878.subset.fastq \
    --db /cvmfs/data.galaxyproject.org/managed/kraken2_databases/kalamari \
    --report NA12878.subset.k2report \
    --report-minimizer-data \
    --minimum-hit-groups 3 \
    --threads 2 \
    --output NA12878.subset.k2_out.txt
awk '$3 == 9606 { print $2 }' NA12878.subset.k2_out.txt > NA12878.subset.human_reads.txt

# Extract just the classified reads from the FASTQ
awk 'BEGIN { p = 0 } NR == FNR { r[$0] = 1; next } NR > FNR && FNR % 4 == 1 { id = gensub("^@(\\S+)\\s.*$", "\\1", "g", $0); if ( id in r ) { p = 1 } else { p = 0 } } NR > FNR && p == 1 { print $0 }' NA12878.subset.human_reads.txt NA12878.subset.fastq > NA12878.subset.human.fastq

# Run Kraken2 again to test
kraken2 \
    NA12878.subset.human.fastq \
    --db /cvmfs/data.galaxyproject.org/managed/kraken2_databases/kalamari \
    --report NA12878.subset.human.k2report \
    --report-minimizer-data \
    --minimum-hit-groups 3 \
    --threads 2 \
    --output NA12878.subset.human.k2_out.txt
```

### Creating subsampled FASTQ files

The FASTQs were subsampled to extract 10% of their reads for use within the workshop:

```bash
for f in *.fastq; do
    seqkit sample -p 0.1 ${f} -o $(basename ${f} .fastq).subset_0.1.fastq
done
```

**Note** that this command relies on the `seqkit` tool. We used the Docker image `quay.io/biocontainers/seqkit:2.13.0--he881be0_0` to run the above command.

### Creating pre-filtered FASTQ files

#### Prussing et al., 2020

The full FASTQ files from the Prussing et al. study were pre-filtered for use in the assembly lesson, since filtering the sub-sampled FASTQs would result in too few reads for full assembly. The FASTQs were filtered similar to in the original study, using `filtlong` with a minimum read length of 1kb (the same as the study) and a target total number of bases per sample of 250 million (half the value used in the study).

```bash
for f in *.fastq; do
    filtlong --min_length 1kb --target_bases 250mb $f > $(basename $f .fastq).filtered.fastq
done
```

#### Roberts et al., 2023

The Roberts et al. FASTQs were also pre-filtered in a similar manner to their paper. The authors filtered to a minimum read length of 1kb and a minimum quality score of Q7. For the workshop, the FASTQs were filtered to a minimum read length of 1kb and a minimum quality score of Q13 (~95% base accuracy):

```bash
for f in *.fastq.gz; do
    filtlong --min_length 1kb --min_mean_q 95 $f > $(basename $f .fastq.gz).filtered.fastq
done
```

### Creating contaminated FASTQ files

#### Mixing bacterial species

The following is an example of creating a simulated contaminated sample containing two separate species:

```bash
# Make a copy of SRR11909877
cp SRR11909877.filtered.fastq SRR11909877.filtered.contam.bac.fastq
# Get num. reads in SRR11909877
wc -l SRR11909877.filtered.contam.bac.fastq  # 72412
# Approx. 18k reads.

# Take 1000 reads from SRR11909879 (~5%)
# and append to SRR11909877.filtered.contam.bac.fastq
head -n 4000 SRR11909879.filtered.fastq >> SRR11909877.filtered.contam.bac.fastq
```

#### Contamination with human reads

The following is an example of creating a simulated contaminated sample containing bacterial and human reads:

```bash
# Make a copy of SRR11909877
cp SRR11909877.filtered.fastq SRR11909877.filtered.contam.human.fastq

# Append the NA12878 trimmed reads
cat NA12878.subset.human.fastq >> SRR11909877.filtered.contam.human.fastq
```

## Software

### Creating modules on training VMs

The training VMs use BioShell, which makes a wide variety of container images available via the CVMFS shared file system. BioShell provides the `shelley-bio` tool which allows you to build modules that can expose these container images as regular commands. We used the following general steps to use `shelley-bio` and create modules for the various tools used in the workshop:

1. Find the tools you want to create a module for:

```bash
shelley-bio find multiqc
```

2. Determine the version of the tool that you want and install:

```bash
shelley-bio build multiqc/1.35--pyhdfd78af_1
```

**Note** that the tool needs to be provided as `<TOOL_NAME>/<VERSION>`, separated by a forward slash. The `<VERSION>` component can be either the full version tag (e.g. `1.35--pyhdfd78af_1`) or the short-hand version number (e.g. `1.35`). Note that the short-hand version may match with mutliple builds, in which case you will be asked to interactively select the build that you want.

3. Check that the new module is available:

```bash
module avail
```

The new module should appear under the `/apps/Modules/modulefiles` section:

```console
----------------------- /apps/Modules/modulefiles ------------------------
   R/4.3.3
   ansible/2.16.3
   jupyter/2026.04
   multiqc/1.35--pyhdfd78af_1
   nextflow/25.10.4
   nf-core/3.5.2
   rstudio/2023.12.1
   snakemake/7.32.4

------------------------ /opt/Modules/modulefiles ------------------------
   shpc (L)    singularity (L)

  Where:
   L:  Module is loaded

If the avail list is too long consider trying:

"module --default avail" or "ml -d av" to just list the default modules.
"module overview" or "ml ov" to display the number of modules for each
name.

Use "module spider" to find all possible modules and extensions.
Use "module keyword key1 key2 ..." to search for all possible modules
matching any of the "keys".
```

**Note** that the `module use` command will need to be re-run for every session. Place these in `~/.bashrc` to make sure they always run on login.

4. Load and test the module:

```bash
module load multiqc/1.35--pyhdfd78af_1

multiqc --help
```

Output

```console
/// MultiQC 🔍 v1.35

Usage: multiqc [OPTIONS] [ANALYSIS DIRECTORY]

...
```
