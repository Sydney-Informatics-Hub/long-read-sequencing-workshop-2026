# Workshop development notes

## Downloading FASTQ files

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

## Creating subsampled FASTQ files

The FASTQs were subsampled to extract 10% of their reads for use within the workshop:

```bash
for f in *.fastq; do
    seqkit sample -p 0.1 ${f} -o $(basename ${f} .fastq).subset_0.1.fastq
done
```

**Note** that this command relies on the `seqkit` tool. We used the Docker image `quay.io/biocontainers/seqkit:2.13.0--he881be0_0` to run the above command.

## Creating pre-filtered FASTQ files

The full FASTQ files were also pre-filtered for use in the assembly lesson, since filtering the sub-sampled FASTQs would result in too few reads for full assembly. The FASTQs were filtered similar to in the original study, using `filtlong` with a minimum read length of 1kb (the same as the study) and a target total number of bases per sample of 250 million (half the value used in the study).

```bash
for f in *.fastq; do
    filtlong --min_length 1kb --target_bases 250mb $f > $(basename $f .fastq).filtered.fastq
done
```

## Creating modules on training VMs

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
