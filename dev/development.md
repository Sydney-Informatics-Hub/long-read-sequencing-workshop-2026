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

The training VMs use BioShell, which makes a wide variety of container images available via the CVMFS shared file system. BioShell also provides `shpc` which allows you to build modules that can expose these container images as regular commands. We used the following general steps to use `shpc` and create modules for the various tools used in the workshop:

1. Enable shpc:

```bash
module load shpc
```

2. Find the tools you want to create a module for:

```bash
shpc show -f multiqc
```

Output:

```console
quay.io/biocontainers/multiqc-bcbio
quay.io/biocontainers/multiqc-xenium-extra
quay.io/biocontainers/multiqc
quay.io/biocontainers/multiqc_sav
quay.io/biocontainers/pmultiqc
quay.io/biocontainers/zavolan-multiqc-plugins
```

3. Find the version of the tool you want:

```bash
shpc show quay.io/biocontainers/multiqc
```

```console
url: https://biocontainers.pro/tools/multiqc
maintainer: '@vsoch'
description: shpc-registry automated BioContainers addition for multiqc
latest:
  1.35--pyhdfd78af_1: 
    sha256:b65e3fe879df27b92334dda0fd987a6e21bdee09a2848551d4f287099a93b7ac
tags:
  1.9--py_1: 
    sha256:67cc651cb350b1ee2fc0929bd6bcd5189ec8c17f09566a3cd54cde7479e48a09
  1.10.1--pyhdfd78af_1: 
    sha256:c64ea8fcaf49dfc4b0594bc7349e6d1a662eb4484f5aac3252f4eea86cad164c

...

  1.34--pyhdfd78af_0: 
    sha256:c6cbb73af77cc2eb59a926e41967dd2a7bdbb12cefe7e325c130ecde7e8730af
  1.35--pyhdfd78af_1: 
    sha256:b65e3fe879df27b92334dda0fd987a6e21bdee09a2848551d4f287099a93b7ac
docker: quay.io/biocontainers/multiqc
aliases:
  multiqc: /usr/local/bin/multiqc
```

4. Determine the container image in the CVMFS file system. All images will be located at `/cvmfs/singularity.galaxyproject.org/all/<TOOL_NAME>:<VERSION>`. For example: `/cvmfs/singularity.galaxyproject.org/all/multiqc:1.35--pyhdfd78af_1`.

5. Install the module with `shpc install <REGISTRY>/<TOOL>:<TAG> <CVMFS_PATH> --keep-path`.

- `<REGISTRY>/<TOOL>` should be the full image name shown in the `docker` field in the `shpc show` output. In the above example, this is `quay.io/biocontainers/multiqc`.
    - `<REGISTRY>` will be everything before the last part of this name, e.g. `quay.io/biocontainers`. This is the Docker registry that hosts the image.
    - `<TOOL>` will the the last part of this name, e.g. `multiqc`.
- `<TAG>` is the image tag that you want to use. The latest version of `multiqc` in the example above has the tag `1.35--pyhdfd78af_1`.
- `<CVMFS_PATH>` is the path you determined in step 4. In this example, this is `/cvmfs/singularity.galaxyproject.org/all/multiqc:1.35--pyhdfd78af_1`.
- `--keep-path` ensures that `shpc` **does not** download the container image, and instead uses the image at the `<CVMFS_PATH>`.

```bash
shpc install \
    quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1 \
    /cvmfs/singularity.galaxyproject.org/all/multiqc:1.35--pyhdfd78af_1 \
    --keep-path
```

Output:

```console
Module quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1 was created.
```

6. Make the new module available for use:

```bash
module use ~/shpc/modules/<REGISTRY>
```

**Note** that `<REGISTRY>` here is the same as in step 5 when running `shpc install`.

Check that the new module is available:

```bash
module avail
```

Output:

```console
---------------------------------------------------------------------- /home/user/shpc/modules/quay.io/biocontainers -----------------------------------------------------------------------
   multiqc/1.35--pyhdfd78af_1/module

--------------------------------------------------------------------------------- /apps/Modules/modulefiles ---------------------------------------------------------------------------------
   R/4.3.3    ansible/2.16.3    jupyter/2026.04    nextflow/25.10.4    nf-core/3.5.2    rstudio/2023.12.1    snakemake/7.32.4

--------------------------------------------------------------------------------- /opt/Modules/modulefiles ----------------------------------------------------------------------------------
   shpc (L)    singularity (L)

  Where:
   L:  Module is loaded

Module defaults are chosen based on Find First Rules due to Name/Version/Version modules found in the module tree.
See https://lmod.readthedocs.io/en/latest/060_locating.html for details.

If the avail list is too long consider trying:

"module --default avail" or "ml -d av" to just list the default modules.
"module overview" or "ml ov" to display the number of modules for each name.

Use "module spider" to find all possible modules and extensions.
Use "module keyword key1 key2 ..." to search for all possible modules matching any of the "keys".
```

**Note** that the `module use` command will need to be re-run for every session. Place these in `~/.bashrc` to make sure they always run on login.

7. Load and test the module:

```bash
module load multiqc/1.35--pyhdfd78af_1/module

multiqc --help
```

Output

```console
/// MultiQC 🔍 v1.35

Usage: multiqc [OPTIONS] [ANALYSIS DIRECTORY]

...
```
