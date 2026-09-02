# Workshop scripts

This directory contains all of the final scripts we will build during the workshop. The scripts make up a complete workflow from FASTQ quality control, through genome assembly, and to antimicrobial resistance gene detection and plasmid analysis. The scripts are:

0. `00.multiqc.sh`: MultiQC script for compiling all outputs into one report.
1. `01.nanoplot_summary.sh`: Initial sequencing run QC using a sequencing summary file as input.
2. Initial FASTQ quality control:
    1. `02.a.nanoplot_fastq.sh`: Initial FASTQ quality control with NanoPlot.
    2. `02.b.fastqc.sh`: Initial FASTQ quality control with FASTQC.
3. `03.trim_filter.sh`: Read trimming and filtering to remove adapters and low-quality bases/reads using fastplong.
4. `04.a.contamination.sh`: Contamination assessment with Kraken2.
    1. `04.b.remove_contamination.sh`: Host read removal by aligning to human genome and removing aligned reads.
5. Genome assembly:
    1. `05.a.assembly.flye.sh`: Main genome assembly with Flye.
    2. `05.b.assembly.plassembler.sh`: Plasmid assembly with Plassembler.
    3. `05.c.assembly.raven.sh`: Alternate genome assembly with Raven.
6. `06.consensus.sh`: Consensus assembly creation with Autocycler.
7. `07.assembly_qc.sh`: Assembly quality control and visualisation with QUAST, BUSCO, MultiQC, and Bandage.
8. `08.polish.sh`: Assembly polishing with Medaka.
9. `09.annotate.sh`: Perform genome annotation with Bakta.
10. `10.amr_gene_detection.sh`: Detect antimicrobial resistance (AMR) and virulence genes with AMRFinderPlus and ABRicate.
11. `11.plasmid_comparison.sh`: Perform a multiple sequence alignment on the plasmid sequences to detect similarities.
12. `12.strain_id.sh`: Identify bacterial strains with cgMLST.

## Workshop outline

1. Day 1
    1. Introduction
        - Introduction to workshop & learning objectives
        - General overview of long read sequencing and applications to bacterial genome analysis
        - Introduce samples we're using and the context of hospital-aquired infections and antibiotic resistance
        - Format: lecture with slides
    2. Long read sequencing QC
        - Introductory lecture with slides
            - Discuss tools and metrics for QC analysis
        - Exercises:
            - Performing qequencing run and FASTQ quality control
            - Inspecting and interpreting QC reports
    3. Species identification and contaimination analysis
        - Introductory lecture with slides
            - Discuss methods for detecting species of origin of reads, and discuss identifying and removing contaminating reads
        - Exercises:
            - Run Kraken2 to detect species of origin of reads
            - Remove host genes by aligning to human reference and discarding mapped genes
                - **Note:** could be skipped, left to participants to try out for themselves
    4. Bacterial genome assembly
        - Introductory lecture with slides
            - Discuss how assembly works, esp. in context of bacterial genomics
                - E.g. trying to identify circularised genomes
            - Discuss main metrics for assessing assembly quality
                - N50, L50, etc.
        - Exercises:
            - Run Flye to assemble bacterial genomes
            - Run Plassembler to identify potential plasmids that Flye missed
            - Run another assembler (e.g. Raven) to create an alternate assembly for consensus assembly
            - Perform assembly QC and visualisation
            - Inspect and interpret QC reports
2. Day 2
    1. Introduction
        - Recap day 1
        - Overview of Day 2 learning objectives
        - Format: lecture with slides
    2. Consensus assemblies
        - Introductory lecture with slides
            - Discuss the benefits of constructing a consensus assembly
        - Exercises:
            - Run Autocycler to generate a consensus assembly
            - Perform assembly QC and visualisation
            - Inspect and interpret QC reports
        - **Note:** could be skipped, left to participants to try out for themselves
    3. Assembly polishing
        - Introductory lecture with slides
            - Discuss the need for assembly polishing and methods available
        - Exercises:
            - Run medaka to polish the genome assembly
            - Perform assembly QC and visualisation
            - Inspect and interpret QC reports
    4. Genome annotation
        - Introductory lecture with slides
            - Discuss annotation options
            - Discuss the goal for the lesson: annotating AMR genes
        - Exercises:
            - Run AMRFinderPlus on genome assembly
            - Run Abricate on genome assembly
            - Integrate and inspect results of AMR gene annotation
            - Identify location of AMR genes in genome: chromosome vs plasmids
    5. Comparitive genomics
        - Introductory lecture with slides
            - Discuss the goal for the lesson: identifying plasmid sequence similarity between samples, identifying horizontal gene transfer events
            - Recap the origin of each sample and the timeline of infections
        - Exercises:
            - Perform phylogenetic analysis of samples
            - Identify bacterial strains with cgMLST
                - **Note:** could be skipped, left to participants to try out for themselves
            - Perform multiple sequence alignment of plasmids
            - Interpret results: almost identical plasmids in two different species - horizontal gene transfer

## Container images

The following container images are being used for each tool:

| Tool | Image URL |
| ---- | --------- |
| MultiQC | docker://quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1 |
| NanoPlot | docker://quay.io/biocontainers/nanoplot:1.47.0--pyhdfd78af_0 |
| FastQC | docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0 |
| fastplong | docker://quay.io/biocontainers/fastplong:0.4.1--h224cc79_0 |
| filtlong | docker://quay.io/biocontainers/filtlong:0.3.1--h077b44d_0 |
| Kraken2 | docker://quay.io/biocontainers/kraken2:2.17.1--pl5321h077b44d_0 |
| Flye | docker://quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1 |
| Raven |  |
| Canu |  |
| seqkit | docker://quay.io/biocontainers/seqkit:2.13.0--he881be0_0 |
| seqtk | docker://quay.io/biocontainers/seqtk:1.4--h577a1d6_3 |
| Plassembler | docker://quay.io/biocontainers/plassembler:1.8.2--pyhdfd78af_0 |
| QUAST | docker://quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2 |
| BUSCO | docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_1 |
| Bandage | docker://quay.io/biocontainers/bandage:0.9.0--h9948957_0 |
| medaka | docker://quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0 |
| bakta | docker://quay.io/biocontainers/bakta:1.12.0--pyhdfd78af_0 |
| ncbi-amrfinderplus | docker://quay.io/biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0 |
| Autocycler | docker://quay.io/biocontainers/autocycler:0.5.2--h3ab6199_0 |
| MUMmer | docker://quay.io/biocontainers/mummer:3.23--pl526_6 |

### Install script

See the installation script [`install.sh`](install.sh) for the commands used to install the above containers via BioShell. Also see the [development documentation](../dev/development.md) for more information on installing BioContainer modules on the VMs.

## Databases

The following tools require databases:

| Tool | Database Location | Notes |
| ---- | ----------------- | ----- |
| Kraken2 | //cvmfs/data.galaxyproject.org/managed/kraken2_databases/kalamari | Kalamari database is curated and only 2GB in size, which is suitable for the workshop, while the standard Kraken2 database is ~60GB in size and can't fit into memory on the VMs. |
| Plassembler | N/A | Download with `plassembler download -d <db directory>`. Size is 437MB, can be pre-loaded on the VMs. |

## Run order

### Ryan Wick samples (`barcode*.fastq.gz`)

1. FASTQC

```bash
for f in data/fastqs/barcode*.fastq.gz; do
    ./long-read-sequencing-workshop-2026/scripts/02.fastqc.sh $f
done
```

2. fastplong with adapter trimming, Q13 filter, and poly-X filter

```bash
for f in data/fastqs/barcode*.fastq.gz; do
    ./long-read-sequencing-workshop-2026/scripts/03.c.fastplong.filter.sh $f
done
```

3. Re-run FASTQC on trimmed/filtered FASTQs

```bash
for f in trim_filter/*.fastq; do
    ./long-read-sequencing-workshop-2026/scripts/02.fastqc.sh $f
done
```

4. Species ID

```bash
for f in trim_filter/*.fastq; do
    ./long-read-sequencing-workshop-2026/scripts/04.a.species_id.sh $f
done
```

5. Flye assembly

```bash
for f in trim_filter/*.fastq; do
    ID=$(basename $f .trimmed.filtered.fastq)
    ./long-read-sequencing-workshop-2026/scripts/05.a.assembly.flye.sh $f $ID
done
```

### Outbreak samples (`{SRR,ERR}*.fastq.gz`)

1. FASTQC

```bash
for f in data/fastqs/{SRR,ERR}*.fastq.gz; do
    ./long-read-sequencing-workshop-2026/scripts/02.fastqc.sh $f
done
```

2. fastplong to remove adapters and poly-X

```bash
for f in data/fastqs/barcode*.fastq.gz; do
    ./long-read-sequencing-workshop-2026/scripts/03.b.fastplong.adapters.sh $f
done
```

3. filtlong with 1kb minimum length and 50% best reads filters

```bash
for f in trim/{SRR,ERR}*.fastq.gz; do
    ./long-read-sequencing-workshop-2026/scripts/03.a.filtlong.sh $f
done
```

3. Re-run FASTQC on trimmed/filtered FASTQs

```bash
for f in {trim,filter}/*.fastq.gz; do
    ./long-read-sequencing-workshop-2026/scripts/02.fastqc.sh $f
done
```

4. Species ID

```bash
for f in filter/*.fastq.gz; do
    ./long-read-sequencing-workshop-2026/scripts/04.a.species_id.sh $f
done
```

5. Flye assembly

```bash
for f in filter/*.fastq.gz; do
    ID=$(basename $f .trimmed.filtered.fastq.gz)
    ./long-read-sequencing-workshop-2026/scripts/05.a.assembly.flye.sh $f $ID
done
```

6. Polishing

```bash
# TODO
```

7. AMRFinder

```bash

```