# run_sample_pipeline.sh

`scripts/run_sample_pipeline.sh` runs the workshop pipeline for a single FASTQ from preprocessing through assembly, annotation, and report generation.

## Usage

```bash
bash scripts/run_sample_pipeline.sh [options] <fastq>
```

Example:

```bash
bash scripts/run_sample_pipeline.sh data/sample.fastq
```

## Stages

The wrapper runs these existing scripts in order unless skipped:

1. `02.fastqc.sh`
2. `03.b.fastplong.sh`
3. `04.a.species_id.sh`
4. `05.a.assembly.flye.sh`
5. `05.b.assembly.plassembler.sh`
6. `06.a.consensus.subsample.sh`
7. `06.b.consensus.autocycler.sh`
8. `08.polish.sh`
9. `07.a.assembly_qc.quast.sh`
10. `07.b.assembly_qc.busco.sh`
11. `07.c.assembly_qc.bandage.sh`
12. `09.a.get_plasmids.sh` when `--chromosome-contig` is provided
13. `09.b.amrfinder.sh`
14. `09.c.bakta.sh`
15. `00.multiqc.sh`

## Options

- `--skip-fastqc`: skip FastQC
- `--skip-trim`: skip fastplong trimming
- `--skip-species-id`: skip Kraken2 species identification
- `--skip-flye`: skip Flye assembly
- `--skip-plassembler`: skip Plassembler assembly
- `--skip-consensus`: skip Autocycler consensus generation
- `--skip-polish`: skip Medaka polishing
- `--skip-assembly-qc`: skip QUAST, BUSCO, and Bandage
- `--skip-amrfinder`: skip AMRFinderPlus
- `--skip-bakta`: skip Bakta annotation
- `--skip-multiqc`: skip final MultiQC summary
- `--chromosome-contig <name>`: extract plasmid contigs by excluding the named chromosome contig
- `-h`, `--help`: show CLI help

## Notes

- The script can be run from any directory; it switches into the repository root before running stage scripts.
- Downstream sample IDs follow the existing stage script behavior and are derived from the current FASTQ basename. After trimming, a sample such as `sample.fastq` becomes `sample.trimmed.fastq`, so later outputs use `sample.trimmed` as the ID.
- Consensus and Medaka output FASTA files are detected automatically from their output directories.
