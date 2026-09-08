#!/bin/bash

# Install all modules required for workshop
shelley build multiqc/1.35--pyhdfd78af_1
shelley build nanoplot/1.47.0--pyhdfd78af_0
shelley build fastqc/0.12.1--hdfd78af_0
shelley build fastplong/0.4.1--h224cc79_0
shelley build filtlong/0.3.1--h077b44d_0
shelley build kraken2/2.17.1--pl5321h077b44d_0
shelley build flye/2.9.6--py313h7fbb527_1
shelley build raven-assembler/1.8.3--h5ca1c30_3
shelley build canu/2.3--h636b4d1_3
shelley build seqkit/2.13.0--he881be0_0
shelley build seqtk/1.4--h577a1d6_3
shelley build plassembler/1.8.2--pyhdfd78af_0
shelley build quast/5.3.0--py313pl5321h5ca1c30_2
shelley build busco/6.1.0--pyhdfd78af_1
shelley build bandage/0.9.0--h9948957_0
shelley build medaka/2.2.2--py312h3050eb1_0
shelley build medaka/1.3.3--py38h130def0_0
shelley build bakta/1.12.0--pyhdfd78af_0
shelley build ncbi-amrfinderplus/4.2.7--hf69ffd2_0
shelley build autocycler/0.5.2--h3ab6199_0

# Create data directory
mkdir -p data/ref
cd data/ref

# Link CVMFS databases
ln -s /cvmfs/data.galaxyproject.org/managed/kraken2_databases/kalamari
ln -s /cvmfs/data.galaxyproject.org/byhand/bakta_database/10522951 bakta

# Download non-CVMFS databases
shopt -s expand_aliases
module load plassembler
module load ncbi-amrfinderplus
module load busco

## Plassembler
plassembler-exec plassembler download -d plasmid_db_plassembler

## BUSCO
mkdir busco
cd busco
wget https://busco-data.ezlab.org/v6/data/lineages/bacteria_odb12.2.2026-05-22.tar.gz
tar -xzf bacteria_odb12.2.2026-05-22.tar.gz
rm bacteria_odb12.2.2026-05-22.tar.gz
cd ..

## AMRFinderPlus
amrfinder_update -d amrfinderplus_db

# Download FASTQs
cd ..
mkdir fastqs
cd fastqs

wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR828/002/ERR8282742/ERR8282742.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR828/001/ERR8282751/ERR8282751.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR828/002/ERR8282752/ERR8282752.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR828/003/ERR8282753/ERR8282753.fastq.gz

# Subset large FASTQs
module load seqkit
