#!/bin/bash

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
