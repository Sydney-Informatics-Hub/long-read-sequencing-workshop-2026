# Workshop data

The data used in this workshop comes from [Prussing et al. *Nanopore MinION Sequencing Reveals Possible Transfer of blaKPC-2 Plasmid Across Bacterial Species in Two Healthcare Facilities.* **Front Microbiol.** 2020 Aug 19;11:2007.](https://doi.org/10.3389/fmicb.2020.02007)

In this study, five bacterial isolates from two healthcare facilities were sequenced using both ONT long read sequencing and Illumina short read sequencing. The authors found evidence for horizontal transfer of a plasmid between three species of carbapenemase-producing *Enterobacteriaceae* bacteria, including *Klebsiella pneumoniae*, *Citrobacter freundii*, and *Escherichia coli*. Four of the five isolates contained a plasmid that was entirely, or almost entirely identical to a known plasmid (pKPC_UVA01) that contains the AMR gene *bla_KPC-2*.

All of the data for this study was made publicly available on the [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject) archive (accession ID: [PRJNA636827](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA636827)). The FASTQs used in this study were created by subsetting the original ONT long read FASTQs downloaded from the Sequence Read Archive (SRA) entries associated with this study. The table below lists important information about each sample:

| Sample | Species | Isolate/patient number from paper | Facility | AMR gene *bla_KPC-2* present? | Plasmid pKPC_UVA01 present? |
| ------ | ------- | --------------------------------- | -------- | ----------------------------- | --------------------------- |
| SRR11909885 | C. freundii | 1 | A | Yes | Yes |
| SRR11909883 | K. pneumoniae | 2 | A, B | Yes | Yes |
| SRR11909881 | K. pneumoniae | 3 | A, B | Yes | Yes |
| SRR11909879 | C. freundii | 4 | A | Yes | No |
| SRR11909877 | E. coli | 5 | B | Yes | Yes |

## Subsetting data

The raw FASTQs from the original study are large (~2GB and ~200k reads). For the workshop, we will try to use a subset of reads to achieve a similar results.

Initial tests with SRR11909877 show that subsetting the FASTQ file to ~10% of the original reads still allows for assembly and circularisation of the main chromosome (~4Mb) and the AMR-containing plasmid (~43kb).