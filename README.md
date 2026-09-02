# Australian BioCommons Long Read Sequencing Workshop 2026

These are the training materials for the 2026 Australian BioCommons workshop for long read sequencing data analysis (title TBD 🚧). This workshop was developed and presented by [Sydney Informatics Hub](https://www.sydney.edu.au/research/facilities/sydney-informatics-hub.html) (SIH) at the University of Sydney, [Queensland Cyber Infrastructure Foundation](https://www.qcif.edu.au/) (QCIF), and [Australian BioCommons](https://www.biocommons.org.au/).

This workshop is designed to take you through the process of analysing Oxford Nanoport Technologies (ONT) long read sequencing data within the context of bacterial genome assembly and analysis. It is intended for researchers and bioinformaticians working within the microbiology field who are interested in learning more about long read sequencing analysis and how it can be applied to their research.

## Pre-requisites

🚧

## Learning objectives

🚧

## For developers

These materials are written in Markdown and are intended to be rendered using [Material for MkDocs](https://github.com/squidfunk/mkdocs-material). All of the pre-requisite software for developing and locally rendering the materials are defined within a [uv](https://docs.astral.sh/uv/) project. To render the materials locally, follow these steps:

1. Install uv [as per the documentation](https://docs.astral.sh/uv/getting-started/installation/):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

2. Syncronise the project's virtual environment

```bash
uv sync
```

3. Run `mkdocs serve`:

```bash
uv run mkdocs serve
```

4. Open the materials in a web browser at: `http://localhost:8000`
