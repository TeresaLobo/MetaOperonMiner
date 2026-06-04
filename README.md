# MetaOperonMiner

**Metagenome-resolved reconstruction of metal resistance operon-like modules in aquatic environments**

[![GitHub release](https://img.shields.io/github/v/release/TeresaLobo/MetaOperonMiner)](https://github.com/TeresaLobo/MetaOperonMiner/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20544294)](https://doi.org/10.5281/zenodo.20544294)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Reproducible workflow](https://img.shields.io/badge/workflow-reproducible-brightgreen.svg)](README.md)

## Overview

MetaOperonMiner is a reproducible bioinformatics workflow designed to identify and reconstruct candidate metal resistance operon-like modules from public environmental metagenomes.

The workflow integrates:

- metagenomic read preprocessing;
- quality filtering;
- read subsampling;
- de novo metagenomic assembly;
- gene prediction;
- HMM-based domain screening;
- reconstruction of operon-like modules;
- mobilome-associated screening;
- statistical analysis;
- publication-ready figure generation.

This repository supports the manuscript:

**Metagenomic reconstruction of mobile-associated metal resistance operon-like modules in aquatic environments**

## Scientific rationale

Environmental microbiomes exposed to anthropogenic contamination are shaped by multiple selective pressures, including heavy metals, biocides, antibiotics, and industrial pollutants. Most environmental resistome studies focus on the abundance of individual resistance genes. MetaOperonMiner instead focuses on genomic organization by reconstructing candidate operon-like modules from assembled metagenomic contigs.

## Workflow

```mermaid
flowchart TD
A[Public SRA metagenomes] --> B[FASTQ conversion]
B --> C[Quality filtering with fastp]
C --> D[Read subsampling with seqtk]
D --> E[Metagenomic assembly with MEGAHIT]
E --> F[Gene prediction with Prodigal]
F --> G[HMM screening with Pfam/HMMER]
G --> H[Metal-resistance domain filtering]
H --> I[Operon-like module reconstruction]
I --> J[Mobilome screening]
J --> K[Statistics and visualization]
