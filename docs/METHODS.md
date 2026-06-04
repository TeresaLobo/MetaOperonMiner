# Methods

## Dataset selection

Public paired-end Illumina shotgun metagenomes from aquatic environments were selected from the NCBI Sequence Read Archive. The analysis focused on wastewater-associated and polluted river-associated microbiomes.

## Read preprocessing

Raw reads were quality-filtered using fastp. Reads with low-quality regions, excessive ambiguous bases, or insufficient length were removed.

## Subsampling

Filtered reads were randomly subsampled using seqtk with a fixed random seed to ensure reproducibility.

## Assembly

Subsampled reads were assembled using MEGAHIT. Contigs shorter than 1,000 bp were excluded from downstream genomic-context analyses.

## Gene prediction

Open reading frames were predicted using Prodigal in metagenomic mode.

## HMM-based detection

Predicted proteins were screened against Pfam-A HMM profiles using HMMER. Domains associated with metal resistance systems were retained for downstream analysis.

## Operon-like module reconstruction

Candidate modules were reconstructed by integrating HMMER hits with GFF-derived gene coordinates. Genes were grouped when they occurred on the same contig, shared strand orientation, were located within 10 kb, and included at least two metal-associated domains.

## Mobilome screening

Predicted proteins were screened for MGE-associated Pfam domains. MGE signatures within 10 kb of reconstructed modules were evaluated.

## Visualization

Publication-ready figures were generated using Python, Matplotlib, and NetworkX.
