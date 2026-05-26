# MetaOperonMiner: metagenome-resolved reconstruction of metal resistance operon-like modules

This repository contains the reproducible workflow used to identify and reconstruct candidate metal resistance operon-like modules from public environmental metagenomes.

## Study overview

This study investigates the genomic organization of metal-resistance-associated systems in environmental metagenomes. Instead of quantifying isolated resistance genes only, the workflow reconstructs candidate operon-like modules by integrating HMM-based domain detection with contig-level genomic coordinates.

## Environmental groups

Two environmental groups were analyzed:

- Wastewater metagenomes
- Polluted river metagenomes

## Main workflow

1. Download public SRA datasets
2. Convert SRA files to FASTQ
3. Perform quality control with fastp
4. Subsample reads using seqtk
5. Assemble metagenomes with MEGAHIT
6. Predict genes with Prodigal
7. Screen predicted proteins using HMMER and Pfam
8. Filter metal-resistance-associated domains
9. Reconstruct candidate operon-like modules
10. Generate statistical summaries and figures

## Metal resistance systems screened

| Pfam | System | Interpretation |
|---|---|---|
| PF02614 | Cus/Sil-like | Copper/silver efflux |
| PF00394 | Pco/Cop | Multicopper oxidase |
| PF03960 | ArsC | Arsenic resistance |
| PF13450 | Mer-like | Mercury resistance |
| PF01545 | Czc-like | Cation efflux |

## Main outputs

- `metal_genes_detected.tsv`
- `metal_operon_candidates.tsv`
- assembly statistics
- operon distribution figures
- heatmap of resistance systems
- operon network figure

## Reproducibility

The workflow was executed in Linux Ubuntu using Bash and Python scripts. Raw sequencing files are not deposited in this repository because they are publicly available in the NCBI SRA.

## Citation

If you use this workflow, please cite this repository and the associated manuscript.
