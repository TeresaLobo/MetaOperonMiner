# Methods

## Dataset selection

Public shotgun metagenomic datasets were selected from the NCBI Sequence Read Archive. Datasets were included if they met the following criteria: Illumina paired-end sequencing, shotgun metagenomic strategy, environmental origin, and sufficient sequencing depth for assembly-based genomic context analysis.

## Read preprocessing

Raw reads were converted to FASTQ using SRA Toolkit. Quality filtering and adapter trimming were performed using fastp with a minimum Phred score of Q20, a maximum low-quality base percentage of 30%, and a minimum post-trimming length of 50 bp.

## Read subsampling

Filtered reads were randomly subsampled using seqtk with a fixed seed. Five million paired reads per sample were retained to ensure reproducibility and reduce computational memory usage.

## Metagenomic assembly

Subsampled reads were assembled using MEGAHIT. Only contigs equal to or longer than 1000 bp were retained for downstream analyses.

## Gene prediction

Coding sequences were predicted using Prodigal in metagenomic mode. Protein sequences, nucleotide sequences, and GFF files containing genomic coordinates were generated for each sample.

## HMM-based metal resistance detection

Predicted proteins were screened against Pfam-A HMM profiles using HMMER. The analysis focused on Pfam domains associated with metal resistance systems, including Cus/Sil-like efflux, Pco/Cop multicopper oxidases, ArsC arsenate reductases, Mer-like reductases, and Czc-like cation efflux systems.

## Candidate operon-like module reconstruction

Candidate operon-like modules were reconstructed by combining HMMER hits with GFF-derived genomic coordinates. Genes were grouped as candidate modules when they occurred on the same contig, shared the same strand orientation, were separated by no more than 10 kb, and contained at least two metal-resistance-associated domains.

## Statistical analysis

Descriptive statistics were calculated for gene counts, operon counts, operon classes, and operon lengths. The Mann–Whitney U test was used as an exploratory comparison of operon length distributions between environmental groups.

## Visualization

Figures were generated using Python libraries including Pandas, Matplotlib, and NetworkX.
