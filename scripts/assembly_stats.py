from pathlib import Path
import pandas as pd

def fasta_lengths(fasta):
    lengths = []
    seq = []
    with open(fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    lengths.append(sum(len(x) for x in seq))
                seq = []
            else:
                seq.append(line)
        if seq:
            lengths.append(sum(len(x) for x in seq))
    return lengths

def n50(lengths):
    total = sum(lengths)
    csum = 0
    for l in sorted(lengths, reverse=True):
        csum += l
        if csum >= total / 2:
            return l
    return 0

rows = []

for fasta in Path("results/assembly_subsample").glob("*/*_megahit/final.contigs.fa"):
    group = fasta.parts[2]
    sample = fasta.parts[3].replace("_megahit", "")
    lengths = fasta_lengths(fasta)
    rows.append({
        "group": group,
        "sample": sample,
        "n_contigs": len(lengths),
        "total_bp": sum(lengths),
        "largest_contig": max(lengths) if lengths else 0,
        "mean_length": round(sum(lengths)/len(lengths), 2) if lengths else 0,
        "N50": n50(lengths)
    })

df = pd.DataFrame(rows)
df.to_csv("manuscript/tables/Table2_assembly_statistics.tsv", sep="\t", index=False)
df.to_csv("manuscript/supplementary/tables/Table_S3_assembly_statistics.tsv", sep="\t", index=False)

print(df)
