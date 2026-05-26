import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact

df = pd.read_csv("results/operons/metal_operon_candidates.tsv", sep="\t")
genes = pd.read_csv("results/operons/metal_genes_detected.tsv", sep="\t")

print("\n=== Operons por grupo ===")
print(df["group"].value_counts())

print("\n=== Genes metálicos por grupo ===")
print(genes["group"].value_counts())

print("\n=== Operons por amostra ===")
operons_sample = df.groupby(["group", "sample"]).size().reset_index(name="n_operons")
print(operons_sample)

print("\n=== Comprimento dos operons por grupo ===")
print(df.groupby("group")["length_bp"].describe())

print("\n=== Classes de operons ===")
print(df["operon_class"].value_counts())

print("\n=== Sistemas detectados em operons ===")
systems = df.assign(system=df["systems"].str.split(",")).explode("system")
print(pd.crosstab(systems["group"], systems["system"]))

# Teste exploratório: comprimento dos operons entre grupos
groups = df["group"].unique()

if len(groups) == 2:
    g1 = df[df["group"] == groups[0]]["length_bp"]
    g2 = df[df["group"] == groups[1]]["length_bp"]

    if len(g1) >= 2 and len(g2) >= 2:
        stat, p = mannwhitneyu(g1, g2, alternative="two-sided")
        print("\n=== Mann-Whitney U: comprimento dos operons ===")
        print(f"{groups[0]} vs {groups[1]}: U={stat}, p={p}")
    else:
        print("\nPoucos operons por grupo para teste Mann-Whitney confiável.")

print("\nObservação: devido ao baixo número de operons candidatos, os testes devem ser interpretados como exploratórios.")
