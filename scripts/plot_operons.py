import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results/operons/metal_operon_candidates.tsv", sep="\t")

counts = df['group'].value_counts()

plt.figure()
counts.plot(kind='bar')

plt.ylabel("Number of operon-like modules")
plt.xlabel("Environment")
plt.title("Distribution of metal resistance operons")

plt.tight_layout()
plt.savefig("results/operons/operon_barplot.png", dpi=300)
