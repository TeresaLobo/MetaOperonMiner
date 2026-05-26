import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("results/operons/metal_operon_candidates.tsv", sep="\t")

df_expanded = df.assign(system=df['systems'].str.split(',')).explode('system')

matrix = pd.crosstab(df_expanded['sample'], df_expanded['system'])

plt.figure()
sns.heatmap(matrix, annot=True)

plt.title("Metal resistance operon systems")
plt.xlabel("System")
plt.ylabel("Sample")

plt.tight_layout()
plt.savefig("results/operons/operon_heatmap.png", dpi=300)
