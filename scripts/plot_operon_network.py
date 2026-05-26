import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

df = pd.read_csv("results/operons/metal_operon_candidates.tsv", sep="\t")

G = nx.Graph()

for _, row in df.iterrows():
    operon = row["operon_id"]
    sample = row["sample"]
    group = row["group"]
    systems = str(row["systems"]).split(",")

    G.add_node(operon, type="operon", group=group)
    G.add_node(sample, type="sample", group=group)
    G.add_edge(sample, operon)

    for system in systems:
        system = system.strip()
        G.add_node(system, type="system", group="system")
        G.add_edge(operon, system)

plt.figure(figsize=(12, 8))

pos = nx.spring_layout(G, seed=42, k=0.8)

node_sizes = []
for node in G.nodes():
    if G.nodes[node]["type"] == "sample":
        node_sizes.append(900)
    elif G.nodes[node]["type"] == "operon":
        node_sizes.append(650)
    else:
        node_sizes.append(500)

nx.draw_networkx_nodes(G, pos, node_size=node_sizes)
nx.draw_networkx_edges(G, pos, width=1.2, alpha=0.7)
nx.draw_networkx_labels(G, pos, font_size=7)

plt.title("Network of candidate metal resistance operon-like modules")
plt.axis("off")
plt.tight_layout()
plt.savefig("results/operons/operon_network.png", dpi=300)
plt.savefig("results/operons/operon_network.pdf")
