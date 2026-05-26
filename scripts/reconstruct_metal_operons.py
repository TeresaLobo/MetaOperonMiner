#!/usr/bin/env python3

import os
import glob
import re
import csv
from collections import defaultdict

ANNOT_DIR = "results/annotation_gff"
HMM_DIR = "results/hmmer"
OUT_DIR = "results/operons"

EVAL_CUTOFF = 0.05
SCORE_CUTOFF = 12.0
MAX_GAP = 10000
MIN_GENES_CLUSTER = 2

PFAM_MAP = {
    "PF02614": "sil_cus_RND_efflux",
    "PF00394": "pco_cop_multicopper_oxidase",
    "PF03960": "ars_arsenate_reductase",
    "PF13450": "mer_mercuric_reductase",
    "PF01545": "czc_cation_efflux",
}

os.makedirs(OUT_DIR, exist_ok=True)


def normalize_pfam(acc):
    return acc.split(".")[0]


def infer_contig_from_gene(gene_id):
    parts = gene_id.split("_")
    if len(parts) < 2:
        return gene_id
    return "_".join(parts[:-1])


def parse_gff(gff_file):
    """
    Prodigal protein IDs usually follow:
    contig_1, contig_2, contig_3...
    We reconstruct this by ordering CDS per contig.
    """
    cds_by_contig = defaultdict(list)

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            if cols[2] != "CDS":
                continue

            contig = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]

            cds_by_contig[contig].append({
                "contig": contig,
                "start": start,
                "end": end,
                "strand": strand,
            })

    gene_table = {}

    for contig, genes in cds_by_contig.items():
        genes_sorted = sorted(genes, key=lambda x: x["start"])
        for idx, gene in enumerate(genes_sorted, start=1):
            gene_id = f"{contig}_{idx}"
            gene["gene_id"] = gene_id
            gene["gene_number"] = idx
            gene_table[gene_id] = gene

    return gene_table


def parse_hmmer(tbl_file):
    hits = {}

    with open(tbl_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.strip().split()
            if len(cols) < 6:
                continue

            target = cols[0]
            pfam_name = cols[2]
            pfam_acc = normalize_pfam(cols[3])

            try:
                evalue = float(cols[4])
                score = float(cols[5])
            except ValueError:
                continue

            if pfam_acc not in PFAM_MAP:
                continue

            if evalue <= EVAL_CUTOFF and score >= SCORE_CUTOFF:
                hits[target] = {
                    "gene_id": target,
                    "pfam_name": pfam_name,
                    "pfam_acc": pfam_acc,
                    "system": PFAM_MAP[pfam_acc],
                    "evalue": evalue,
                    "score": score,
                }

    return hits


def build_operons(sample, gene_table, hits):
    annotated_hits = []

    for gene_id, hit in hits.items():
        if gene_id not in gene_table:
            continue

        g = gene_table[gene_id]
        annotated_hits.append({
            **hit,
            **g,
            "sample": sample,
        })

    by_contig = defaultdict(list)
    for h in annotated_hits:
        by_contig[h["contig"]].append(h)

    clusters = []
    cluster_id = 0

    for contig, genes in by_contig.items():
        genes = sorted(genes, key=lambda x: x["start"])

        current = []

        for gene in genes:
            if not current:
                current = [gene]
                continue

            previous = current[-1]
            gap = gene["start"] - previous["end"]

            same_strand = gene["strand"] == previous["strand"]

            if gap <= MAX_GAP and same_strand:
                current.append(gene)
            else:
                if len(current) >= MIN_GENES_CLUSTER:
                    cluster_id += 1
                    clusters.append((cluster_id, current))
                current = [gene]

        if len(current) >= MIN_GENES_CLUSTER:
            cluster_id += 1
            clusters.append((cluster_id, current))

    return annotated_hits, clusters


def classify_operon(systems):
    systems = set(systems)

    if any("ars" in s for s in systems):
        if len(systems) >= 2:
            return "arsenic_mosaic_operon"
        return "arsenic_operon_like"

    if any("czc" in s for s in systems):
        if any("sil" in s or "cus" in s for s in systems):
            return "cation_efflux_mosaic_operon"
        return "czc_like_operon"

    if any("sil" in s or "cus" in s for s in systems):
        return "silver_copper_efflux_operon"

    if any("pco" in s or "cop" in s for s in systems):
        return "copper_oxidation_operon"

    if any("mer" in s for s in systems):
        return "mercury_operon_like"

    return "metal_resistance_operon_like"


all_genes_out = []
all_operons_out = []

for gff_file in glob.glob(f"{ANNOT_DIR}/*/*.gff"):
    group = gff_file.split("/")[-2]
    sample = os.path.basename(gff_file).replace(".gff", "")
    tbl_file = f"{HMM_DIR}/{sample}_pfam.tbl"

    if not os.path.exists(tbl_file):
        print(f"[AVISO] HMMER não encontrado para {sample}: {tbl_file}")
        continue

    print(f"[INFO] Processando {sample}")

    gene_table = parse_gff(gff_file)
    hits = parse_hmmer(tbl_file)

    annotated_hits, clusters = build_operons(sample, gene_table, hits)

    for h in annotated_hits:
        all_genes_out.append([
            group,
            sample,
            h["contig"],
            h["gene_id"],
            h["start"],
            h["end"],
            h["strand"],
            h["pfam_acc"],
            h["pfam_name"],
            h["system"],
            h["evalue"],
            h["score"],
        ])

    for cid, cluster in clusters:
        systems = [g["system"] for g in cluster]
        pfams = [g["pfam_acc"] for g in cluster]
        genes = [g["gene_id"] for g in cluster]

        start = min(g["start"] for g in cluster)
        end = max(g["end"] for g in cluster)
        strand = cluster[0]["strand"]
        contig = cluster[0]["contig"]

        all_operons_out.append([
            group,
            sample,
            f"{sample}_operon_{cid}",
            contig,
            start,
            end,
            end - start + 1,
            strand,
            len(cluster),
            ",".join(genes),
            ",".join(pfams),
            ",".join(systems),
            classify_operon(systems),
        ])


with open(f"{OUT_DIR}/metal_genes_detected.tsv", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "group", "sample", "contig", "gene_id", "start", "end", "strand",
        "pfam_acc", "pfam_name", "system", "evalue", "score"
    ])
    writer.writerows(all_genes_out)


with open(f"{OUT_DIR}/metal_operon_candidates.tsv", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "group", "sample", "operon_id", "contig", "start", "end",
        "length_bp", "strand", "n_genes", "genes", "pfams",
        "systems", "operon_class"
    ])
    writer.writerows(all_operons_out)

print("\n[OK] Arquivos gerados:")
print(f"{OUT_DIR}/metal_genes_detected.tsv")
print(f"{OUT_DIR}/metal_operon_candidates.tsv")
