import json
import os

from src.simulate_pangenome import simulate
from src.antigen_scoring import score_candidates
from src.phylogenetic_clustering import cluster_genes
from src.validation import pairwise_clustering_precision


def main():
    presence_df, gene_meta = simulate()
    ranked = score_candidates(presence_df, gene_meta)
    cluster_labels = cluster_genes(presence_df, n_clusters=18)
    precision = pairwise_clustering_precision(cluster_labels, gene_meta)

    summary = {"n_genomes": len(presence_df), "n_candidate_genes": len(gene_meta),
               "top10_candidates": ranked.head(10)["gene_id"].tolist(),
               "clustering_validation": precision}

    os.makedirs("results", exist_ok=True)
    with open("results/summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
