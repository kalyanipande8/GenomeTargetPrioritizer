"""Validates clustering against the 18 true gene families via
pairwise clustering precision: of all gene pairs placed in the same
discovered cluster, what fraction are genuinely from the same true
family -- the standard pair-counting precision metric for clustering
quality, computed against ground truth the clustering never saw."""
from __future__ import annotations

import numpy as np
import pandas as pd
from itertools import combinations


def pairwise_clustering_precision(cluster_labels: pd.Series, gene_meta: pd.DataFrame, max_pairs: int = 200000) -> dict:
    meta = gene_meta.set_index("gene_id")
    clusters = cluster_labels.groupby(cluster_labels).groups

    same_cluster_pairs = 0
    same_family_and_cluster = 0
    rng = np.random.default_rng(809)

    for cluster_id, genes in clusters.items():
        genes = list(genes)
        if len(genes) < 2:
            continue
        pairs = list(combinations(genes, 2))
        if len(pairs) > max_pairs // max(1, len(clusters)):
            idx = rng.choice(len(pairs), max_pairs // max(1, len(clusters)), replace=False)
            pairs = [pairs[i] for i in idx]
        for g1, g2 in pairs:
            same_cluster_pairs += 1
            if meta.loc[g1, "true_family"] == meta.loc[g2, "true_family"]:
                same_family_and_cluster += 1

    precision = round(same_family_and_cluster / same_cluster_pairs, 4) if same_cluster_pairs else None
    return {"n_discovered_clusters": cluster_labels.nunique(), "n_true_families": meta["true_family"].nunique(),
            "n_same_cluster_pairs_evaluated": same_cluster_pairs, "pairwise_precision": precision}
