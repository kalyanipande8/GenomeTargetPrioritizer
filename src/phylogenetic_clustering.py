"""Real automated phylogenetic/functional clustering: hierarchical
clustering (average linkage) on the Jaccard distance between genes'
presence/absence patterns across genomes -- genes that co-occur in the
same genomes cluster together, the real proxy a pangenome pipeline
uses for orthologous-family assignment without full tree-building."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist


def cluster_genes(presence_df: pd.DataFrame, n_clusters: int = 18) -> pd.Series:
    gene_matrix = presence_df.values.T  # genes x genomes
    distances = pdist(gene_matrix, metric="jaccard")
    Z = linkage(distances, method="average")
    cluster_labels = fcluster(Z, t=n_clusters, criterion="maxclust")
    return pd.Series(cluster_labels, index=presence_df.columns, name="cluster")
