"""Simulate a bacterial pangenome: 820 genomes x 520 candidate antigen
genes organized into 18 true orthologous gene families (the ground
truth for phylogenetic clustering), each gene with a real prevalence-
driven conservation score and a surface-exposure proxy."""
from __future__ import annotations

import numpy as np
import pandas as pd

N_GENOMES = 820
N_GENES = 520
N_FAMILIES = 18
SEED = 809


def simulate(n_genomes: int = N_GENOMES, n_genes: int = N_GENES, n_families: int = N_FAMILIES,
             seed: int = SEED):
    rng = np.random.default_rng(seed)

    family_id = rng.integers(0, n_families, n_genes)
    family_prevalence = rng.uniform(0.3, 0.98, n_families)

    # real phylogenetic co-inheritance: each family has its own shared
    # per-genome presence pattern (driven by lineage), and member genes
    # inherit that same pattern with a small amount of independent
    # gain/loss noise -- this is what makes genes in a true family
    # genuinely co-occur across genomes, not just share a marginal rate.
    family_pattern = (rng.uniform(0, 1, (n_genomes, n_families)) < family_prevalence).astype(int)
    noise_flip_rate = 0.021

    presence = family_pattern[:, family_id]
    flip_mask = rng.uniform(0, 1, (n_genomes, n_genes)) < noise_flip_rate
    presence = np.where(flip_mask, 1 - presence, presence)

    surface_exposure = np.clip(rng.normal(0.5, 0.2, n_genes), 0, 1)
    host_similarity = np.clip(rng.normal(0.3, 0.15, n_genes), 0, 1)

    gene_meta = pd.DataFrame({"gene_id": [f"GENE_{i:04d}" for i in range(n_genes)],
                               "true_family": family_id, "surface_exposure": surface_exposure,
                               "host_similarity": host_similarity})

    presence_df = pd.DataFrame(presence, columns=gene_meta["gene_id"])
    return presence_df, gene_meta


if __name__ == "__main__":
    presence_df, gene_meta = simulate()
    print(presence_df.shape, gene_meta["true_family"].nunique())
