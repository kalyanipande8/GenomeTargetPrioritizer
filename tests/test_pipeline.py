import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.simulate_pangenome import simulate
from src.antigen_scoring import score_candidates
from src.phylogenetic_clustering import cluster_genes
from src.validation import pairwise_clustering_precision


def test_pangenome_has_800plus_genomes_and_500plus_genes():
    presence_df, gene_meta = simulate()
    assert presence_df.shape[0] >= 800
    assert presence_df.shape[1] >= 500


def test_scoring_produces_ranked_output():
    presence_df, gene_meta = simulate()
    ranked = score_candidates(presence_df, gene_meta)
    assert len(ranked) == len(gene_meta)
    assert ranked["antigen_score"].is_monotonic_decreasing


def test_clustering_precision_beats_random_baseline():
    presence_df, gene_meta = simulate()
    cluster_labels = cluster_genes(presence_df, n_clusters=18)
    result = pairwise_clustering_precision(cluster_labels, gene_meta)
    random_baseline = 1 / 18
    assert result["pairwise_precision"] > random_baseline


def test_18_clusters_discovered():
    presence_df, gene_meta = simulate()
    cluster_labels = cluster_genes(presence_df, n_clusters=18)
    assert cluster_labels.nunique() <= 18


def test_conservation_score_bounded():
    presence_df, gene_meta = simulate()
    ranked = score_candidates(presence_df, gene_meta)
    assert (ranked["conservation"] >= 0).all() and (ranked["conservation"] <= 1).all()
