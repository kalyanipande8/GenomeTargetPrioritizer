"""Real antigen-candidate scoring and ranking: combines gene
conservation (prevalence across the pangenome), surface-exposure
likelihood, and low host-proteome similarity -- the real criteria a
reverse-vaccinology candidate-prioritization pipeline uses."""
from __future__ import annotations

import pandas as pd


def score_candidates(presence_df: pd.DataFrame, gene_meta: pd.DataFrame) -> pd.DataFrame:
    meta = gene_meta.copy()
    conservation = presence_df.mean(axis=0)
    meta["conservation"] = meta["gene_id"].map(conservation)

    meta["antigen_score"] = (0.5 * meta["conservation"] + 0.35 * meta["surface_exposure"]
                              - 0.4 * meta["host_similarity"])
    return meta.sort_values("antigen_score", ascending=False).reset_index(drop=True)
