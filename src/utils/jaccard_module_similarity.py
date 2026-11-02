import numpy as np
import pandas as pd
from pathlib import Path
from itertools import product, combinations


def jaccard_modules(list1: list, list2: list) -> float:
    """
    Compute the Jaccard index between two modules.

    Parameters
    ----------
    list1 : list
        Genes of the first module.
    list2 : list
        Genes of the second module.

    Returns
    -------
    float
        Shared genes divided by total unique genes. 0 if both lists are empty.
    """
    set1, set2 = set(list1), set(list2)
    union = set1 | set2
    if not union:
        return 0.0
    return len(set1 & set2) / len(union)


def summarize_similarity(similarity_matrix: np.ndarray, match_threshold: float = 0.5) -> dict:
    """
    Summarize a module similarity matrix into four scores.

    Rows are reference modules and columns are compared modules. For each
    reference module, the nonzero Jaccard values are averaged; modules with
    no overlap count as 0.

    Parameters
    ----------
    similarity_matrix : np.ndarray
        Jaccard index between every pair of modules (reference x compared).
    match_threshold : float
        Minimum Jaccard index to count a module as a strong match. Default 0.5.

    Returns
    -------
    dict
        similarity_score (mean over reference modules), matching_modules
        (reference modules with any overlap), strong_matches (reference
        modules with Jaccard >= match_threshold) and perfect_matches
        (reference modules with an identical module).
    """
    row_similarity_scores = []

    for row in similarity_matrix:
        overlap_values = row[row > 0]
        if overlap_values.size == 0:
            row_similarity_scores.append(0)
        else:
            row_similarity_scores.append(np.mean(overlap_values))

    similarity_score = np.mean(row_similarity_scores)

    matching_modules = (similarity_matrix > 0).any(axis=1).sum()
    strong_matches = (similarity_matrix >= match_threshold).any(axis=1).sum()
    perfect_matches = (similarity_matrix == 1).any(axis=1).sum()

    return {
        "similarity_score": similarity_score,
        "matching_modules": matching_modules,
        "strong_matches": strong_matches,
        "perfect_matches": perfect_matches,
    }


def similarity_matrix(reference_modules: list, compared_modules: list) -> np.ndarray:
    """
    Compute the Jaccard index between every pair of modules of two groups.

    Parameters
    ----------
    reference_modules : list
        Modules of the reference group (matrix rows), as sets of genes.
    compared_modules : list
        Modules of the compared group (matrix columns), as sets of genes.

    Returns
    -------
    np.ndarray
        Jaccard matrix of shape (len(reference_modules), len(compared_modules)).
    """
    n_reference = len(reference_modules)
    n_compare = len(compared_modules)

    matrix = np.zeros((n_reference, n_compare))

    for i, j in product(range(n_reference), range(n_compare)):
        matrix[i, j] = jaccard_modules(reference_modules[i], compared_modules[j])

    return matrix


def module_similarity_between_groups(groups: list, communities: dict, min_size: int,
                                     outpath: Path, match_threshold: float = 0.5) -> pd.DataFrame:
    """
    Compare the modules of every group against every other group and save them.

    Modules with fewer than min_size genes are dropped, so small modules do
    not weigh as much as large ones in similarity_score. Each pair of groups
    is computed once and summarized in both directions.

    Parameters
    ----------
    groups : list
        Groups to compare. Each one is used once as reference.
    communities : dict
        {group: {community_id: [genes]}}.
    min_size : int
        Minimum number of genes for a module to be compared.
    outpath : Path
        Output CSV path.
    match_threshold : float
        Minimum Jaccard index for a strong match. Default 0.5.

    Returns
    -------
    pd.DataFrame
        One row per (reference, compared) pair: reference, compared,
        similarity_score, matching_modules, strong_matches, perfect_matches.
    """
    modules = {}

    for group in groups:
        modules[group] = [
            set(genes) for genes in communities[group].values()
            if len(genes) >= min_size
        ]

    rows = []

    for group1, group2 in combinations(groups, 2):

        matrix = similarity_matrix(
            reference_modules=modules[group1],
            compared_modules=modules[group2]
        )

        rows.append({
            "reference": group1,
            "compared": group2,
            **summarize_similarity(matrix, match_threshold)
        })

        rows.append({
            "reference": group2,
            "compared": group1,
            **summarize_similarity(matrix.T, match_threshold)
        })

    df_similarity = pd.DataFrame(rows)

    df_similarity.to_csv(outpath, index=False)

    return df_similarity