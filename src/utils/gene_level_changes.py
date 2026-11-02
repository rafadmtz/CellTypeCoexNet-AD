import pandas as pd
import numpy as np
import polars as pl
import networkx as nx
from pathlib import Path
from scipy.stats import hypergeom, false_discovery_control

def read_communities(path: Path) -> dict:
    """
    Read the communities detected with Infomap into a dictionary.

    Parameters
    ----------
    path : Path
        CSV file with columns gene and community.

    Returns
    -------
    dict
        {community_id: [genes]}.
    """
    df_communities = pd.read_csv(path)

    n_communities= df_communities["community"].unique()

    communities = {i: df_communities.loc[df_communities["community"] ==i, "gene"].tolist()
                   for i in n_communities}

    return communities

def neighbors_by_state(gene: str, G: nx.Graph) -> set:
    """
    Get the neighbors of a gene in a network.

    Parameters
    ----------
    gene : str
        Gene name.
    G : nx.Graph
        Coexpression network.

    Returns
    -------
    set
        Neighbors of gene. Empty if gene is not in G.
    """
    if G.has_node(gene):
        return set(G.neighbors(gene))
    return set()

def gene_neighborhood_change(graphs: dict, contrasts: list, outpath: Path) -> pd.DataFrame:
    """
    Compare each gene's neighborhood between pairs of networks and test
    whether its neighbors are conserved more than expected by chance.

    For genes in both networks, a one-sided hypergeometric test is computed on
    the neighbors restricted to nodes shared by both networks, with
    Benjamini-Hochberg correction within each contrast. A low p_adj means the
    gene keeps more neighbors than expected by chance.

    Parameters
    ----------
    graphs : dict
        {group: nx.Graph}.
    contrasts : list
        [(group1, group2)]. group1 is the reference.
    outpath : Path
        Output CSV path.

    Returns
    -------
    pd.DataFrame
        One row per gene and contrast. status is in_both, in (only in group2)
        or out (only in group1). expected, p_value and p_adj are NaN for genes
        not in both networks.
    """
    results = []

    for group1, group2 in contrasts:
        G1, G2 = graphs[group1], graphs[group2]
        nodes1, nodes2 = set(G1.nodes), set(G2.nodes)
        shared_nodes_G = nodes1 & nodes2
        genes = sorted(nodes1 | nodes2)

        for gene in genes:
            gene_set1 = neighbors_by_state(gene, G1)
            gene_set2 = neighbors_by_state(gene, G2)

            if gene in nodes1 and gene in nodes2:
                status = "in_both"
            elif gene in nodes2:
                status = "in"
            else:
                status = "out"

            conserved = gene_set1 & gene_set2
            gained = gene_set2 - gene_set1
            lost = gene_set1 - gene_set2
            union = gene_set1 | gene_set2

            results.append({
                "gene": gene,
                "group1": group1,
                "group2": group2,
                "status": status,
                "degree1": len(gene_set1),
                "degree2": len(gene_set2),
                "n_universe": len(genes),
                "conserved_neighbors": sorted(conserved),
                "gained": sorted(gained),
                "lost": sorted(lost),
                "n_conserved_neighbors": len(conserved),
                "net_change_degree": len(gained) - len(lost),
                "edges_changed": len(gained) + len(lost),
                "jaccard": len(conserved) / len(union) if union else None,
                "degree1_shared_G": len(gene_set1 & shared_nodes_G),
                "degree2_shared_G": len(gene_set2 & shared_nodes_G),
                "n_shared_nodes_G": len(shared_nodes_G),
            })

    df = pd.DataFrame(results)

    both = df["status"] == "in_both"
    M = df.loc[both, "n_shared_nodes_G"] - 1
    n = df.loc[both, "degree1_shared_G"]
    N = df.loc[both, "degree2_shared_G"]
    k = df.loc[both, "n_conserved_neighbors"]

    df["expected"] = np.nan
    df["p_value"] = np.nan
    df.loc[both, "expected"] = n * N / M
    df.loc[both, "p_value"] = hypergeom.sf(k - 1, M, n, N)

    df["p_adj"] = np.nan
    df.loc[both, "p_adj"] = (
        df[both]
        .groupby(["group1", "group2"])["p_value"]
        .transform(false_discovery_control)
    )

    df.to_csv(outpath, index=False)
    print(f"Gene neighbors change analysis completed. Results saved to {outpath}.")
    return df

def get_intracommunity_metrics(communities: dict, G: nx.Graph, min_size: int) -> pd.DataFrame:
    """
    Calculate intracommunity degree statistics for each community.

    Parameters
    ----------
    communities : dict
        {community_id: [genes]}.
    G : nx.Graph
        Coexpression network.
    min_size : int
        Minimum number of genes for a community to be included.

    Returns
    -------
    pd.DataFrame
        One row per community: community, size, mean_degree_in, std_degree_in.

    Notes
    -----
    Used by gene_roles to compute the within-module z-score (Guimerà & Amaral, 2005).
    """
    communities_minsz = {
        community: set(genes)
        for community, genes in communities.items()
        if len(genes) >= min_size
    }

    results = []

    for community, com_genes in communities_minsz.items():

        intracommunity_degrees = [
            sum(neighbor in com_genes for neighbor in G.neighbors(com_gene))
            for com_gene in com_genes
        ]

        results.append({
            "community": community,
            "size": len(com_genes),
            "mean_degree_in": np.mean(intracommunity_degrees),
            "std_degree_in": np.std(intracommunity_degrees, ddof=0)
        })

    return pd.DataFrame(results)


def gene_roles(G: nx.Graph, communities: dict, group: str, min_size: int) -> pd.DataFrame:
    """
    Calculate the within-module degree z-score (z) and the participation
    coefficient (P) of each gene.

    z measures how connected a gene is inside its own community compared with
    the other genes of that community; a high z marks a module hub. P measures
    how its links are spread across communities: near 0 if all its neighbors
    are in its own community, near 1 if it connects many communities.

    Parameters
    ----------
    G : nx.Graph
        Coexpression network.
    communities : dict
        {community_id: [genes]}. Every gene in G must belong to one community.
    group : str
        Name of the group the network belongs to.
    min_size : int
        Minimum community size to compute z.

    Returns
    -------
    pd.DataFrame
        One row per gene: group, gene, degree, community, community_size,
        community_degree, z, P. z is NaN if the community is smaller than
        min_size or has zero standard deviation. P is NaN if degree is 0.

    Raises
    ------
    ValueError
        If a gene in G has no assigned community.

    Notes
    -----
    Roles as defined by Guimerà & Amaral (2005).
    """

    intracommunity_metrics = get_intracommunity_metrics(communities, G, min_size)

    intracommunity_metrics = intracommunity_metrics.set_index("community")

    gene_to_community = {
        gene: community
        for community, community_genes in communities.items()
        for gene in community_genes
    }

    for gene in G.nodes():
        if gene not in gene_to_community:
            raise ValueError(f"Gene {gene} has no assigned community.")

    results = []

    for gene in G.nodes():

        community = gene_to_community[gene]

        degree = G.degree(gene)

        neighbors_by_community = {}

        for neighbor in G.neighbors(gene):
            neighbor_community = gene_to_community[neighbor]
            neighbors_by_community[neighbor_community] = (
                neighbors_by_community.get(neighbor_community, 0) + 1
            )

        community_degree = neighbors_by_community.get(community, 0)

        z = np.nan

        if community in intracommunity_metrics.index:

            mean_community_degree = intracommunity_metrics.loc[community, "mean_degree_in"]

            std_community_degree = intracommunity_metrics.loc[community, "std_degree_in"]

            if std_community_degree > 0:
                z = (community_degree - mean_community_degree) / std_community_degree

        P = np.nan

        if degree > 0:
            squared_sum = 0

            for neighbors_in_community in neighbors_by_community.values():
                squared_sum += (neighbors_in_community / degree) ** 2

            P = 1 - squared_sum

        results.append({
            "group": group,
            "gene": gene,
            "degree": degree,
            "community": community,
            "community_size": len(communities[community]),
            "community_degree": community_degree,
            "z": z,
            "P": P
        })

    return pd.DataFrame(results)

def roles_all_groups(graphs: dict, communities: dict, min_size: int, outpath: Path) -> pd.DataFrame:
    """
    Calculate gene roles (z and P) for every group and save them.

    Parameters
    ----------
    graphs : dict
        {group: nx.Graph}.
    communities : dict
        {group: {community_id: [genes]}}. Same keys as graphs.
    min_size : int
        Minimum community size to compute z.
    outpath : Path
        Output CSV path.

    Returns
    -------
    pd.DataFrame
        gene_roles output for all groups concatenated.
    """
    all_roles = []

    for group, G in graphs.items():

        df_roles = gene_roles(
            G=G,
            communities=communities[group],
            group=group,
            min_size=min_size
        )

        all_roles.append(df_roles)

    df_roles_all_groups = pd.concat(all_roles, ignore_index=True)

    df_roles_all_groups.to_csv(outpath, index=False)

    return df_roles_all_groups

def merge_gene_metrics(df_neighbors_change: pd.DataFrame, df_roles: pd.DataFrame, outpath: Path) -> pd.DataFrame:
    """
    Add the z and P of each gene in group1 and group2 to the neighborhood
    change table.

    Parameters
    ----------
    df_neighbors_change : pd.DataFrame
        Output of gene_neighborhood_change.
    df_roles : pd.DataFrame
        Output of roles_all_groups.
    outpath : Path
        Output CSV path.

    Returns
    -------
    pd.DataFrame
        One row per gene and contrast with z1, P1 (group1) and z2, P2 (group2).
        z1 and P1 are NaN for in genes; z2 and P2 are NaN for out genes.
    """
    roles = df_roles[["group", "gene", "z", "P"]]

    df = (
        df_neighbors_change
        .merge(roles.rename(columns={"group": "group1", "z": "z1", "P": "P1"}),
               on=["gene", "group1"], how="left")
        .merge(roles.rename(columns={"group": "group2", "z": "z2", "P": "P2"}),
               on=["gene", "group2"], how="left")
    )

    df.to_csv(outpath, index=False)

    return df