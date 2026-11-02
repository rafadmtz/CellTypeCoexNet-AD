import pandas as pd
from pathlib import Path
import polars as pl
import networkx as nx
import igraph as ig
import matplotlib.pyplot as plt
from itertools import combinations


def read_graphs(paths: dict) -> dict:
    """
    Read MI edgelists and build one undirected graph per group.

    Parameters
    ----------
    paths : dict
        {group: Path} to parquet edgelists with columns source, target and mi.

    Returns
    -------
    dict
        {group: nx.Graph} with mi stored as an edge attribute.
    """
    edgelists = {
        group: pl.read_parquet(
            path
        )
        for group, path in paths.items()
        }
    graphs = {
        group: nx.from_pandas_edgelist(
            df.to_pandas(), source="source", target="target", edge_attr="mi"
        )
        for group, df in edgelists.items()
        }
    return graphs

def network_metrics(G: nx.Graph) -> dict:
    """
    Calculate basic topological metrics of a network.

    Parameters
    ----------
    G : nx.Graph
        Coexpression network.

    Returns
    -------
    dict
        nodes, edges, connected_components and clustering_coefficient
        (average clustering, nodes with degree < 2 count as 0).
    """
    components = list(nx.connected_components(G))
    return {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "connected_components": len(components),
        "clustering_coefficient": nx.average_clustering(G),
    }

def jaccard_index(G1: nx.Graph, G2: nx.Graph) -> float:
    """
    Calculate the Jaccard index between the edge sets of two networks.

    Edges are compared as unordered node pairs. 1 means identical edge sets,
    0 means no shared edges.

    Parameters
    ----------
    G1, G2 : nx.Graph
        Networks to compare.

    Returns
    -------
    float
        Shared edges divided by the union of edges.
    """
    e1 = {frozenset(e) for e in G1.edges()}
    e2 = {frozenset(e) for e in G2.edges()}
    return len(e1 & e2) / len(e1 | e2)


def networks_properties(graphs: dict, groups: list, properties_outpath: Path, jaccard_outpath: Path) -> None:
    """
    Calculate network properties per group and pairwise edge Jaccard indices.

    The average shortest path length is computed on the giant component
    with igraph. Both tables are written to CSV.

    Parameters
    ----------
    graphs : dict
        {group: nx.Graph}.
    groups : list
        Groups to analyze, in output order.
    properties_outpath : Path
        Output CSV with the network properties per group.
    jaccard_outpath : Path
        Output CSV with the pairwise Jaccard indices.
    """
    results = []

    for group in groups:
        G = graphs[group]

        largest_cc = max(nx.connected_components(G), key=len)
        G_giant = G.subgraph(largest_cc)

        g_giant = ig.Graph.from_networkx(G_giant)
        avg_path_length = g_giant.average_path_length(directed=False)

        row = {
            "group": group,
            **network_metrics(G),
            "giant_component_nodes": G_giant.number_of_nodes(),
            "avg_shortest_path_length": avg_path_length,
        }
        results.append(row)

    df_results = pd.DataFrame(results)
    df_results.to_csv(properties_outpath, index=False,)

    jaccard_results= []

    for i, k in combinations(groups, 2):
        jac = jaccard_index(graphs[i], graphs[k])
        jaccard_results.append({
            "group_1": i,
            "group_2": k,
            "jaccard_index": jac,
        })
        print(f"{i} vs {k}: jaccard = {jac:.4f}")

    df_jaccard = pd.DataFrame(jaccard_results)
    df_jaccard.to_csv(jaccard_outpath, index=False)



def degree_change(groups: list, outdatadir: Path, outpath: Path) -> None:
    """
    Merge gene degrees across groups and rank genes by degree change.

    Genes missing from a group get degree 0. Genes are sorted by max_range,
    the difference between their highest and lowest degree across groups.

    Parameters
    ----------
    groups : list
        Groups to merge, in column order.
    outdatadir : Path
        Folder with the gene_degrees_{group}_protein_coding.csv files
        (columns gene and degree).
    outpath : Path
        Output CSV: gene, one degree column per group, max_range.
    """
    degree_tables = {
        group: pd.read_csv(
            outdatadir / f"gene_degrees_{group}_protein_coding.csv"
        )
        for group in groups
    }

    all_genes = set()
    for group in groups:
        all_genes.update(degree_tables[group]["gene"])
    print(f"all genes len {len(all_genes)}")

    df_degree = pd.DataFrame({"gene": sorted(all_genes)})

    for group in groups:
        df_group = degree_tables[group].rename(columns={"degree": group})
        df_degree = df_degree.merge(df_group, on="gene", how="left")

    df_degree[groups] = df_degree[groups].fillna(0)

    df_degree["max_range"] = df_degree[groups].max(axis=1) - df_degree[groups].min(axis=1)
    df_degree = df_degree.sort_values("max_range", ascending=False).reset_index(drop=True)

    df_degree.to_csv(outpath, index=False)



def mi_distribution_plot(groups: list, paths: dict, outfigdir: Path) -> None:
    """
    Plot the MI distribution of each group as log-scale histograms.

    Draws a 2x2 grid with shared axes, so only the first four groups are
    plotted.

    Parameters
    ----------
    groups : list
        Groups to plot, in panel order.
    paths : dict
        {group: Path} to parquet edgelists with an mi column.
    outfigdir : Path
        Folder where mi_distribution_4groups_protein_coding.png is saved.
    """
    colors = ["#abd9e9", "#2c7bb6", "#fdae61", "#d7191c"]

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
    axes = axes.flatten()

    for ax, group, color in zip(axes, groups, colors):
        edgelist = pl.read_parquet(paths[group])
        ax.hist(edgelist["mi"].to_numpy(), bins=100, color=color, edgecolor="black", linewidth=0.3)
        ax.set_yscale("log")
        ax.set_title(group)
        ax.set_xlabel("MI")
        ax.set_ylabel("Frequency")

    fig.suptitle("MI distribution by severity group")
    fig.tight_layout()
    fig.savefig(outfigdir / "mi_distribution_4groups_protein_coding.png", dpi=300)
    plt.close(fig)