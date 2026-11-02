import time
import polars as pl
from numba_mi import pipeline
import pandas as pd


def filter_protein_coding_counts(count_path: str, protein_coding: list, outpath: str) -> None:
    """
    Keep only protein-coding genes in a count matrix and save it as TSV.

    Versioned duplicates (e.g. GENE.1 when GENE also exists) are collapsed
    into the base name by summing their counts. Genes follow the order of
    protein_coding.

    Parameters
    ----------
    count_path : str
        TSV count matrix with genes as rows.
    protein_coding : list
        Protein-coding gene names to keep.
    outpath : str
        Output TSV with the filtered count matrix.
    """
    gene_count = pd.read_csv(count_path, sep="\t", index_col=0)

    base_names = gene_count.index.str.replace(r"\.\d+$", "", regex=True)
    is_duplicate = (base_names != gene_count.index) & base_names.isin(gene_count.index)

    if is_duplicate.any():
        print(f"Collapsing duplicated genes in {count_path}: {list(gene_count.index[is_duplicate])}")
        gene_count.index = gene_count.index.where(~is_duplicate, base_names)
        gene_count = gene_count.groupby(level=0, sort=False).sum()

    genes = [g for g in protein_coding if g in gene_count.index]

    gene_count = gene_count.loc[genes]
    gene_count.index.name = "gene"

    print(f"Shape of {count_path}:")
    print(gene_count.shape)

    gene_count.to_csv(outpath, sep="\t", index_label="gene")

    print(f"Saved in: {outpath}")


def matrix_to_edgelist(path: str, outpath: str) -> pl.DataFrame:
    """
    Convert a square MI matrix into an edgelist and save it as parquet.

    Each gene pair is kept once (source < target). Rows coming from index
    labels ("Unnamed: 0", "gene", "") are removed.

    Parameters
    ----------
    path : str
        Parquet MI matrix whose first column holds the gene names.
    outpath : str
        Output parquet edgelist.

    Returns
    -------
    pl.DataFrame
        Edgelist with columns source, target and mi.
    """
    adj_matrix = pl.read_parquet(path)
    col_genes = adj_matrix.columns[0]

    edgelist = adj_matrix.unpivot(index=col_genes, variable_name="target", value_name="mi")
    edgelist = edgelist.rename({col_genes: "source"})
    edgelist = edgelist.filter(pl.col("source") < pl.col("target"))

    n_edges_before = len(edgelist)

    index_labels = ["Unnamed: 0", "gene", ""]
    edgelist = edgelist.filter(
        ~pl.col("source").is_in(index_labels) & ~pl.col("target").is_in(index_labels)
    )
    if n_edges_before != len(edgelist):
        print(f" {n_edges_before - len(edgelist)} edges from index labels removed")

    edgelist.write_parquet(outpath, compression="zstd")

    print(f"Edgelist saved in: {outpath}")
    return edgelist


def calculate_mi(count_path: str, mi_outpath: str, edgelist_outpath: str) -> None:
    """
    Calculate the gene-gene MI matrix of a count matrix and its edgelist
    with the numba_mi pipeline.
 
    Parameters
    ----------
    count_path : str
        TSV count matrix with genes as rows.
    mi_outpath : str
        Output parquet with the MI matrix.
    edgelist_outpath : str
        Output parquet with the MI edgelist.
    """
    print(f"{count_path}: calculating MI...")
    start = time.time()
    mi_matrix = pipeline.whole_process(
        file_path=count_path,
        axis=1,
        sep="\t",
        index_col=0,
        num_threads=40,
    )
    print(f"  done in {time.time() - start:.1f}s — shape: {mi_matrix.shape}")
 
    mi_matrix.write_parquet(mi_outpath, compression="zstd")
    print(f"  saved in: {mi_outpath}")
 
    matrix_to_edgelist(
        path=mi_outpath,
        outpath=edgelist_outpath
    )
 