
"""
Download the SEA-AD MTG snRNA-seq dataset from CELLxGENE and export its
cell metadata and gene annotation as CSV files.

Source: CELLxGENE Discover, SEA-AD collection, whole MTG taxonomy (h5ad).

Outputs:
    data/raw/MTG_snRNAseq.h5ad
    data/raw/MTG_snRNAseq_metadata.csv
    data/raw/gene_annotation_from_anndata.csv

Usage:
    python get_MTG_snRNAseq_anndata.py
"""

from pathlib import Path
import urllib.request
import anndata
import pandas as pd

MTG_ANNDATA="https://datasets.cellxgene.cziscience.com/732b4aa7-e441-4196-b05e-127433425613.h5ad"

def get_MTG_snRNAseq_anndata():
    """
    Download the MTG snRNA-seq anndata object from CELLxGENE.

    The file is downloaded only if it does not already exist. The download
    is written to a temporary .part file and renamed when complete, so an
    interrupted download is never mistaken for a finished one.

    Returns:

    Path
        Path to the downloaded h5ad file (data/raw/MTG_snRNAseq.h5ad).
    """
    
    anndata_path = Path("data/raw/MTG_snRNAseq.h5ad")
    if not anndata_path.exists():
        print(f"Downloading MTG snRNA-seq anndata from {MTG_ANNDATA}...")
        anndata_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = anndata_path.with_suffix(".h5ad.part")
        urllib.request.urlretrieve(MTG_ANNDATA, tmp_path)
        tmp_path.rename(anndata_path)
        print(f"Downloaded to {anndata_path}.")
    else:
        print(f"MTG snRNA-seq anndata already exists at {anndata_path}.")
    return anndata_path

def get_metadata_from_anndata(anndata_path):
    """
    Extract cell metadata from the anndata object and save it as a CSV file.

    The h5ad file is opened in backed mode, so the expression matrix is not
    loaded into memory.

    Parameters:

    anndata_path : Path
        Path to the MTG snRNA-seq h5ad file.

    Returns:

    Path
        Path to the CSV with one row per nucleus and its metadata columns
        (e.g. donor_id, ADNC, Subclass).
    """

    adata = anndata.read_h5ad(anndata_path, backed= "r")
    
    adata = adata.obs.loc[
        adata.obs["ADNC"] != "Reference"
        ]
    
    metadata = adata.obs
    metadata_path = Path("data/raw/MTG_snRNAseq_metadata.csv")
    metadata.to_csv(metadata_path)
    print(f"Metadata saved to {metadata_path}.")
    return metadata_path

def get_gene_annotation_from_anndata(anndata_path):
    """
    Extract gene annotation from the anndata object and save it as a CSV file.

    The h5ad file is opened in backed mode, so the expression matrix is not
    loaded into memory.

    Parameters:

    anndata_path : Path
        Path to the MTG snRNA-seq h5ad file.

    Returns:
    
    Path
        Path to the CSV with Ensembl ID (gene_id), gene symbol (gene_name)
        and biotype (gene_biotype) for each gene.
    """
    
    adata = anndata.read_h5ad(anndata_path, backed="r")
    gene_annotation = pd.DataFrame({
        "gene_id": adata.var_names,
        "gene_name": adata.var["feature_name"].values,
        "gene_biotype": adata.var["feature_type"].values
    })
    gene_annotation_path = Path("data/raw/gene_annotation_from_anndata.csv")
    gene_annotation_path.parent.mkdir(parents=True, exist_ok=True)
    gene_annotation.to_csv(gene_annotation_path, index=False)
    print(f"Gene annotation saved to {gene_annotation_path}.")
    return gene_annotation_path

def main():
    """
    Download the MTG snRNA-seq dataset and export its cell metadata and
    gene annotation as CSV files.
    """
    anndata_path= get_MTG_snRNAseq_anndata()
    get_metadata_from_anndata(anndata_path)
    get_gene_annotation_from_anndata(anndata_path)
    
if __name__ == "__main__":
    main()