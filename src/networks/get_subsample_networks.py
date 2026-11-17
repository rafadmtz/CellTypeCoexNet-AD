
from pathlib import Path
import pandas as pd
import polars as pl
import numpy as np
from pathlib import Path
from anndata import AnnData, read_h5ad
from src.utils.numba_mi_calculation import filter_protein_coding_counts, calculate_mi


def get_subsample(adata, donors, n_nuclei, subclass, min_UMIs, adnc):
    
    try:
        subsample_adata = adata[
            (
                adata.obs["ADNC"]
                == adnc
            )
            & (adata.obs["Subclass"] == subclass)
            & (adata.obs["Number of UMIs"] >= min_UMIs),
            :
        ].to_memory()
    finally:
        adata.file.close() 
        
    rng = np.random.default_rng(2)
    
    subsample_adata = subsample_adata[
        subsample_adata.obs["donor_id"].isin(donors)
        ].copy()
     
    nuclei = rng.choice(subsample_adata.n_obs, size=n_nuclei, replace=False)
    subsample_adata = subsample_adata[nuclei].copy()
    
    counts = pd.DataFrame(
        subsample_adata.X.toarray().T,
        index=subsample_adata.var["feature_name"],
        columns=subsample_adata.obs_names,
    )
    
    return counts

def main():
    
    count_path = Path("output/vip_sub_test.tsv")
    mi_path = Path("output/network/mi_matrix_vip_test.parquet")
    edgelist_path = Path("output/network/edgelist_vip_test.parquet")
    
    for path in (count_path, mi_path, edgelist_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    
    ANNDATA = read_h5ad("data/raw/MTG_snRNAseq.h5ad", backed="r")
    
    annotation = pd.read_csv("data/raw/gene_annotation_from_anndata.csv")
    protein_coding_list = annotation.loc[
        annotation["gene_biotype"] == "protein_coding", "gene_name"
    ].to_list()
    
    subclass = "Vip"
    n_nuclei = 100

    rng = np.random.default_rng(2)

    donors = ANNDATA.obs["donor_id"].groupby(
        ANNDATA.obs["ADNC"], observed=True
    ).unique().to_dict()
        
    subsample_donors = rng.choice(donors["High"], size=6, replace=False)
    
    counts_per_donor = []

    for sub_donor in subsample_donors:
        counts_donor = get_subsample(
            ANNDATA, [sub_donor], n_nuclei, subclass,
            min_UMIs=2000, adnc="High"
        )
        counts_per_donor.append(counts_donor)

    counts = pd.concat(counts_per_donor, axis=1, verify_integrity=True)
    
    filter_protein_coding_counts(
        counts,
        protein_coding_list,
        outpath=count_path,
        subsample="vip_sub_test"
        )
    
    calculate_mi(
        count_path=count_path,
        mi_outpath=mi_path,
        edgelist_outpath=edgelist_path
        )
    
if __name__ == "__main__":
    main()