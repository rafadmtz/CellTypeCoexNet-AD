from pathlib import Path
import urllib.request
import pandas as pd
from anndata import AnnData, read_h5ad

ANNDATA = Path("data/raw/MTG_snRNAseq.h5ad")

DATA_PROCESSED_DIR= Path("data/processed/counts")
DATA_PROCESSED_DIR.parent.mkdir(parents=True, exist_ok=True)

def subclasses_by_donor(adata: AnnData) -> Path:
    """
    Count nuclei per cell subclass and donor in the MTG snRNA-seq dataset.

    Nuclei from donors with ADNC == "Reference" are excluded. Two CSV files
    are written to "DATA_PROCESSED_DIR":

    Parameters

    adata : AnnData
        MTG snRNA-seq object.

    Returns

    Path
        Path to "subclasses_counts_by_donor.csv".
    """

    subclasses = adata.obs.loc[
        adata.obs["ADNC"] != "Reference",
        ["Subclass", "donor_id"]
    ]

    matrix_subclasses = pd.crosstab(
        subclasses["Subclass"],
        subclasses["donor_id"]
    )

    stats = matrix_subclasses.agg(
        ["count", "mean", "median", "std", "min", "max"],
        axis=1
    )

    stats.to_csv(DATA_PROCESSED_DIR / "subclasses_counts_stats.csv", index=True)

    matrix_subclasses.to_csv(DATA_PROCESSED_DIR / "subclasses_counts_by_donor.csv", index=True)
    
    return Path(DATA_PROCESSED_DIR / "subclasses_counts_by_donor.csv")
    
def supertypes_by_donor(adata: AnnData) -> Path:
    """
    Count nuclei per cell supertypes and donor in the MTG snRNA-seq dataset.

    Nuclei from donors with ADNC == "Reference" are excluded. Two CSV files
    are written to "DATA_PROCESSED_DIR":

    Parameters

    adata : AnnData
        MTG snRNA-seq object.

    Returns

    Path
        Path to "supertypes_counts_by_donor.csv".
    """
    
    supertypes = adata.obs.loc[
        adata.obs["ADNC"] != "Reference",
        ["Supertypes", "donor_id"]
    ]
    
    matrix_supertypes = pd.crosstab(
        supertypes["Supertypes"],
        supertypes["donor_id"]
    )
    
    stats = matrix_supertypes.agg(
        ["count", "mean", "median", "std", "min", "max"],
        axis=1
    )
    
    stats.to_csv(DATA_PROCESSED_DIR / "supertypes_counts_stats.csv", index=True)
    
    matrix_supertypes.to_csv(DATA_PROCESSED_DIR / "supertypes_counts_by_donor.csv", index=True)
    

    
def main():
    adata = read_h5ad(ANNDATA, backed="r")
    subclasses_by_donor(adata)
    supertypes_by_donor(adata)