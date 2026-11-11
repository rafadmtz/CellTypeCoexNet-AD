"""
Describe the sequencing depth of every nucleus in the SEA-AD MTG snRNA-seq
dataset and summarize it by donor, by donor and subclass, and by donor and
supertype, so that depth distributions can be plotted and compared across
cell types and disease groups.

Sequencing depth is the total number of UMIs (unique RNA molecules) detected
in a nucleus. The number of detected genes and the mitochondrial fraction are
reported as well. All three come from the cell metadata ("Number of UMIs",
"Genes detected", "Fraction mitochrondrial UMIs"), so the expression matrix
is not needed. With --verify N, N random nuclei are checked against the raw
count matrix (raw/X), reading only those rows from the h5ad file. The metadata
values are slightly higher than raw/X row sums (about 0.5% for UMIs and 1% for
genes) because SEA-AD computed them on its full gene set, while the CELLxGENE
matrix keeps 35,483 genes.

All nuclei are kept, including those from Reference donors. Use the
"Neurotypical reference" or "ADNC" column to exclude them downstream.

Inputs (paths set in config/config.yaml, section sequencing_depth):
    data/raw/MTG_snRNAseq_metadata.parquet  (or .csv) one row per nucleus

Outputs (in results/sequencing_depth/tables/):
    depth_by_nucleus.parquet          one row per nucleus
    depth_by_donor.csv                one row per donor
    depth_by_donor_subclass.csv       one row per donor and subclass
    depth_by_donor_supertype.csv      one row per donor and supertype
    depth_column_dictionary.csv       description of every output column
    run_info.txt                      date, input, counts and checks
    config_used.yaml                  configuration used in this run
    verification_raw_counts.csv       (only with --verify) metadata vs raw/X

Usage (from the repository root):
    python src/preprocessing/sequencing_depth.py
    python src/preprocessing/sequencing_depth.py --verify 20
"""

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

CONFIG = Path("config/config.yaml")

# Depth columns in the CELLxGENE metadata -> short names used in the outputs
DEPTH_COLUMNS = {
    "Number of UMIs": "n_umis",
    "Genes detected": "n_genes",
    "Fraction mitochrondrial UMIs": "frac_mito",  # spelled this way in SEA-AD
}

# Nucleus-level annotations kept for grouping
NUCLEUS_COLUMNS = ["Class", "Subclass", "Supertype"]

# Donor-level annotations (constant across the nuclei of a donor)
DONOR_COLUMNS = [
    "Neurotypical reference", "ADNC", "Cognitive status", "Braak stage",
    "Thal phase", "CERAD score", "APOE4 status", "LATE-NC stage",
    "Age at death", "sex", "PMI",
]

DONOR_ID = "donor_id"
SPECIMEN_ID = "Specimen ID"


def read_metadata(path: Path) -> pd.DataFrame:
    """
    Read the cell metadata from a parquet or CSV file.

    Parameters:

    path : Path
        Path to the metadata file (.parquet or .csv, one row per nucleus).

    Returns:

    pd.DataFrame
        Cell metadata indexed by nucleus ID.
    """
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    if path.suffix == ".csv":
        return pd.read_csv(path, index_col=0)
    raise ValueError(f"Unsupported format (expected .parquet or .csv): {path}")


def depth_by_nucleus(obs: pd.DataFrame) -> pd.DataFrame:
    """
    Build one row per nucleus with its identifiers, annotations and depth metrics.

    Parameters:

    obs : pd.DataFrame
        Cell metadata from the SEA-AD MTG dataset.

    Returns:

    pd.DataFrame
        Columns: donor_id, specimen_id, nucleus and donor annotations,
        n_umis, n_genes, frac_mito, log10_umis, log10_genes, genes_per_umi.
    """
    required = [DONOR_ID, SPECIMEN_ID, *DEPTH_COLUMNS, *NUCLEUS_COLUMNS, *DONOR_COLUMNS]
    missing = [c for c in required if c not in obs.columns]
    if missing:
        raise KeyError(f"Columns not found in the metadata: {missing}")

    depth = pd.DataFrame(index=obs.index)
    depth.index.name = "nucleus_id"
    depth["donor_id"] = obs[DONOR_ID].astype(str)
    depth["specimen_id"] = obs[SPECIMEN_ID].astype(str)
    for column in NUCLEUS_COLUMNS + DONOR_COLUMNS:
        depth[column] = obs[column]
    for original, short in DEPTH_COLUMNS.items():
        depth[short] = pd.to_numeric(obs[original], errors="coerce")

    # Depth is right-skewed; log10 values are easier to plot and compare
    depth["log10_umis"] = np.log10(depth["n_umis"])
    depth["log10_genes"] = np.log10(depth["n_genes"])
    # Library complexity: distinct genes contributed per UMI
    depth["genes_per_umi"] = depth["n_genes"] / depth["n_umis"]
    return depth


def _describe(values: pd.Series, prefix: str) -> dict:
    """Descriptive statistics of one depth metric."""
    return {
        f"{prefix}_mean": values.mean(),
        f"{prefix}_sd": values.std(),
        f"{prefix}_min": values.min(),
        f"{prefix}_q25": values.quantile(0.25),
        f"{prefix}_median": values.median(),
        f"{prefix}_q75": values.quantile(0.75),
        f"{prefix}_max": values.max(),
    }


def summarize(depth: pd.DataFrame, group_by: list) -> pd.DataFrame:
    """
    Summarize depth metrics for each group, with the donor annotations attached.

    Parameters:

    depth : pd.DataFrame
        Output of depth_by_nucleus().
    group_by : list
        Columns that define the groups; must start with "donor_id".

    Returns:

    pd.DataFrame
        One row per group with n_nuclei, n_specimens, UMI and gene statistics
        and the median mitochondrial fraction.
    """
    rows = []
    for key, group in depth.groupby(group_by, observed=True, sort=True):
        key = key if isinstance(key, tuple) else (key,)
        row = dict(zip(group_by, key))
        row["n_nuclei"] = len(group)
        row["n_specimens"] = group["specimen_id"].nunique()
        row.update(_describe(group["n_umis"], "umis"))
        row.update(_describe(group["n_genes"], "genes"))
        row["frac_mito_median"] = group["frac_mito"].median()
        rows.append(row)
    summary = pd.DataFrame(rows)

    donors = depth.groupby("donor_id", observed=True)[DONOR_COLUMNS].first().reset_index()
    summary = summary.merge(donors, on="donor_id", how="left")
    order = group_by + DONOR_COLUMNS + [c for c in summary.columns
                                        if c not in group_by + DONOR_COLUMNS]
    return summary[order]


def quality_checks(depth: pd.DataFrame) -> list:
    """
    Check that the depth metrics and donor annotations are consistent.

    Returns:

    list
        One message per check.
    """
    messages = [f"Missing values in {c}: {int(depth[c].isna().sum())}"
                for c in ["n_umis", "n_genes", "frac_mito"]]
    messages.append(f"Nuclei with 0 UMIs: {int((depth['n_umis'] <= 0).sum())}")
    messages.append("Nuclei with more genes than UMIs (impossible): "
                    f"{int((depth['n_genes'] > depth['n_umis']).sum())}")
    outside = int(((depth["frac_mito"] < 0) | (depth["frac_mito"] > 1)).sum())
    messages.append(f"Nuclei with mitochondrial fraction outside [0, 1]: {outside}")

    per_donor = depth.groupby("donor_id", observed=True)[DONOR_COLUMNS].nunique(dropna=False)
    inconsistent = per_donor.columns[(per_donor > 1).any()].tolist()
    messages.append("Donor annotations that vary within a donor: "
                    f"{inconsistent if inconsistent else 'none'}")
    return messages


def verify_with_raw_counts(depth: pd.DataFrame, n: int, url: str, seed: int = 0) -> pd.DataFrame:
    """
    Compare the metadata depth of n random nuclei with their row in raw/X.

    Only the selected rows are read from the h5ad file.

    Parameters:

    depth : pd.DataFrame
        Output of depth_by_nucleus(), in the same row order as the h5ad file.
    n : int
        Number of nuclei to check.
    url : str
        URL (or local path) of the h5ad file.
    seed : int
        Random seed for the selection.

    Returns:

    pd.DataFrame
        One row per checked nucleus with UMIs and genes from both sources.
    """
    import fsspec
    import h5py

    rng = np.random.default_rng(seed)
    rows = np.sort(rng.choice(len(depth), size=n, replace=False))
    checks = []
    with fsspec.open(url, "rb", block_size=2**20, cache_type="blockcache") as f:
        with h5py.File(f, "r") as h5:
            raw = h5["raw/X"]
            indptr = raw["indptr"]
            for i in rows:
                start, end = int(indptr[i]), int(indptr[i + 1])
                values = raw["data"][start:end]
                checks.append({
                    "nucleus_id": depth.index[i],
                    "umis_metadata": depth["n_umis"].iloc[i],
                    "umis_raw_counts": float(values.sum()),
                    "genes_metadata": depth["n_genes"].iloc[i],
                    "genes_raw_counts": int((values > 0).sum()),
                })
    return pd.DataFrame(checks)


def summarize_verification(verification: pd.DataFrame, tolerance: float = 0.02) -> list:
    """
    Summarize how closely the metadata depth agrees with the raw/X row sums.

    Exact equality is not expected: SEA-AD computed "Number of UMIs" and
    "Genes detected" on its full gene set, while the CELLxGENE matrix keeps
    35,483 genes, so raw/X row sums are slightly lower. Rows are considered
    aligned when the correlation is close to 1 and ratios are within tolerance.

    Parameters:

    verification : pd.DataFrame
        Output of verify_with_raw_counts(). Ratio columns are added in place.
    tolerance : float
        Maximum relative difference counted as agreement (default 2%).

    Returns:

    list
        One message per summary line.
    """
    verification["ratio_umis"] = verification["umis_raw_counts"] / verification["umis_metadata"]
    verification["ratio_genes"] = verification["genes_raw_counts"] / verification["genes_metadata"]
    messages = []
    for metric in ["umis", "genes"]:
        ratio = verification[f"ratio_{metric}"]
        within = ((ratio - 1).abs() <= tolerance).mean()
        corr = verification[f"{metric}_metadata"].corr(verification[f"{metric}_raw_counts"])
        messages.append(
            f"  {metric}: raw/metadata ratio median {ratio.median():.3f} "
            f"(range {ratio.min():.3f}-{ratio.max():.3f}), Pearson r = {corr:.4f}, "
            f"{within:.0%} of {len(ratio)} nuclei within {tolerance:.0%}")
    return messages


COLUMN_DICTIONARY = [
    ("nucleus_id", "Unique nucleus identifier"),
    ("donor_id", "Donor identifier"),
    ("specimen_id", "Specimen (library) the nucleus comes from"),
    ("Class / Subclass / Supertype", "Cell type, from coarse to fine"),
    ("Neurotypical reference", "True for the 5 neurotypical reference donors"),
    ("ADNC / Cognitive status / Braak stage / ...", "Donor clinical and neuropathological annotations"),
    ("n_umis", "SEQUENCING DEPTH: total UMIs (unique RNA molecules) in the nucleus"),
    ("n_genes", "Number of genes with at least one UMI in the nucleus"),
    ("frac_mito", "Fraction of UMIs from mitochondrial genes (0 to 1)"),
    ("log10_umis / log10_genes", "log10 of n_umis / n_genes, for plotting"),
    ("genes_per_umi", "n_genes / n_umis: library complexity of the nucleus"),
    ("n_nuclei", "(summaries) Number of nuclei in the group"),
    ("n_specimens", "(summaries) Number of distinct specimens in the group"),
    ("umis_* / genes_*", "(summaries) mean, sd, min, q25, median, q75 and max"),
    ("frac_mito_median", "(summaries) Median mitochondrial fraction in the group"),
]


def main():
    parser = argparse.ArgumentParser(description="Sequencing depth by nucleus and donor")
    parser.add_argument("--verify", type=int, default=None,
                        help="nuclei to check against raw/X (overrides config; 0 = no check)")
    args = parser.parse_args()

    config = yaml.safe_load(CONFIG.read_text())["sequencing_depth"]
    metadata_path = Path(config["metadata"])
    out_dir = Path(config["output_dir"])
    n_verify = args.verify if args.verify is not None else config.get("verify_nuclei", 0)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "config_used.yaml").write_text(yaml.safe_dump(config, sort_keys=False))

    print(f"Reading metadata from {metadata_path}...")
    depth = depth_by_nucleus(read_metadata(metadata_path))

    print("Summarizing by donor, donor x subclass and donor x supertype...")
    by_donor = summarize(depth, ["donor_id"])
    by_subclass = summarize(depth, ["donor_id", "Subclass"])
    by_supertype = summarize(depth, ["donor_id", "Supertype"])

    checks = quality_checks(depth)
    for message in checks:
        print("  [check]", message)

    depth.to_parquet(out_dir / "depth_by_nucleus.parquet")
    by_donor.to_csv(out_dir / "depth_by_donor.csv", index=False)
    by_subclass.to_csv(out_dir / "depth_by_donor_subclass.csv", index=False)
    by_supertype.to_csv(out_dir / "depth_by_donor_supertype.csv", index=False)
    pd.DataFrame(COLUMN_DICTIONARY, columns=["column", "description"]).to_csv(
        out_dir / "depth_column_dictionary.csv", index=False)

    info = [
        f"Run: {datetime.now():%Y-%m-%d %H:%M}",
        f"Input: {metadata_path}",
        f"Nuclei: {len(depth):,}",
        f"Donors: {depth['donor_id'].nunique()}",
        f"Specimens: {depth['specimen_id'].nunique()}",
        f"Subclasses: {depth['Subclass'].nunique()}",
        f"Median UMIs per nucleus (all): {depth['n_umis'].median():,.0f}",
        f"Median genes per nucleus (all): {depth['n_genes'].median():,.0f}",
        "Checks:", *[f"  {m}" for m in checks],
    ]

    if n_verify > 0:
        print(f"Checking {n_verify} nuclei against raw/X...")
        verification = verify_with_raw_counts(depth, n_verify, config["anndata_url"])
        info += ["Verification against raw/X:", *summarize_verification(verification)]
        verification.to_csv(out_dir / "verification_raw_counts.csv", index=False)

    (out_dir / "run_info.txt").write_text("\n".join(info) + "\n")
    print("\n".join(info))
    print(f"\nDone. Results in {out_dir}")


if __name__ == "__main__":
    main()
