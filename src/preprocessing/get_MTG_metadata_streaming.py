"""
Extract the cell metadata (obs) and gene annotation (var) of the SEA-AD MTG
snRNA-seq dataset from CELLxGENE without downloading the expression matrix.

The h5ad file is read over HTTP, and only the obs and var groups are fetched.
This is an alternative to get_MTG_snRNAseq_anndata.py when the full file
cannot be downloaded (for example, on a laptop or in Colab).

Unlike get_MTG_snRNAseq_anndata.py, all nuclei are kept, including those from
Reference donors. Use the "Neurotypical reference" or "ADNC" column to
exclude them downstream.

Outputs (paths set in config/config.yaml, section sequencing_depth):
    data/raw/MTG_snRNAseq_metadata.parquet       one row per nucleus
    data/raw/MTG_snRNAseq_var.parquet            one row per gene
    data/raw/MTG_snRNAseq_metadata_info.txt      source URL, date and shape

Usage (from the repository root):
    python src/preprocessing/get_MTG_metadata_streaming.py
"""

from datetime import datetime
from pathlib import Path

import fsspec
import h5py
import yaml

try:
    from anndata.io import read_elem
except ImportError:  # older anndata versions
    from anndata.experimental import read_elem

CONFIG = Path("config/config.yaml")


def get_metadata_streaming(url: str, obs_path: Path, var_path: Path) -> Path:
    """
    Read the obs and var groups of a remote h5ad file and save them as parquet.

    Parameters:

    url : str
        URL (or local path) of the h5ad file.
    obs_path : Path
        Output path for the cell metadata.
    var_path : Path
        Output path for the gene annotation.

    Returns:

    Path
        Path to the cell metadata parquet file.
    """
    obs_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Reading metadata from {url}...")
    with fsspec.open(url, "rb", block_size=16 * 2**20, cache_type="blockcache") as f:
        with h5py.File(f, "r") as h5:
            obs = read_elem(h5["obs"])
            var = read_elem(h5["var"])

    obs.to_parquet(obs_path)
    var.to_parquet(var_path)

    info = [
        f"Extracted: {datetime.now():%Y-%m-%d %H:%M}",
        f"Source: {url}",
        f"obs: {obs.shape[0]:,} nuclei x {obs.shape[1]} columns -> {obs_path}",
        f"var: {var.shape[0]:,} genes x {var.shape[1]} columns -> {var_path}",
        "obs columns:", *[f"  - {c}" for c in obs.columns],
    ]
    info_path = obs_path.with_name(obs_path.stem + "_info.txt")
    info_path.write_text("\n".join(info) + "\n")
    print("\n".join(info))
    return obs_path


def main():
    config = yaml.safe_load(CONFIG.read_text())["sequencing_depth"]
    get_metadata_streaming(
        config["anndata_url"],
        Path(config["metadata"]),
        Path(config["gene_annotation"]),
    )


if __name__ == "__main__":
    main()
