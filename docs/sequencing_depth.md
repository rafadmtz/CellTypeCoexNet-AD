# Sequencing depth

Sequencing depth is the total number of UMIs (unique RNA molecules) detected in
a nucleus. Nuclei with lower depth have more zero counts and fewer detected
genes, which affects mutual information estimates. Depth distributions are
therefore compared across donors, cell types and ADNC groups before sampling
nuclei for network inference.

## Inputs

The cell metadata of the SEA-AD MTG dataset (CELLxGENE), one row per nucleus.
It already contains `Number of UMIs`, `Genes detected` and
`Fraction mitochrondrial UMIs`, so the expression matrix is not needed.

These values were computed by SEA-AD on its full gene set. The CELLxGENE
expression matrix keeps 35,483 genes, so summing a nucleus's row of `raw/X`
gives slightly lower values. A check of 20 random nuclei gave raw/metadata
ratios of 0.992–0.997 for UMIs and 0.989–0.996 for detected genes, with a
correlation of 1.0, so the metadata and the matrix rows are aligned. The
depth reported here is the SEA-AD value.

The metadata can be obtained in two ways:

1. `src/preprocessing/get_MTG_metadata_streaming.py` reads only the metadata
   from the CELLxGENE URL, without downloading the expression matrix.
2. `src/preprocessing/get_MTG_snRNAseq_anndata.py` downloads the full h5ad and
   exports the metadata as CSV (Reference donors excluded).

The path to the metadata is set in `config/config.yaml` (`sequencing_depth.metadata`).

## Running

From the repository root:

```bash
python src/preprocessing/get_MTG_metadata_streaming.py   # only if data/raw has no metadata
python src/preprocessing/sequencing_depth.py
python src/preprocessing/sequencing_depth.py --verify 20 # optional: check 20 nuclei against raw/X
```

In Google Colab, with the metadata already in Google Drive:

```python
from google.colab import drive
drive.mount('/content/drive')

!git clone https://github.com/rafadmtz/CellTypeCoexNet-AD.git
%cd CellTypeCoexNet-AD
!pip install -q "pandas==2.2.3" "anndata==0.12.19" h5py fsspec aiohttp pyarrow pyyaml

# Link the metadata stored in Drive to the path expected by config/config.yaml
!mkdir -p data/raw
!ln -s "/content/drive/MyDrive/SEA-AD/seaad_mtg_obs.parquet" data/raw/MTG_snRNAseq_metadata.parquet

!python src/preprocessing/sequencing_depth.py
```

## Outputs

Written to `results/sequencing_depth/tables/`.

| File | Content |
|---|---|
| `depth_by_nucleus.parquet` | One row per nucleus: donor, specimen, cell type, donor annotations and depth metrics. Not committed (too large); shared separately. |
| `depth_by_donor.csv` | One row per donor: number of nuclei and mean, SD, min, quartiles and max of UMIs and detected genes |
| `depth_by_donor_subclass.csv` | Same, for each donor and subclass |
| `depth_by_donor_supertype.csv` | Same, for each donor and supertype |
| `depth_column_dictionary.csv` | Description of every column |
| `run_info.txt` | Date, input, counts and quality checks |
| `config_used.yaml` | Configuration used in the run |
| `verification_raw_counts.csv` | Only with `--verify`: metadata depth vs. raw counts, with raw/metadata ratios |

All nuclei are kept, including Reference donors. Filter with `ADNC` or
`Neurotypical reference` when plotting.

Read the nucleus-level table in Python with `pd.read_parquet()` and in R with
`arrow::read_parquet()`.
