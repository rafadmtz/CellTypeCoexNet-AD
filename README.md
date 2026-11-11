# CellTypeCoexNet-AD

> **Status:** work in progress. Preliminary results are available for Vip interneurons.

How does gene coexpression differ across levels of Alzheimer's disease neuropathologic change (ADNC) within a single neuronal population? This project uses single-nucleus RNA sequencing (snRNA-seq) from the Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD) to compare cell-type-specific networks in the human middle temporal gyrus (MTG). The first analysis focuses on Vip interneurons, a subclass that includes subtypes depleted in late-stage pathology (Gabitto et al., 2024).

The preliminary Vip analysis was developed in the [NeuroNet_AD](https://github.com/silvanaycristo/NeuroNet_AD) repository. This repository continues that work with a donor-balanced resampling design and extends it to other cell populations.

## Vip cohort and network properties

The analysis includes 75,942 Vip nuclei across four donor ADNC groups. Donors labeled Reference were excluded. Each row summarizes one network inferred from the Vip nuclei in that group. A connected gene has at least one edge above the group's mutual information (MI) threshold.

| ADNC group | Vip nuclei | MI threshold | Connected genes | Edges | Infomap communities | Modularity | Clustering coefficient |
|---|---:|---:|---:|---:|---:|---:|---:|
| Not AD | 10,711 | 0.209 | 1,005 | 10,016 | 64 | 0.530 | 0.404 |
| Low | 12,218 | 0.216 | 989 | 9,938 | 59 | 0.513 | 0.436 |
| Intermediate | 19,633 | 0.240 | 838 | 8,371 | 51 | 0.500 | 0.483 |
| High | 33,380 | 0.277 | 773 | 7,740 | 50 | 0.413 | 0.560 |

## Analysis

1. Select Vip nuclei from the SEA-AD whole-taxonomy MTG dataset (CELLxGENE), split them by donor ADNC group, and retain protein-coding genes. Genes with zero counts in a group are removed before MI estimation.
2. Estimate MI between gene pairs within each group with `numba_mi`, an MI estimator developed in the de Anda-Jáuregui lab. Rank gene pairs by MI and choose a group-specific cutoff that gives a connected-gene-to-edge ratio close to 0.1, corresponding to a mean degree of about 20 among connected genes (de Anda-Jáuregui et al., 2016). Genes with no edge above the cutoff are assigned degree 0.
3. Identify communities with Infomap (1,000 trials). Infomap was chosen among seven algorithms because its partitions were the most stable across random seeds (AMI 0.93–0.96). Test each community for over-representation of Gene Ontology Biological Process terms with clusterProfiler (Benjamini–Hochberg adjusted *p* < 0.05; universe = all genes measured in the MI matrix). Compare groups using Jaccard overlap of retained edges, Jaccard overlap of enriched terms, and module-to-module similarity (de Anda-Jáuregui, 2019). Cross-group comparisons use unsimplified GO results; simplified terms are used only for figures.
4. Track each gene's degree, its connectivity within its community (within-module degree z-score), and how its edges are distributed across communities (participation coefficient) (Guimerà & Amaral, 2005).

## Preliminary findings

- **Network structure.** High had the fewest connected genes and communities, the lowest modularity, and the highest clustering coefficient. The MI threshold rose with the number of nuclei, from 0.209 in Not AD to 0.277 in High. Not AD and High shared 5.8% of the distinct edges found in either network (Jaccard index 0.058). Module similarity to Not AD decreased with ADNC level (Low 0.147, Intermediate 0.119, High 0.077).
- **Community function.** Oxidative phosphorylation appeared in all four groups in a small community of 9–11 genes. RNA processing terms also appeared in all four groups, in progressively larger communities. In High, mRNA regulation and mRNA destabilization terms were concentrated in one community of 178 genes. Synaptic vesicle recycling no longer formed its own community and appeared in a 184-gene community together with nucleocytoplasmic transport. Protein quality control terms were enriched in another High community of 98 genes. These patterns come from simplified GO terms and will be checked against the unsimplified results.
- **Splicing-related genes.** From Not AD to High, the number of retained neighbors rose from 17 to 201 for SUGP2 and from 8 to 119 for LUC7L. These network changes motivate direct tests of whether splice-site usage differs across ADNC groups.

## Limitations

- Each ADNC group is represented by a single network that pools nuclei from all donors, so donor effects are not separated from group effects.
- The number of Vip nuclei differs widely between groups (High has about three times as many as Not AD), and the MI threshold rises with it.
- Networks with the same connected-gene-to-edge ratio but fewer genes are denser, so part of the higher clustering in High may reflect network size.
- The analysis is descriptive; differences between networks have not yet been tested statistically.

## Validation and next steps

To be eligible for resampling, a donor must contribute at least 100 Vip nuclei. Eligibility was calculated from the whole-taxonomy MTG dataset by counting Vip nuclei for each donor within each ADNC group. The 75,942 nuclei in the network table above are those retained in the four analysis matrices; reproducing that exact subset from the public file requires the list of retained nuclei.

| ADNC group | Donors in cohort | Donors with ≥100 Vip nuclei |
|---|---:|---:|
| Not AD | 9 | 9 |
| Low | 12 | 12 |
| Intermediate | 21 | 21 |
| High | 42 | 40 |

The two High donors below the threshold contributed 93 and 15 Vip nuclei. Each also contributed fewer nuclei than the median of the 40 eligible High donors in 21 of the 23 other subclasses, which suggests low nuclei yield for these donors rather than a selective loss of Vip nuclei.

As a first approach, we will run 10 resampling iterations. In each iteration, we will select 6 eligible donors per ADNC group and sample 100 Vip nuclei from each donor, taking into account differences in sequencing depth distributions between groups. We will infer one network per group from 600 nuclei, for a total of 2,400 nuclei across the four networks in each iteration. Because Not AD has only 9 eligible donors, only 84 donor combinations are possible in that group, so its resamples will overlap more than those of other groups. We will assess whether edge overlap, community structure, and the degree changes in SUGP2 and LUC7L persist across resamples.

We then plan to apply the same workflow to Sst interneurons and layer 2/3 intratelencephalic (L2/3 IT) excitatory neurons. These choices follow reports of selective vulnerability among Sst and superficial-layer excitatory subtypes, as described by Gabitto et al. (2024) and reviewed by Emerson et al. (2026). This comparison will test which patterns recur across populations and which are specific to Vip interneurons.

## Repository structure

```
CellTypeCoexNet-AD/
├── config/      # Analysis configuration
├── docs/        # Documentation
├── envs/        # Conda environments (Python and R)
├── notebooks/   # Exploratory notebooks
├── resources/   # Auxiliary resources
├── results/Vip/ # Vip results
└── src/
    ├── preprocessing/  # Data retrieval and nuclei counts
    └── utils/          # Analysis functions
```

## Data and code

| Resource | Purpose |
|---|---|
| SEA-AD MTG expression matrix (.h5ad) | Whole-taxonomy snRNA-seq input from [CELLxGENE](https://cellxgene.cziscience.com/), including donor and cell-type annotations. Dataset ID: `[TODO]`. |
| [Data retrieval and preprocessing](src/preprocessing/) | Scripts for retrieving the MTG dataset and counting nuclei by donor, subclass, and supertype. |
| [Sequencing depth](docs/sequencing_depth.md) | Depth (UMIs), detected genes and mitochondrial fraction for every nucleus, summarized by donor, subclass and supertype. |
| [Analysis code](src/utils/) | Python and R functions for mutual information estimation, network topology, gene-level changes, community detection, and functional enrichment. |

Expression data are not redistributed in this repository. When using them, cite Gabitto et al. (2024).

## How to cite

For the preliminary Vip results, cite the poster using its original Spanish title.

Díaz-Martínez, R., Cristo-Martínez, S. Y., Miyazono-Ushijima, A. K., & de Anda-Jáuregui, G. (2026, November 12). *Análisis topológico y funcional de redes de coexpresión génica en interneuronas Vip entre estadios neuropatológicos de Alzheimer* [Conference poster]. 12.º Encuentro Estudiantil «Del gen al genoma: 100 años de genética», Instituto Nacional de Medicina Genómica. [Event page](TODO)

When reusing the analysis code, cite the repository and include the commit used. The code is available under the MIT License.

rafadmtz. (2026). *CellTypeCoexNet-AD* [Source code]. [GitHub](https://github.com/rafadmtz/CellTypeCoexNet-AD)

## References

- de Anda-Jáuregui, G. (2019). Guideline for comparing functional enrichment of biological network modular structures. *Applied Network Science*. [DOI](https://doi.org/10.1007/s41109-019-0129-0)
- de Anda-Jáuregui, G., Velázquez-Caldelas, T. E., Espinal-Enríquez, J., & Hernández-Lemus, E. (2016). Transcriptional network architecture of breast cancer molecular subtypes. *Frontiers in Physiology*. [DOI](https://doi.org/10.3389/fphys.2016.00568)
- Emerson, N., Vallejo, K., Rahimzadeh, N., Kwang, N. E., & Swarup, V. (2026). Single-cell and spatial transcriptomics of cellular heterogeneity in Alzheimer disease. *Nature Reviews Neuroscience*. [DOI](TODO)
- Gabitto, M. I., Travaglini, K. J., Rachleff, V. M., et al. (2024). Integrated multimodal cell atlas of Alzheimer's disease. *Nature Neuroscience*. [DOI](https://doi.org/10.1038/s41593-024-01774-5)
- Guimerà, R., & Amaral, L. A. N. (2005). Functional cartography of complex metabolic networks. *Nature*. [DOI](https://doi.org/10.1038/nature03288)
- Rosvall, M., & Bergstrom, C. T. (2008). Maps of random walks on complex networks reveal community structure. *PNAS*. [DOI](https://doi.org/10.1073/pnas.0706851105)
- Wu, T., Hu, E., Xu, S., et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation*. [DOI](https://doi.org/10.1016/j.xinn.2021.100141)
