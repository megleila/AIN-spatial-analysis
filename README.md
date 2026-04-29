# AIN-spatial-analysis

Documentation of the analysis workflow for **Baker, Kakade et al.**, *Spatial analysis reveals the cellular microenvironments and mechanisms of inflammation and kidney injury in acute interstitial nephritis*, *Nature Communications* (in press).

## Purpose

These scripts document the Xenium spatial transcriptomic workflow used in the manuscript: data loading, quality control, normalization, integration, clustering, cell-type annotation, and proximal tubule subclustering. They are written as **annotated documentation** of the analytical steps — not as a self-contained executable pipeline. They mirror the procedures described in the manuscript Methods and reflect the actual workflow used by the authors. File paths and per-tissue specifics are placeholders that would be replaced when applying the workflow to a new dataset; cluster-to-cell-type mappings depend on the stochastic clustering output of any given run and so are documented in Supplementary Tables 6 and 7 rather than hardcoded here.

## Contents

| Script | Description |
|---|---|
| `scripts/01_xenium_load_qc_merge.R` | Loads Xenium acquisitions with `LoadXenium`, attaches per-tissue patient metadata exported from Xenium Explorer, applies barcode prefixing for slide uniqueness, and merges into a single Seurat v5 object. |
| `scripts/02_xenium_normalize_integrate_cluster.R` | QC filter (`nFeature_Xenium ≥ 10`), patient-level layer splitting, SCTransform, PCA, Harmony integration via `IntegrateLayers`, UMAP, multi-resolution Louvain clustering with silhouette evaluation (resolution 0.8 selected per Methods), `FindAllMarkers` for cell-type annotation. |
| `scripts/03_xenium_pt_subclustering.R` | Targeted PT subclustering (21 subclusters → 3 PT states), `HAVCR1`/`VCAM1` stratification of AIN cortex PT cells (4 marker groups), and cross-platform validation against the KPMP snRNA-seq atlas (Lake et al. *Nature* 2023). |

## Software environment

- R 4.3.0
- Seurat 5.0
- Harmony 0.1.1
- batchelor 1.12.3
- imcRtools 1.14.0 (used in the IMC analyses, not in this repository)
- anndata + reticulate (for KPMP `.h5ad`)
- Python 3.9 with scanpy 1.9.3 (for the same `.h5ad` interface)

## Data availability

The processed Seurat object (post-clustering, post-annotation, with cell-type labels and patient metadata) used to produce the manuscript's analyses and figures is deposited at Zenodo: `10.5281/zenodo.YYYYYYY`. Raw Xenium output is deposited at GEO: `GSE…`. The KPMP atlas used for cross-platform validation is publicly available at the [CZ CELLxGENE Discover collection](https://cellxgene.cziscience.com/collections/bcb61471-2a44-4d00-a0af-ff085512674c).

Analysis code is available at GitHub: [github.com/megleila/AIN-spatial-analysis](https://github.com/megleila/AIN-spatial-analysis)

The companion IMC analyses are documented separately and rely on `imcRtools`, the `Rphenograph` clustering implementation, and the SORBET niche-classification framework ([github.com/KlugerLab/SORBET](https://github.com/KlugerLab/SORBET)).

## Citation

If you use these analyses or the deposited data, please cite:

> Baker M.L.\*, Kakade V.R.\*, Budiman T., Weiß M., Cunningham J.M., Sadarangani S., Lerner G., Moeckel G., Rosenberg A.Z., Parikh C.R., Kluger Y., Moledina D.G., Cantley L.G. Spatial analysis reveals the cellular microenvironments and mechanisms of inflammation and kidney injury in acute interstitial nephritis. *Nature Communications* (in press).

(\*co-first authors)

## Contact

Megan L. Baker (megan.baker@yale.edu) and Lloyd G. Cantley.

Supported by NIH grants R01DK126815, R01DK128087, T32DK007276, U01DK133768, P30DK045735, U54DK137331, and ADA Postdoctoral Fellowship 11-23-PDF-63. Specimens were obtained under Johns Hopkins IRB protocols IRB00221958 and IRB00090103.
