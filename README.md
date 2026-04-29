# AIN-spatial-analysis

Documentation of the analysis workflow for **Baker, Kakade et al.**, *Spatial analysis reveals the cellular microenvironments and mechanisms of inflammation and kidney injury in acute interstitial nephritis*. Manuscript under revision at *Nature Communications*.

## Purpose

These scripts document the two spatial profiling workflows used in the manuscript:

- **Imaging mass cytometry (IMC)** with a 31-marker antibody panel applied to 106 kidney tissues across the discovery (Yale, n=65) and validation (Johns Hopkins, n=41) cohorts. Used for protein-level spatial cell-type mapping, PT subtype identification (PT-H, PT-INJ1, PT-INJ2), and neighborhood analyses.
- **Single-cell spatial transcriptomics (ST)** with the 10x Genomics Xenium Prime 5K Human Pan Tissue & Pathways Panel plus a 100-gene custom panel, applied to 20 tissues from the discovery cohort. Used for transcript-level cell typing, PT injury state cross-validation against the KPMP atlas, ligand-receptor interaction analysis, and lymphoid aggregate identification.

The scripts are written as **annotated documentation** of the analytical steps, not as a self-contained executable pipeline. They mirror the procedures described in the manuscript Methods. File paths and per-tissue specifics are placeholders; cluster-to-cell-type mappings depend on the stochastic clustering output of any given run and so are documented in Supplementary Tables 4, 6, and 7 rather than hardcoded here.

Raw data ingestion steps (steinbock pre-processing of `.mcd` files for IMC; Xenium Analyzer cell segmentation) are described in the manuscript Methods but not included here, as they operate on raw imaging files rather than on the cell-level data structures used for the downstream analyses.

## Contents

### IMC scripts (`imc/`)

| Script | Description |
|---|---|
| `imc/01_imc_normalize_cluster_annotate.R` | Spillover correction (Chevrier et al. workflow), arcsinh transformation (cofactor 1), `fastMNN` batch correction across slides via `batchelor`, Rphenograph clustering on the corrected embedding, targeted subclustering using cluster-enriched markers, manual cell-type annotation. |
| `imc/02_imc_neighborhoods_density.R` | k-nearest-neighbor spatial graph construction (`k=5`, calibrated to capture immediately adjacent neighbors), per-patient cell densities (cells/mm²) using manually measured tissue areas, Kruskal-Wallis comparisons across AIN/ATI/reference, generation of index-neighbor pair tables for GLMM analyses (the GLMMs themselves were fit in Stata 18). |

### Xenium spatial transcriptomics scripts (`xenium/`)

| Script | Description |
|---|---|
| `xenium/01_xenium_load_qc_merge.R` | Loads Xenium acquisitions with `LoadXenium`, attaches per-tissue patient metadata, applies barcode prefixing for slide uniqueness, and merges into a single Seurat v5 object. |
| `xenium/02_xenium_normalize_integrate_cluster.R` | QC filter (`nFeature_Xenium >= 10`), patient-level layer splitting, SCTransform, PCA, Harmony integration via `IntegrateLayers`, UMAP, multi-resolution Louvain clustering with silhouette evaluation (resolution 0.8 selected per Methods), `FindAllMarkers` for cell-type annotation. |
| `xenium/03_xenium_pt_subclustering.R` | Targeted PT subclustering (21 subclusters → 3 PT states: PT-H, PT-dd, PT-INJ), `HAVCR1`/`VCAM1` stratification of AIN cortex PT cells (4 marker groups), and cross-platform validation against the KPMP snRNA-seq atlas (Lake et al. *Nature* 2023). |

## Software environment

- R 4.3.0
- Seurat 5.0
- Harmony 0.1.1
- batchelor 1.12.3
- imcRtools 1.6.0
- CATALYST 1.24.0 (signal spillover correction)
- Rphenograph (`JinmiaoChenLab/Rphenograph` on GitHub)
- SpatialExperiment 1.10.0
- anndata + reticulate (for KPMP `.h5ad`)
- Python 3.9 with scanpy 1.9.3 (for the same `.h5ad` interface)
- Stata 18 (for IMC GLMM neighbor-pair statistics, per Methods)
- Fiji/ImageJ (for IMC tissue area measurements)

## Data availability

Processed data objects used to produce the manuscript's analyses and figures are deposited at Zenodo: `10.5281/zenodo.YYYYYYY`. The deposit contains:

- `xenium_AIN_for_deposit.rds` — Seurat v5 object with 321,333 cells from 20 discovery cohort tissues, with cell-type annotations, Harmony-corrected embeddings, and per-patient clinical metadata.
- `imc_yale_discovery_for_deposit.rds` — `SpatialExperiment` object with 1,331,664 cells from the 65-patient discovery cohort.
- `imc_jhu_validation_for_deposit.rds` — `SpatialExperiment` object with 1,216,118 cells from the 41-patient validation cohort.
- `xenium_clinical_metadata_for_deposit.csv` — per-patient clinical fields for the Xenium cohort.
- `imc_clinical_metadata_for_deposit.csv` — per-patient clinical fields for the combined IMC cohort.

Patient identifiers are de-identified throughout: `AIN_##_D` / `ATI_##_D` / `REF_##_D` for discovery cohort patients; `AIN_##_V` / `ATI_##_V` / `REF_##_V` for validation cohort patients. The 20 patients shared between the Xenium and IMC discovery cohorts use matching identifiers across both modalities. Where multiple tissue sections were profiled from the same biological reference patient in the validation cohort, those sections share a single `REF_##_V` patient identifier with section identity preserved via `sample_id`.

Cell-type label columns in the IMC objects:

- `simplified_annotations` and `core_annotations` are the internal working labels used during clustering and analysis (preserved for backward compatibility with previously written analysis scripts).
- `manuscript_label` and `manuscript_core_label` give the published cluster names from Supplementary Table 4 (e.g., `TAL-H` rather than the internal `TAL-Umod-hi`; `MP-CD206+` rather than `Macs-CD206+`). The IMC scripts in this repository use `manuscript_label` so analyses align with the names readers see in the manuscript figures.

In the Xenium object, T cell subtype columns (`final_cd4`, `final_cd8`) reflect the published nomenclature in Extended Data Fig. 3c,d (e.g., `CD8_Exhausted`, `CD8_Checkpoint_high_TRM`, `CD4_naive`, `CD4_Memory`, `CD4_Treg`).

Raw Xenium output will be deposited at GEO upon manuscript acceptance. Raw IMC `.mcd` files are available from the corresponding authors on reasonable request.

The KPMP atlas used for cross-platform PT injury state validation is publicly available at the [CZ CELLxGENE Discover collection](https://cellxgene.cziscience.com/collections/bcb61471-2a44-4d00-a0af-ff085512674c).

## Code availability

This GitHub repository is archived on Zenodo at `10.5281/zenodo.ZZZZZZZ`.

Companion methods rely on:
- The Bodenmiller [IMC Workflow](https://bodenmillergroup.github.io/IMCWorkflow/) for IMC pre-processing and analysis
- [steinbock](https://github.com/BodenmillerGroup/steinbock) for IMC raw image pre-processing
- [Mesmer](https://github.com/vanvalenlab/intro-to-deepcell) for IMC cell segmentation, with custom kidney-tuned parameters (resolution 0.8, small object threshold 53, maxima threshold 0.3, interior threshold 0.47, 1µm pixel expansion)
- [SORBET](https://github.com/KlugerLab/SORBET) for ST PT injury microenvironment classification

## Citation

A formal citation will be added here once the manuscript is accepted for publication. Until then, please contact the corresponding authors before using these analyses.

## Contact

Megan L. Baker (megan.baker@yale.edu) and Lloyd G. Cantley.

Supported by NIH grants R01DK126815, R01DK128087, T32DK007276, U01DK133768, P30DK045735, U54DK137331, and ADA Postdoctoral Fellowship 11-23-PDF-63. Specimens were obtained under Yale IRB protocol 11110009286 and Johns Hopkins IRB protocols IRB00221958 and IRB00090103.
