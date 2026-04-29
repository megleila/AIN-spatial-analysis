# =============================================================================
# 02_imc_neighborhoods_density.R
# -----------------------------------------------------------------------------
# Code for: Baker, Kakade et al. Manuscript under revision at Nature
# Communications.
#
# This script documents the spatial-neighborhood and cell-density analyses
# applied to the annotated, batch-corrected SpatialExperiment from script 01.
# It mirrors the procedures described in the Methods sections "Spatial
# neighborhood analysis" and "Statistical analysis of cell densities and
# spatial neighbors for IMC data," and corresponds to Fig. 2d-g, Fig. 3d,
# Fig. 5a, Extended Data Figs. 4-5, and the GLMM-based PT subtype neighbor
# enrichments shown in Fig. 2f-g.
#
# Inputs:
#   - Annotated SpatialExperiment from script 01 (imc_annotated_batch_corrected.rds)
#   - Per-ROI ablated tissue area data (Excel) from manual ImageJ measurements
#
# Outputs:
#   - Spatial kNN interaction graph stored as a colPair on the object
#   - Per-patient cell density tables (cells/mm^2) used for Kruskal-Wallis
#     comparisons of cell type abundance across AIN, ATI, reference
#   - Per-patient index-neighbor pair counts used for GLMM neighbor analyses
#     (the GLMM models themselves were fit in Stata 18 per the Methods)
# =============================================================================

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(imcRtools)
    library(readxl)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(stats)
})

# -----------------------------------------------------------------------------
# 1. Load annotated, batch-corrected SpatialExperiment
# -----------------------------------------------------------------------------
spe <- readRDS("imc_annotated_batch_corrected.rds")

# -----------------------------------------------------------------------------
# 2. Build kNN spatial interaction graph (k = 5)
# -----------------------------------------------------------------------------
# Per the Methods: "Using 2 living donor biopsy tissues, 1 AKI biopsy tissue,
# and 1 CKD biopsy tissue, we manually visually counted the number of
# immediately adjacent neighboring cells for a subset of 100 random cells
# per tissue. Across 400 cells and 3 tissue types, we found that individual
# cells in our kidney tissues have 5 immediately adjacent neighboring cells,
# on average, with a mean centroid-to-centroid distance of 16 µm. Using a
# k-nearest neighbors approach and a threshold of 5, we generated a spatial
# interaction graph which was used to determine the identity of each cell's
# immediately adjacent neighbors across the entire project."

spe <- buildSpatialGraph(spe,
                         img_id = "sample_id",   # graph computed within ROIs
                         type   = "knn",
                         k      = 5)
# The interaction graph is stored as a colPair (`knn_interaction_graph`).
# Each row of the colPair is one (from cell, to cell) edge.

# -----------------------------------------------------------------------------
# 3. Per-patient cell densities (cells/mm^2)
# -----------------------------------------------------------------------------
# Per the Methods: "Cell type densities (cells/mm²) in IMC data were
# calculated per patient by dividing total cell count by tissue area,
# yielding one density value per cell type per patient. Tissue area was
# measured manually using Fiji/ImageJ."
#
# Tissue areas live in an external Excel file (one row per ROI, with patient
# and category labels). Areas are summed per patient before computing
# densities to give one value per patient per cell type.

areas <- read_excel("path/to/Ablated_Area_per_ROI_Data.xlsx")
# Expected columns: Sample_id, Patient, Indication, Size (in mm^2)

cd <- as.data.frame(colData(spe))

# Sum tissue area per patient
patient_area <- areas %>%
    group_by(Patient) %>%
    summarise(area_mm2 = sum(Size), .groups = "drop")

# Cell counts per patient per published cell type (manuscript_label)
cell_counts <- cd %>%
    group_by(Patient, manuscript_label) %>%
    summarise(n_cells = n(), .groups = "drop")

# Densities (cells/mm^2)
densities <- cell_counts %>%
    left_join(patient_area, by = "Patient") %>%
    mutate(density = n_cells / area_mm2)

# -----------------------------------------------------------------------------
# 4. Density comparisons across diagnostic categories
# -----------------------------------------------------------------------------
# Per the Methods: "Densities were compared between AIN, ATI, and reference
# tissue groups using by Kruskal-Wallis test with post-hoc pairwise Wilcoxon
# rank-sum tests with Bonferroni correction for multiple comparisons."
#
# Densities are compared separately within each tissue region (cortex /
# medulla) -- e.g., Fig. 2d compares PT subtype densities across AIN, ATI,
# and reference cortex tissues only. Patients contribute one density per
# cell type per region.

# Attach patient_category and tissue_region from colData (one row per patient
# per region, with patient_category derived from indication)
patient_meta <- cd %>%
    distinct(Patient, indication) %>%
    mutate(patient_category = case_when(
                indication %in% c("AC","AM") ~ "AIN",
                indication %in% c("BC","BM") ~ "ATI",
                indication %in% c("RC","RM") ~ "REF"),
           tissue_region    = ifelse(indication %in% c("AC","BC","RC"),
                                     "Cortex","Medulla"))

dens_with_meta <- densities %>%
    left_join(patient_meta, by = "Patient")

# For one cell type and one region, run the Kruskal-Wallis omnibus test
# followed by pairwise Wilcoxon with Bonferroni correction.
test_one_celltype <- function(df_one_celltype) {
    kw <- kruskal.test(density ~ patient_category, data = df_one_celltype)
    pairs <- pairwise.wilcox.test(df_one_celltype$density,
                                  df_one_celltype$patient_category,
                                  p.adjust.method = "bonferroni")
    list(kw_p = kw$p.value, pairwise_p = pairs$p.value)
}

# Apply across all published cell types within each region:
density_stats <- dens_with_meta %>%
    filter(tissue_region == "Cortex") %>%
    group_by(manuscript_label) %>%
    do(test = test_one_celltype(.))

# -----------------------------------------------------------------------------
# 5. Index-neighbor pair tables for GLMM neighbor analyses
# -----------------------------------------------------------------------------
# Per the Methods: "For analysis testing the association of diagnosis with
# number of index and surrounding cell pairs, we allowed each type of cell
# to be surrounded by up to 5 nearest neighbors. To test the association of
# each diagnosis with index cell-surrounding cell pair, we fit a generalized
# linear mixed model with outcome as number of index cell surrounding cell
# pairs and predictor as diagnosis using log link. We clustered the analysis
# at the level of participant. We present fold difference in number of cells
# between various diagnoses (β coefficients) as well as multiple comparison
# adjusted q-values using Benjamini Hochberg procedure (Fig. 2f; Extended
# Data Fig. 5a)."
#
# Per the Methods, the GLMMs themselves were fit in Stata 18. Here we
# generate the input table that was exported from R and read into Stata.
# Each row is one (index cell type, neighbor cell type, patient, diagnosis)
# combination with the number of cell pairs.

# Note on cell type label columns:
# The deposit object includes both the original `simplified_annotations`
# (internal/working labels used during analysis) and `manuscript_label`
# (the published cluster names from Supplementary Table 4). The published
# labels are used here so that index/neighbor types match the names readers
# see in the manuscript figures and tables. The two columns are 1:1 mappable.

knn <- colPair(spe, "knn_interaction_graph")

pair_df <- data.frame(
    from_cell      = knn@from,
    to_cell        = knn@to,
    from_celltype  = colData(spe)$manuscript_label[knn@from],
    to_celltype    = colData(spe)$manuscript_label[knn@to],
    Patient        = colData(spe)$Patient[knn@from],
    indication     = colData(spe)$indication[knn@from]
)

# Pairs per (patient, indication, from_celltype, to_celltype)
pair_counts <- pair_df %>%
    group_by(Patient, indication, from_celltype, to_celltype) %>%
    summarise(n_pairs = n(), .groups = "drop")

# Restrict to tubular index cells with immune / endothelial / vascular /
# stromal neighbors as described in the Methods:
#   "For tubular-interstitial interaction analysis, we used tubular epithelial
#    clusters as index cell types and immune, endothelial, vascular, and
#    stromal clusters as neighboring cell types."
# Cluster names below match Supplementary Table 4 (manuscript_label).
tubular_indices <- c("PT-H","PT-dd","PT-INJ1","PT-INJ2","PT-Ki67+",
                     "TAL-H","TAL-dd","TAL-FACL4+","TAL-Ki67+",
                     "DCT-Calb-hi","DCT-Calb-low","DCT-FACL4+","DCT-Ki67+",
                     "CD-CK7-hi","CD-CK7-low","CD-FACL4+","CD-Ki67+",
                     "CNT","tDL","ddTE","ddTE-FACL4+",
                     "DCT")             # validation cohort only
interstitial_neighbors <- c(
    "CD4+ T cells","CD4+ T-IL9+","CD4+ T-Ki67+",
    "CD8+ T cells","CD8+ T-IL9+","CD8+ T-Ki67+",
    "CD4/8-low/- T cells","CD4/8-low/- T-Ki67+",
    "B cells","Mixed B & T cells","Mixed B & T-Ki67+ cells",
    "MP-CD206+","MP-CD206+-Ki67+","MP-CD206-",
    "MP & T cells","MP & T cells-IL9+","MP & T cells-Ki67+","MP & CD8+ T cells",
    "Mast","Eos",
    "ECs","LECs","Vascular",
    "Stromal","Stromal & Vascular"
)

pair_counts_tubular <- pair_counts %>%
    filter(from_celltype %in% tubular_indices,
           to_celltype   %in% interstitial_neighbors)

# Export to .csv for analysis in Stata 18 per the Methods
write.csv(pair_counts_tubular,
          "imc_pair_counts_for_glmm.csv", row.names = FALSE)

# The GLMM was fit in Stata using:
#    meglm n_pairs i.indication, family(poisson) link(log) || Patient:
#    margins indication
# with multiple-comparison correction across cell-type pairs by Benjamini-
# Hochberg.

# -----------------------------------------------------------------------------
# 6. Within-AIN PT subtype neighbor analysis (Fig. 2g, Extended Data Fig. 5b)
# -----------------------------------------------------------------------------
# Per the Methods: "To examine whether different PT injury states within AIN
# tissues exhibited distinct patterns of immune cell recruitment, we analyzed
# the immediate cell neighbors of PT-H, PT-INJ1, and PT-INJ2 cells exclusively
# within AIN cortex samples. For each PT cell, we identified up to 5 nearest
# neighbors. We fit a generalized linear mixed model with log link for each
# surrounding cell type, with the number of index-surrounding cell pairs as
# the outcome and PT cell type as the predictor, using PT-H as the reference
# category."

ain_cortex_pairs <- pair_counts %>%
    filter(indication == "AC",
           from_celltype %in% c("PT-H","PT-INJ1","PT-INJ2"))

write.csv(ain_cortex_pairs,
          "imc_AIN_cortex_PT_neighbor_pairs_for_glmm.csv", row.names = FALSE)

# Fit in Stata with PT-H as the reference category.

# -----------------------------------------------------------------------------
# 7. Save updated SpatialExperiment with kNN graph
# -----------------------------------------------------------------------------
saveRDS(spe, "imc_with_knn_graph.rds")
