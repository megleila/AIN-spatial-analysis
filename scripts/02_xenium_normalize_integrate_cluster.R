# =============================================================================
# 02_xenium_normalize_integrate_cluster.R
# -----------------------------------------------------------------------------
# Code for: Baker, Kakade et al., Nature Communications.
#
# This script documents the normalization, integration, and clustering steps
# applied to the merged AIN Xenium object. It mirrors what is described in
# the Methods section under "Spatial transcriptomic analysis":
#
#   "Single-cell spatial transcriptomics data were processed using Seurat v5.
#    Following quality control to remove cells expressing fewer than 10 unique
#    genes, we created a unified Seurat object containing patient identifiers,
#    region, and category information. [...] Data layers were split by patient
#    groups to enable sample-specific processing before integration. Expression
#    matrices underwent SCTransform normalization [...] To address batch
#    effects while preserving biological variation, we applied Harmony
#    integration through Seurat's IntegrateLayers function [...] Shared
#    nearest neighbor graphs were constructed using the first 30 Harmony
#    dimensions, followed by Louvain clustering. Optimal clustering resolution
#    was determined by evaluating silhouette scores across multiple parameters
#    (0.4-1.2), with resolution 0.8 providing the highest score and optimal
#    biological separation."
#
# Inputs:  the raw merged object from script 01.
# Outputs: a clustered, batch-corrected Seurat v5 object with cell-type
#          annotations in `final_annotations2`.
# =============================================================================

suppressPackageStartupMessages({
    library(Seurat)         # v5.0
    library(harmony)        # v0.1.1
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(cluster)
})

set.seed(42)

# Required for large datasets (321k cells)
options(future.globals.maxSize = 50 * 1024^3)   # 50 GB

# -----------------------------------------------------------------------------
# 1. Load raw merged object
# -----------------------------------------------------------------------------
xenium <- readRDS("xenium_AIN_raw_merged.rds")

# -----------------------------------------------------------------------------
# 2. QC: remove cells with fewer than 10 unique genes
# -----------------------------------------------------------------------------
# Per the Methods. Manually-annotated rare cell types (Eosinophils, Parietal
# cells) identified by H&E morphology are exempted from this filter, as some
# of these cells were identifiable by appearance but had below-threshold
# transcript counts.
keep <- xenium$nFeature_Xenium >= 10 |
        xenium$final_annotations2 %in% c("Eosinophils", "Parietal")
xenium <- xenium[, keep]
cat("Cells after QC:", ncol(xenium), "\n")

# -----------------------------------------------------------------------------
# 3. Split layers by patient
# -----------------------------------------------------------------------------
# Splitting by patient enables per-patient SCTransform normalization, which
# accounts for tissue-quality and probe-efficiency differences between
# individual patient samples. This subsumes any slide-level batch effects
# (each patient is acquired on exactly one slide).
xenium[["Xenium"]] <- JoinLayers(xenium[["Xenium"]])
xenium[["Xenium"]] <- split(xenium[["Xenium"]], f = xenium$base_patient_id)

cat("Number of layers (patients):", length(Layers(xenium)), "\n")

# -----------------------------------------------------------------------------
# 4. SCTransform normalization
# -----------------------------------------------------------------------------
# SCTransform applies a regularized negative binomial regression per layer,
# producing a per-patient normalized count matrix with the top 3000 variable
# genes used for downstream analysis. The active assay becomes `SCT`
# automatically.
xenium <- SCTransform(xenium, assay = "Xenium", verbose = FALSE)

# -----------------------------------------------------------------------------
# 5. PCA
# -----------------------------------------------------------------------------
# 50 PCs computed on the SCT-normalized data.
xenium <- RunPCA(xenium, verbose = FALSE)

# Optional checkpoint
# saveRDS(xenium, "xenium_AIN_sct_pca.rds")

# -----------------------------------------------------------------------------
# 6. Harmony integration
# -----------------------------------------------------------------------------
# IntegrateLayers with method = HarmonyIntegration produces a corrected
# embedding stored as `harmony`. The original SCT expression data is not
# modified -- only the reduced-dimensional space used for clustering and
# visualization changes.
xenium <- IntegrateLayers(
    object         = xenium,
    method         = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction  = "harmony",
    verbose        = FALSE
)

# Re-join the per-patient layers (clustering and DE work on a unified layer)
xenium[["Xenium"]] <- JoinLayers(xenium[["Xenium"]])

# -----------------------------------------------------------------------------
# 7. UMAP embeddings
# -----------------------------------------------------------------------------
# Two UMAPs computed: one from the uncorrected PCA (for diagnostic comparison)
# and one from the Harmony-corrected embedding (used for all visualization).
xenium <- RunUMAP(xenium, reduction = "pca",     dims = 1:30,
                  reduction.name = "umap.pca",     reduction.key = "UMAPPCA_")
xenium <- RunUMAP(xenium, reduction = "harmony", dims = 1:30,
                  reduction.name = "umap.harmony", reduction.key = "UMAPHARMONY_")

# -----------------------------------------------------------------------------
# 8. Clustering at multiple resolutions
# -----------------------------------------------------------------------------
# SNN graph from the first 30 Harmony dimensions, then Louvain clustering at
# five resolutions to support resolution-selection diagnostics.
xenium <- FindNeighbors(xenium, reduction = "harmony", dims = 1:30)

for (res in c(0.4, 0.6, 0.8, 1.0, 1.2)) {
    xenium <- FindClusters(xenium, resolution = res, verbose = FALSE)
}

# -----------------------------------------------------------------------------
# 9. Silhouette evaluation
# -----------------------------------------------------------------------------
# Silhouette score requires a pairwise distance matrix; computing this for all
# 321k cells is infeasible, so we subsample 50,000 cells. Resolution 0.8 was
# selected as it provided the highest mean silhouette while preserving
# biological separations relevant to AIN (T cell subsets, MP polarization
# states, PT injury states).

set.seed(42)
sil_cells <- sample(colnames(xenium), 50000)
embed     <- Embeddings(xenium, "harmony")[sil_cells, 1:30]
d_mat     <- dist(embed)

sil_results <- data.frame(resolution      = numeric(),
                          n_clusters      = integer(),
                          mean_silhouette = numeric())

for (res in c(0.4, 0.6, 0.8, 1.0, 1.2)) {
    col_name <- paste0("SCT_snn_res.", res)
    clusters <- as.integer(xenium@meta.data[sil_cells, col_name])
    sil      <- silhouette(clusters, d_mat)
    sil_results <- rbind(sil_results,
                         data.frame(resolution      = res,
                                    n_clusters      = length(unique(clusters)),
                                    mean_silhouette = mean(sil[, "sil_width"])))
}
print(sil_results)

# Set the chosen resolution as the active identity
Idents(xenium) <- "SCT_snn_res.0.8"

# -----------------------------------------------------------------------------
# 10. Find cluster markers and annotate cell types
# -----------------------------------------------------------------------------
# FindAllMarkers identifies positive markers for each unsupervised cluster.
# These DEGs, combined with spatial localization in Xenium Explorer (see
# script's commented section below), drove the manual annotation. The
# unsupervised clusters were collapsed into 26 core cell types and one mixed
# population, recorded in `final_annotations2`.

all_markers <- FindAllMarkers(xenium,
                              assay = "Xenium",
                              only.pos       = TRUE,
                              min.pct        = 0.25,
                              logfc.threshold = 0.25)
write.csv(all_markers, "xenium_cluster_markers_res0.8.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 11. Cell type annotation
# -----------------------------------------------------------------------------
# Annotation was performed manually based on (i) cluster DEGs from step 10,
# (ii) spatial localization in Xenium Explorer (clusters were exported as
# CSVs and overlaid on the H&E images), and (iii) canonical markers for
# kidney-resident populations:
#
#    PT (healthy)         LRP2, CUBN, SLC34A1, SLC22A6
#    PT (dedifferentiated) LRP2 low, CUBN low, no injury markers
#    PT (injured)         HAVCR1, VCAM1, LCN2
#    TAL                  UMOD, SLC12A1
#    DCT                  SLC12A3
#    Collecting duct      AQP2, SLC4A1, ATP6V0D2
#    Endothelial          PECAM1, PLVAP
#    Stromal/Fibroblast   ACTA2, COL1A1
#    Macrophage           CD163, CD68, C1QA
#    T cells              CD3E, CD4, CD8A
#    B / Plasma           MS4A1, MZB1
#
# The collapsed annotation is stored as `final_annotations2`. See
# Supplementary Table 6 for the full mapping of unsupervised clusters
# to the 26 core cell types and one mixed population.

# Example mapping (illustrative; the exact cluster IDs from a fresh run will
# differ since clustering is stochastic):
# annotation_map <- c(
#     "0" = "TAL", "1" = "PT-H", "2" = "PT-INJ", ... etc
# )
# xenium$final_annotations2 <- annotation_map[as.character(Idents(xenium))]

# -----------------------------------------------------------------------------
# 12. Targeted subclustering of immune lineages
# -----------------------------------------------------------------------------
# T cell, macrophage, and dendritic cell clusters were each separately subset
# from the full object and re-clustered using lineage-specific marker panels
# (see Methods, "Immune cell subclustering"):
#   - T cells:       CD4, CD8A, FOXP3, CCR7, CD69, LAG3, TIGIT, PDCD1,
#                    PRF1, GZMB, STAT1, GBP5
#   - Macrophages:   STAT1, IRF1, CXCL9, CXCL10, GBP5, CD74, LIPA, SPP1,
#                    APOE, TREM2, LYVE1, FOLR2, CD163, MRC1, COL1A1, HIF1A,
#                    LDHA, MKI67
#
# Subcluster annotations are stored as `final_cd4`, `final_cd8`, and
# `macrophage_published_subtype`.

# -----------------------------------------------------------------------------
# 13. Save final clustered + annotated object
# -----------------------------------------------------------------------------
saveRDS(xenium, "xenium_AIN_clustered_annotated.rds")
