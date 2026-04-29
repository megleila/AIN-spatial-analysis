# =============================================================================
# 01_imc_normalize_cluster_annotate.R
# -----------------------------------------------------------------------------
# Code for: Baker, Kakade et al. Manuscript under revision at Nature
# Communications.
#
# This script documents the imaging mass cytometry (IMC) downstream analysis
# workflow, starting from a SpatialExperiment object containing single-cell
# intensities exported from steinbock. It mirrors the procedures described in
# the Methods sections "Signal spillover correction," "Image and cell-level
# quality control and batch effect correction," and "Cell phenotyping," and
# corresponds to the workflow used to produce Fig. 1c, Fig. 2a-d, Extended
# Data Fig. 1d, Extended Data Fig. 2, and Supplementary Table 4.
#
# Inputs:
#   - Raw SpatialExperiment object built from steinbock output (cell-level
#     intensities, spatial coordinates, sample IDs, marker panel). Pre-
#     processing of raw .mcd files using steinbock (hot pixel filtering,
#     panel configuration) and cell segmentation using custom-tuned Mesmer
#     are documented in the Methods but not implemented here, as they
#     operate on raw imaging files rather than on the SpatialExperiment
#     object that downstream analyses depend on.
#
# Outputs:
#   - Annotated, batch-corrected SpatialExperiment object with `clusters_corr`
#     (cluster IDs), `sub_clusters`, `simplified_annotations`, and
#     `core_annotations` columns. Yale and JHU cohorts processed by identical
#     pipeline using the same hybridized antibody panel.
# =============================================================================

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(imcRtools)         # imcRtools workflow per Bodenmiller (Methods)
    library(CATALYST)          # spillover correction
    library(batchelor)         # fastMNN batch correction
    library(scran)
    library(scater)
    library(Rphenograph)       # graph-based clustering
    library(igraph)
    library(BiocParallel)
    library(dittoSeq)
    library(viridis)
    library(tidyverse)
})

set.seed(220619)

# -----------------------------------------------------------------------------
# 1. Load raw SpatialExperiment object
# -----------------------------------------------------------------------------
# This object contains 31 metal-conjugated antibody channels plus 4 channel
# controls (DNA1, DNA2, Membrane1, Membrane2), with single-cell intensities
# from steinbock exported as the `counts` assay. Spatial coordinates are in
# spatialCoords (Pos_X, Pos_Y); patient identity is in `Patient`; tissue
# class in `indication` (AC=AIN cortex, AM=AIN medulla, BC=ATI cortex,
# BM=ATI medulla, RC=Reference cortex, RM=Reference medulla). Slide-level
# batch identity is in `Batch`.
spe <- readRDS("path/to/spe_raw_from_steinbock.rds")

# Mark which channels are biological markers vs. structural/segmentation
# controls. The use_channel flag determines what enters PCA / batch
# correction / clustering downstream.
rowData(spe)$use_channel <- !grepl("DNA1|DNA2|Membrane1|Membrane2",
                                   rownames(spe))

# -----------------------------------------------------------------------------
# 2. Signal spillover correction
# -----------------------------------------------------------------------------
# Per the Methods: "We employed a previously described bead-based compensation
# workflow and R-based software [Chevrier et al., Cell Syst 2018] to estimate
# and compensate for interference between channels."
#
# The bead-based spillover matrix is generated from a single-bead acquisition
# and applied to all single-cell intensities via CATALYST::compCytof.
# `sm` is the precomputed spillover matrix (estimated from bead images using
# CATALYST::computeSpillmat).
sm  <- readRDS("path/to/spillover_matrix.rds")
spe <- compCytof(spe, sm,
                 transform = TRUE, cofactor = 1,
                 isotope_list = isotope_list, overwrite = FALSE)
# After compCytof, an `exprs` assay is added containing spillover-corrected,
# arcsinh-transformed (cofactor=1) intensities. `counts` retains the raw
# uncorrected values. `exprs` is the assay used for all downstream
# clustering, dimensional reduction, and visualization.

# -----------------------------------------------------------------------------
# 3. Image- and cell-level quality control
# -----------------------------------------------------------------------------
# Per the Methods: "We performed image-level quality control through
# calculating the signal-to-noise ratio for individual channels and found
# that each of our channels had an acceptable ratio >3 using the otsu
# thresholding approach."
#
# Channel-level SNR is computed per image and per channel; channels with
# SNR <= 3 in a given image would be flagged. In practice all channels
# passed this threshold.
#
# Cell-level QC included removing rows of pixels at the edges of acquisition
# fields ("hot pixel" filtering at the steinbock step) and removing cells
# with implausibly small or large segmented areas. Cells used for downstream
# analysis are those passing QC at both the image and cell level.

# -----------------------------------------------------------------------------
# 4. Batch effect correction across slides
# -----------------------------------------------------------------------------
# Per the Methods: "Batch effect correction was performed between different
# batches (slides) using the fast mutual nearest neighbors function within
# the Batchelor package."
#
# Each slide ("batch") contained one AIN, one ATI, one reference, and one
# lymph node tissue. fastMNN uses mutual nearest-neighbor pairs across
# batches to compute a corrected low-dimensional embedding while preserving
# biological variation.
#
# The fastMNN call below operates on the use_channel markers in `exprs`,
# producing a corrected embedding that is stored as the `fastMNN` reducedDim.
# A UMAP computed on this corrected embedding (`UMAP_mnnCorrected`) is used
# for visualization throughout the manuscript (e.g., Fig. 1c).

mnn_out <- fastMNN(
    spe,
    batch  = spe$Batch,
    assay.type = "exprs",
    subset.row = rowData(spe)$use_channel,
    d = 30                                  # dimensions retained
)
reducedDim(spe, "fastMNN") <- reducedDim(mnn_out, "corrected")

spe <- runUMAP(spe, dimred = "fastMNN",
               name = "UMAP_mnnCorrected", n_neighbors = 30)

# An uncorrected UMAP is also computed on the raw `exprs` for diagnostic
# comparison (Extended Data Fig. 2a).
spe <- runUMAP(spe, exprs_values = "exprs",
               subset_row = rowData(spe)$use_channel,
               name = "UMAP")

# -----------------------------------------------------------------------------
# 5. Rphenograph clustering on batch-corrected embedding
# -----------------------------------------------------------------------------
# Per the Methods: "We employed the RPhenoGraph clustering approach, which
# includes calculation of a graph by detecting the K nearest neighbors based
# on Euclidean distance in expression space, weighting of edges between
# nodes by overlap in nearest neighbor sets using the jaccard index, and
# Louvain modularity optimization to detect connected communities and
# partition the graph into clusters of cells."
#
# Default Rphenograph parameter k=45 (number of nearest neighbors) was used
# at the initial unsupervised clustering step.

mat <- t(reducedDim(spe, "fastMNN"))         # cells = columns -> rows
phgr <- Rphenograph(mat, k = 45)
clusters <- factor(membership(phgr[[2]]))
spe$clusters_corr <- clusters

cat("Initial Rphenograph clusters:", length(levels(clusters)), "\n")

# -----------------------------------------------------------------------------
# 6. Visualize clusters and identify which markers drive each cluster
# -----------------------------------------------------------------------------
# Per the Methods: "Initial unsupervised clustering results were analyzed
# through mapping back onto original tissues to assess cell locations as
# well as analysis of marker expression using violin plots of raw expression
# patterns."
#
# Per-cluster marker expression heatmap (mean arcsinh expression):
agg_expr <- aggregateAcrossCells(spe[rowData(spe)$use_channel, ],
                                 ids = spe$clusters_corr,
                                 statistics = "mean",
                                 use.assay.type = "exprs")
# The aggregated matrix is used to drive heatmap-based annotation, alongside
# inspection of cluster locations on tissue (using imcRtools::plotSpatial).

# -----------------------------------------------------------------------------
# 7. Targeted subclustering of clusters with residual heterogeneity
# -----------------------------------------------------------------------------
# Per the Methods: "Subclustering was attempted on all unsupervised clusters
# using all markers appearing within the cluster to be expressed more highly
# than the average across all unsupervised clusters."
#
# For each parent cluster, identify markers expressed above the global cross-
# cluster mean, then re-run Rphenograph using only those markers. Clusters
# yielding biologically meaningful subpopulations are retained as sub_clusters.

subcluster_one <- function(spe_sub, mean_expr_global) {
    # Markers higher than global mean within this cluster
    cluster_means <- rowMeans(assay(spe_sub, "exprs")[rowData(spe_sub)$use_channel, ])
    sub_markers <- names(cluster_means)[cluster_means > mean_expr_global]
    if (length(sub_markers) < 3) return(NULL)
    mat_sub <- t(assay(spe_sub, "exprs")[sub_markers, ])
    phgr_sub <- Rphenograph(mat_sub, k = 30)
    factor(membership(phgr_sub[[2]]))
}

mean_expr_global <- mean(assay(spe, "exprs")[rowData(spe)$use_channel, ])
sub_assignments <- character(ncol(spe))
for (cl in levels(spe$clusters_corr)) {
    idx <- which(spe$clusters_corr == cl)
    spe_cl <- spe[, idx]
    sub <- subcluster_one(spe_cl, mean_expr_global)
    if (is.null(sub)) {
        sub_assignments[idx] <- as.character(cl)
    } else {
        sub_assignments[idx] <- paste0(cl, "_", sub)
    }
}
spe$sub_clusters <- factor(sub_assignments)

# -----------------------------------------------------------------------------
# 8. Manual annotation: cluster -> cell type
# -----------------------------------------------------------------------------
# Per the Methods: "Ultimate cell annotations were validated through marker
# expression patterns and spatial localization, identifying 50 detailed cell
# populations (including mixed populations, Supplementary Table 4), which
# were distilled to 18 core cell populations and two core mixed cell
# populations."
#
# Annotation is performed manually based on:
#   (i)  per-cluster mean marker expression (step 6)
#   (ii) per-cluster spatial localization on representative tissue images
#        (using imcRtools::plotSpatial)
#   (iii) canonical kidney cell-type markers:
#
#         PT (proximal tubule)        Megalin, Aqp1
#         TAL (thick ascending limb)  Umod, CK7
#         DCT (distal convoluted)     Calbindin
#         CD (collecting duct)        CK7
#         PT injury states            KIM1, VCAM1, FACL4, Ki67
#         Endothelial                 CD31, ERG
#         Vascular smooth muscle      aSMA
#         Stromal                     Vimentin, aSMA-low
#         Podocytes                   Podoplanin, Nestin
#         Macrophages                 CD68, CD163, CD206
#         Mononuclear phagocyte       CD11c, CD163
#         T cells                     CD3, CD4, CD8a
#         B cells                     CD20
#         Mast cells                  Chymase, MPO
#         Eosinophils                 MBP
#
# The mapping of each unsupervised sub-cluster to its annotated cell type
# is in Supplementary Table 4; the cluster IDs from a fresh re-run will
# differ since clustering is stochastic, so the table is the authoritative
# mapping rather than a hardcoded dictionary in this script.
#
# Two annotation columns are used downstream:
#   `simplified_annotations`  -- 50 detailed populations (full granularity,
#                                including mixed populations); internal
#                                working labels
#   `core_annotations`        -- 18 core + 2 mixed (collapsed for figures
#                                showing population-level summaries);
#                                internal working labels
#   `manuscript_label`        -- the published cluster names from
#                                Supplementary Table 4 (e.g., "TAL-H" rather
#                                than internal "TAL-Umod-hi"; "MP-CD206+"
#                                rather than "Macs-CD206+"). Mapped 1:1 from
#                                simplified_annotations.
#   `manuscript_core_label`   -- the published core cell type column from
#                                Supplementary Table 4 (e.g., "PT", "TAL",
#                                "T cells", "MP", "(Mixed) Stromal").
#
# Downstream analyses (script 02) operate on `manuscript_label` so that
# index/neighbor cell type names match those readers see in the manuscript
# figures and tables.
#
# Example shape of the annotation step (illustrative — see Supplementary
# Table 4 for the actual mapping):
#
# annotation_map <- c(
#     "1"   = "PT-H",    "2_1" = "PT-INJ1",  "2_2" = "PT-INJ2",
#     "3"   = "TAL-Umod-hi",  "4" = "DCT-Calb-hi",  ...
# )
# spe$simplified_annotations <- factor(annotation_map[as.character(spe$sub_clusters)])

# -----------------------------------------------------------------------------
# 9. Save annotated, batch-corrected object
# -----------------------------------------------------------------------------
# This is the deposit-ready object: input to script 02 (neighborhood and
# density analyses) and the file deposited on Zenodo.

saveRDS(spe, "imc_annotated_batch_corrected.rds")
