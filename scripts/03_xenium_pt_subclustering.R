# =============================================================================
# 03_xenium_pt_subclustering.R
# -----------------------------------------------------------------------------
# Code for: Baker, Kakade et al., Nature Communications.
#
# This script documents the proximal tubule (PT) injury-state analyses,
# including (i) targeted PT subclustering, (ii) HAVCR1/VCAM1 marker-based
# stratification, and (iii) cross-cohort validation against the Kidney
# Precision Medicine Project (KPMP) snRNA-seq atlas (Lake et al. Nature 2023).
#
# These analyses correspond to the Methods sections "Proximal tubule
# subclustering" and "Cross-platform PT validation against KPMP snRNA-seq."
# The figures supported by this code include Fig. 2, Extended Data Fig. 6,
# and Supplementary Fig. 1.
#
# Inputs:
#   - Clustered Seurat object from script 02
#     (xenium_AIN_clustered_annotated.rds)
#   - KPMP h5ad downloaded from CZ CELLxGENE Discover collection
#     bcb61471-2a44-4d00-a0af-ff085512674c
#
# Outputs:
#   - PT subcluster identities merged back to the full object as
#     `PT_matched_orig_cluster` and `PT_matched_orig_annotation`
#   - Stratification analysis tables and validation against KPMP
# =============================================================================

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(reticulate)   # for reading the KPMP .h5ad
    library(anndata)
})

set.seed(42)

# -----------------------------------------------------------------------------
# 1. Load the full annotated Xenium object
# -----------------------------------------------------------------------------
xenium <- readRDS("xenium_AIN_clustered_annotated.rds")

# -----------------------------------------------------------------------------
# 2. Subset all PT cells from the full object
# -----------------------------------------------------------------------------
# All cells annotated as proximal tubule in `final_annotations2` are extracted
# for dedicated subclustering. Per the Methods: "All cells annotated as PT in
# the primary clustering step were extracted, re-normalized using SCTransform,
# batch-corrected using Harmony, and reclustered using the shared
# nearest-neighbor graph approach with Louvain community detection, identical
# to the parameters used for primary clustering."

pt_cells <- subset(xenium,
                   subset = final_annotations2 %in% c("PT-H", "PT-INJ", "PT-dd"))
cat("PT cells extracted:", ncol(pt_cells), "\n")

# -----------------------------------------------------------------------------
# 3. Re-normalize and re-integrate the PT subset
# -----------------------------------------------------------------------------
# Re-running SCTransform on just the PT compartment finds variable features
# that distinguish PT subtypes from each other rather than PT cells from
# other lineages. Harmony is re-applied to correct any residual patient-level
# variation within the PT subset.
pt_cells[["Xenium"]] <- JoinLayers(pt_cells[["Xenium"]])
pt_cells[["Xenium"]] <- split(pt_cells[["Xenium"]], f = pt_cells$base_patient_id)

pt_cells <- SCTransform(pt_cells, assay = "Xenium", verbose = FALSE)
pt_cells <- RunPCA(pt_cells, verbose = FALSE)

pt_cells <- IntegrateLayers(
    object         = pt_cells,
    method         = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction  = "harmony",
    verbose        = FALSE
)
pt_cells[["Xenium"]] <- JoinLayers(pt_cells[["Xenium"]])

pt_cells <- RunUMAP(pt_cells, reduction = "harmony", dims = 1:30)
pt_cells <- FindNeighbors(pt_cells, reduction = "harmony", dims = 1:30)
pt_cells <- FindClusters(pt_cells, resolution = 0.8)

# Result: 21 PT subclusters (matching Supplementary Table 7).

# -----------------------------------------------------------------------------
# 4. Collapse 21 subclusters to 3 PT states
# -----------------------------------------------------------------------------
# The 21 unsupervised PT subclusters were classified into three states based
# on canonical marker expression. Per the Methods:
#
#    PT-H   (n=7 subclusters): LRP2-high, CUBN-high, preserved transport
#                              machinery, minimal injury markers
#    PT-dd  (n=4 subclusters): LRP2-low, CUBN-low, without HAVCR1, VCAM1,
#                              or LCN2
#    PT-INJ (n=10 subclusters): HAVCR1, VCAM1, and/or LCN2 expression
#
# The exact cluster->state mapping is in Supplementary Table 7. Below is
# the structure of that mapping (cluster IDs from a fresh run will differ
# because clustering is stochastic; the mapping is built by inspecting
# DEGs and marker expression cluster by cluster).

# Example structure (illustrative — see Supplementary Table 7 for the
# actual cluster->state assignments):
# pt_state_map <- c(
#     "0"  = "PT-INJ", "1"  = "PT-H",  "2" = "PT-INJ", "3" = "PT-H",
#     "4"  = "PT-H",   "5"  = "PT-dd", ...
# )
# pt_cells$pt_state <- pt_state_map[as.character(Idents(pt_cells))]

# -----------------------------------------------------------------------------
# 5. Map PT subcluster identities back to the full object
# -----------------------------------------------------------------------------
# Both the granular subcluster ID and the collapsed state are mapped back to
# the full xenium object as `PT_matched_orig_cluster` and
# `PT_matched_orig_annotation`. Non-PT cells get NA in these columns.
xenium$PT_matched_orig_cluster    <- NA_integer_
xenium$PT_matched_orig_annotation <- NA_character_

xenium$PT_matched_orig_cluster[colnames(pt_cells)] <-
    as.integer(as.character(Idents(pt_cells)))
# xenium$PT_matched_orig_annotation[colnames(pt_cells)] <- pt_cells$pt_state

# -----------------------------------------------------------------------------
# 6. HAVCR1/VCAM1 stratification within AIN cortex PT cells
# -----------------------------------------------------------------------------
# Per the Methods: "To parallel the IMC-based PT state definitions (which, in
# part, used KIM-1 and VCAM1 protein expression to define PT-INJ1 and
# PT-INJ2), cortical PT cells from AIN biopsies were further stratified by
# HAVCR1 and VCAM1 mRNA expression status (expression > 0 defining
# positivity), yielding four cross-classified groups (PTneg, HAVCR1+VCAM1-,
# HAVCR1-VCAM1+, and HAVCR1+VCAM1+)."

ain_cortex_pt <- subset(xenium,
                        subset = patient_category   == "AIN" &
                                 tissue_region      == "Cortex" &
                                 final_annotations2 %in% c("PT-H","PT-INJ","PT-dd"))

havcr1_pos <- FetchData(ain_cortex_pt, vars = "HAVCR1", layer = "data")[,1] > 0
vcam1_pos  <- FetchData(ain_cortex_pt, vars = "VCAM1",  layer = "data")[,1] > 0

ain_cortex_pt$marker_group <- paste0(
    ifelse(havcr1_pos, "HAVCR1pos", "HAVCR1neg"), "_",
    ifelse(vcam1_pos,  "VCAM1pos",  "VCAM1neg")
)

cat("AIN cortex PT marker groups:\n")
print(table(ain_cortex_pt$marker_group))

# -----------------------------------------------------------------------------
# 7. Differential expression analyses between marker groups
# -----------------------------------------------------------------------------
# Per the Methods: "DEGs identified between HAVCR1+VCAM1- and HAVCR1-VCAM1-
# (PTneg), and between VCAM1+ and PTneg, using Seurat's FindMarkers with
# Wilcoxon rank-sum test. Parameters: min.pct = 0.10, logfc.threshold = 0.25,
# Bonferroni correction."

Idents(ain_cortex_pt) <- "marker_group"

deg_inj1 <- FindMarkers(ain_cortex_pt,
                        ident.1 = "HAVCR1pos_VCAM1neg",
                        ident.2 = "HAVCR1neg_VCAM1neg",
                        assay = "SCT",
                        min.pct = 0.10,
                        logfc.threshold = 0.25)

deg_inj2 <- FindMarkers(ain_cortex_pt,
                        ident.1 = "VCAM1pos",                # any VCAM1+
                        ident.2 = "HAVCR1neg_VCAM1neg",
                        assay = "SCT",
                        min.pct = 0.10,
                        logfc.threshold = 0.25)

write.csv(deg_inj1, "PT_INJ1_vs_PTneg_AIN_cortex.csv")
write.csv(deg_inj2, "PT_INJ2_vs_PTneg_AIN_cortex.csv")

# -----------------------------------------------------------------------------
# 8. Cross-platform validation against KPMP (Lake et al. Nature 2023)
# -----------------------------------------------------------------------------
# The KPMP atlas was used as an independent validation cohort. Per the
# Methods: "The full annotated dataset was obtained from the CZ CELLxGENE
# Discover collection (collection ID bcb61471-2a44-4d00-a0af-ff085512674c)
# as an AnnData (.h5ad) object and processed in R using the anndata and
# reticulate packages."

kpmp <- read_h5ad("/path/to/kpmp_lake2023_full.h5ad")

# Restrict to PT cells (Lake annotations: PT-S1/2, PT-S3, aPT, dPT)
pt_idx     <- kpmp$obs[["subclass.l2"]] %in% c("PT-S1/2", "PT-S3", "aPT", "dPT")
kpmp_pt    <- kpmp[pt_idx, ]
cat("KPMP PT cells:", nrow(kpmp_pt$obs), "\n")

# Map gene symbols (KPMP var_names are Ensembl IDs; feature_name has symbols)
sym_to_idx <- function(sym) which(kpmp_pt$var$feature_name == sym)[1]

havcr1_kpmp <- as.numeric(kpmp_pt$X[, sym_to_idx("HAVCR1")]) > 0
vcam1_kpmp  <- as.numeric(kpmp_pt$X[, sym_to_idx("VCAM1")])  > 0

# Cross-classify: Lake annotation collapsed to healthy vs injured, plus
# HAVCR1/VCAM1 expression status (parallels the Xenium stratification)
lake_state <- ifelse(kpmp_pt$obs[["subclass.l2"]] %in% c("aPT","dPT"),
                     "injured", "healthy")
kpmp_groups <- paste0(lake_state, "_",
                      ifelse(havcr1_kpmp, "HAVCR1pos", "HAVCR1neg"), "_",
                      ifelse(vcam1_kpmp,  "VCAM1pos",  "VCAM1neg"))

cat("KPMP cross-classification:\n")
print(table(kpmp_groups))

# -----------------------------------------------------------------------------
# 9. Stepwise injury-marker analysis ("staircase" plots)
# -----------------------------------------------------------------------------
# Per Methods: for both Xenium and KPMP, within each of the four cross-
# classified groups, the percentage of cells expressing each injury-
# associated transcript was computed. A cell was counted positive for a
# given gene if normalized expression > 0.

target_genes <- c("C3", "VIM", "STAT1", "NAMPT", "THY1", "CD74",
                  "CUBN", "LRP2")

# Xenium staircase
ain_cortex_pt$pt_state_for_staircase <- ifelse(
    ain_cortex_pt$final_annotations2 == "PT-INJ",
    "PT_INJ",
    "PT_H"   # PT-H and PT-dd grouped here for the four-group analysis;
             # see Methods for the exact n's reported in Supplementary Fig. 1
)
ain_cortex_pt$staircase_group <- paste0(
    ain_cortex_pt$pt_state_for_staircase, "_",
    ain_cortex_pt$marker_group
)

xenium_staircase <- lapply(target_genes, function(g) {
    expr <- FetchData(ain_cortex_pt, vars = g, layer = "data")[,1]
    pos  <- expr > 0
    data.frame(
        gene  = g,
        group = names(table(ain_cortex_pt$staircase_group)),
        pct_positive = sapply(names(table(ain_cortex_pt$staircase_group)),
                              function(grp) {
                                  100 * mean(pos[ain_cortex_pt$staircase_group == grp])
                              })
    )
}) |> bind_rows()

write.csv(xenium_staircase, "staircase_xenium.csv", row.names = FALSE)

# KPMP staircase (same logic, using kpmp_pt expression directly)
kpmp_staircase <- lapply(target_genes, function(g) {
    idx <- sym_to_idx(g)
    if (is.na(idx)) {
        return(NULL)   # SOD2 is absent from the KPMP gene panel
    }
    pos <- as.numeric(kpmp_pt$X[, idx]) > 0
    data.frame(
        gene  = g,
        group = names(table(kpmp_groups)),
        pct_positive = sapply(names(table(kpmp_groups)), function(grp) {
            100 * mean(pos[kpmp_groups == grp])
        })
    )
}) |> bind_rows()

write.csv(kpmp_staircase, "staircase_kpmp.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 10. Independent-axes ("dumbbell") analysis
# -----------------------------------------------------------------------------
# Per the Methods: "To test whether HAVCR1 and VCAM1 define independent
# biological axes within the injured PT compartment, we computed Spearman
# rank correlations between the percent-positivity of each marker gene and
# the percent-positivity of HAVCR1 or VCAM1 across the 10 unsupervised
# PT-INJ-dominant Xenium subclusters."
#
# This uses the granular `PT_matched_orig_cluster` (21 levels), filtered to
# the 10 PT-INJ-dominant subclusters identified in step 4. Spearman
# correlations of marker positivity vs HAVCR1+ and VCAM1+ positivity across
# subclusters reveal which transcripts track preferentially with one marker
# vs the other.
#
# (Implementation omitted here -- the structure is straightforward: compute
# per-subcluster marker positivity for each gene, then cor.test() against
# HAVCR1/VCAM1 positivity.)

# -----------------------------------------------------------------------------
# 11. Save updated full object
# -----------------------------------------------------------------------------
saveRDS(xenium, "xenium_AIN_with_pt_subclusters.rds")
