# =============================================================================
# 01_xenium_load_qc_merge.R
# -----------------------------------------------------------------------------
# Code for: Baker, Kakade et al., "Spatial analysis reveals the cellular
# microenvironments and mechanisms of inflammation and kidney injury in acute
# interstitial nephritis." Nature Communications.
#
# This script documents the workflow used to go from raw Xenium Analyzer
# output to a single merged Seurat v5 object with patient-level metadata
# attached. It mirrors the pipeline described in the Methods section under
# "Spatial transcriptomic analysis."
#
# Inputs (per slide):
#   - Xenium Analyzer output directory containing experiment.xenium,
#     cell_feature_matrix.h5, cells.csv.gz, transcripts.parquet, etc.
#   - Per-patient cell-ID lists exported from 10x Xenium Explorer using
#     lasso selection (one xlsx per tissue / patient).
#
# Output:
#   - A merged Seurat v5 object (Xenium assay only, raw counts) with
#     base_patient_id, patient_category, and tissue_region annotations.
# =============================================================================

suppressPackageStartupMessages({
    library(Seurat)
    library(readxl)
    library(dplyr)
    library(ggplot2)
})

set.seed(42)

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------
# Each slide's Xenium output is in its own directory. The AIN study used
# multiple Xenium acquisitions across the 20 sample tissues (8 AIN, 7 ATI,
# 5 reference) -- the variable below lists those output directories.
# Per-tissue cell-ID xlsx files (created by lasso selection in Xenium Explorer)
# live alongside each acquisition.
#
# File naming convention for cell-ID xlsx files:
#   {batch}_{patient_category}_{tissue_id}.xlsx
# e.g., "Batch3_AIN_23929.xlsx"  ->  patient_category = AIN,
#                                    tissue_id        = 23929,
#                                    base_patient_id  = AIN_3_23929

xenium_outputs <- list(
    list(slide_id  = "Batch1",
         data_dir  = "/path/to/output-XETG..._Region_1__.../",
         id_dir    = "/path/to/Batch1_patient_ids/"),
    list(slide_id  = "Batch2",
         data_dir  = "/path/to/output-XETG..._Region_2__.../",
         id_dir    = "/path/to/Batch2_patient_ids/")
    # ... continue for all Xenium acquisition runs
)

# -----------------------------------------------------------------------------
# 2. Load each acquisition with LoadXenium and tag the slide
# -----------------------------------------------------------------------------
# LoadXenium() reads the Xenium output directly. We tag each loaded object
# with its batch identifier, then later add per-tissue patient annotations
# from the Xenium Explorer-derived xlsx files.

load_one_acquisition <- function(acq) {
    cat("Loading", acq$slide_id, "from", acq$data_dir, "\n")
    obj <- LoadXenium(data.dir = acq$data_dir, fov = "fov")
    obj$slide <- acq$slide_id
    obj
}

xenium_list <- lapply(xenium_outputs, load_one_acquisition)

# -----------------------------------------------------------------------------
# 3. Assign patient metadata from per-tissue xlsx files
# -----------------------------------------------------------------------------
# For each Xenium acquisition, we lasso-selected each tissue in Xenium Explorer
# and exported the cell IDs to xlsx. This function walks those files, parses
# the filenames to extract category + tissue ID, and attaches the metadata
# to the Seurat object.
#
# The xlsx files contain a `Cell_ID` column whose values match the original
# Xenium cell barcodes (e.g. "aaaacgde-1"). Patient assignment must happen
# BEFORE adding any slide prefix to the barcodes (step 5), since the xlsx
# files contain unprefixed barcodes.

assign_patient_metadata <- function(seurat_obj, id_dir) {

    xlsx_files <- list.files(id_dir, pattern = "\\.xlsx$", full.names = TRUE)
    cat(sprintf("  Found %d patient files\n", length(xlsx_files)))

    seurat_obj$patient_category <- NA_character_
    seurat_obj$base_patient_id  <- NA_character_

    for (f in xlsx_files) {
        # Parse: {batch}_{category}_{tissue_id}.xlsx
        parts <- strsplit(gsub("\\.xlsx$", "", basename(f)), "_")[[1]]
        batch        <- parts[1]
        category     <- parts[2]                                 # AIN / ATI / REF
        tissue_id    <- paste(parts[3:length(parts)], collapse = "_")
        patient_id   <- paste(category, batch, tissue_id, sep = "_")

        cell_ids <- read_excel(f, col_types = "text")[[1]]
        match_idx <- which(colnames(seurat_obj) %in% cell_ids)

        if (length(match_idx) > 0) {
            seurat_obj$patient_category[match_idx] <- category
            seurat_obj$base_patient_id[match_idx]  <- patient_id
        }
    }
    seurat_obj
}

# Apply to each acquisition
for (i in seq_along(xenium_list)) {
    xenium_list[[i]] <- assign_patient_metadata(
        xenium_list[[i]],
        xenium_outputs[[i]]$id_dir
    )
}

# -----------------------------------------------------------------------------
# 4. Drop cells that fell between tissues (no patient assignment)
# -----------------------------------------------------------------------------
# Cells in the inter-tissue space don't belong to any patient. After lasso
# selection these typically represent ~0.5-1.5% of total cells and are removed.
xenium_list <- lapply(xenium_list, function(obj) {
    subset(obj, !is.na(base_patient_id))
})

# -----------------------------------------------------------------------------
# 5. Add slide prefix to cell barcodes so they're unique after merging
# -----------------------------------------------------------------------------
# Xenium uses a per-acquisition barcode namespace -- the same barcode like
# "aaaacgde-1" can exist on different slides as completely different cells.
# Prefixing with the slide ID makes barcodes globally unique.
for (i in seq_along(xenium_list)) {
    xenium_list[[i]] <- RenameCells(
        xenium_list[[i]],
        add.cell.id = xenium_outputs[[i]]$slide_id
    )
}

# -----------------------------------------------------------------------------
# 6. Tag tissue_region (Cortex vs Medulla) per tissue
# -----------------------------------------------------------------------------
# Tissue region was annotated on the matched H&E by a renal pathologist.
# This information is stored in a sample sheet keyed by base_patient_id;
# the mapping is applied to each cell.
sample_sheet <- read.csv("/path/to/sample_sheet.csv")
# Expected columns: base_patient_id, tissue_region (Cortex / Medulla)

xenium_list <- lapply(xenium_list, function(obj) {
    region_lookup <- setNames(sample_sheet$tissue_region,
                              sample_sheet$base_patient_id)
    obj$tissue_region <- unname(region_lookup[obj$base_patient_id])
    obj
})

# -----------------------------------------------------------------------------
# 7. Merge into a single Seurat object
# -----------------------------------------------------------------------------
xenium <- merge(xenium_list[[1]],
                y = xenium_list[-1],
                project = "AIN_Xenium")

# Note: When merging multiple Xenium objects, Seurat resolves FOV key naming
# automatically -- you may see warnings like "Key 'Xenium_' taken, using
# 'fov2_' instead". This is expected.

cat("Merged object cells:",  ncol(xenium), "\n")
cat("Merged object genes:",  nrow(xenium), "\n")
cat("\nCells per category:\n")
print(table(xenium$patient_category))
cat("\nCells per tissue region:\n")
print(table(xenium$tissue_region))

# -----------------------------------------------------------------------------
# 8. Save the raw merged object
# -----------------------------------------------------------------------------
# This is the QC-unfiltered merged object that becomes the input to script 02.
saveRDS(xenium, "xenium_AIN_raw_merged.rds")

# -----------------------------------------------------------------------------
# 9. Spatial verification
# -----------------------------------------------------------------------------
# Generates a per-slide PDF showing every cell colored by patient_category
# with patient-id labels at tissue centroids. Used to visually confirm that
# patient assignment from the xlsx files matched the correct anatomic regions
# of each tissue. (Run interactively after the merge.)

verify_patient_assignment <- function(seurat_obj, slide_name, out_pdf) {
    coords <- GetTissueCoordinates(seurat_obj)
    meta   <- seurat_obj@meta.data
    meta$cell <- colnames(seurat_obj)
    meta <- merge(meta, coords, by = "cell")

    label_pos <- meta %>%
        group_by(base_patient_id, patient_category) %>%
        summarise(x_center = median(x),
                  y_center = median(y),
                  .groups = "drop")

    p <- ggplot(meta, aes(x = x, y = -y)) +
        geom_point(aes(color = patient_category),
                   size = 0.01, alpha = 0.1) +
        geom_label(data = label_pos,
                   aes(x = x_center, y = -y_center,
                       label = base_patient_id),
                   size = 2, fontface = "bold") +
        coord_equal() + theme_minimal() +
        labs(title = paste(slide_name, "patient-ID verification"))

    ggsave(out_pdf, p, width = 16, height = 12)
}
# verify_patient_assignment(xenium_list[[1]], "Batch1", "verification_batch1.pdf")
