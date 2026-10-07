###############################################################################
## Fig 4C-E and Fig S5D-F: MHC II+ / MHC II- megakaryocyte bulk RNA-seq
## integrated with the Sun et al. (2021) single-cell MK atlas
##
##   Fig 4C   unified panel across the reference subset pseudobulks
##   Fig 4D   unified panel across the bulk MHC II+ and MHC II- samples
##   Fig 4E   Spearman correlation of each bulk replicate to each subset
##   Fig S5D  validation control (subset pseudo-replicates)
##   Fig S5E  MHC II+ / MHC II- signature scores by subset
##   Fig S5F  the same scores on the reference UMAP
##
## Steps
##   1. Annotate the four reference MK subsets (fixed cluster map + marker QC)
##   2. Load bulk counts and the MHC II+ / MHC II- DE signatures
##   3. Pseudobulk each reference subset (3 pseudo-replicates per subset)
##   4. Build a unified gene panel (subset markers + MHC II machinery)
##   5. Panel heatmaps for the reference pseudobulks and the bulk samples
##   6. Spearman mapping of each bulk replicate to the subsets
##   7. MHC II+ / MHC II- signature scores across reference cells
##   8. Validation control: subset pseudo-replicates through step 6
##
## Run from the repository root, after the Fig 4A/S5A script (it writes the DE table):
##   Rscript figures/Fig4C-E_S5D-F_sun_atlas_integration/Fig4C-E_S5D-F_sun_atlas_integration.R
## Inputs are described in data/README.md. Outputs are written to
## results/Fig4C-E_S5D-F_sun_atlas_integration/.
###############################################################################

## ---------------------------------------------------------------------------
## 0. CONFIG
## Paths are relative to the repository root. Override any of them with the
## environment variables shown, e.g. MK_SC_RDS=/path/to/mouse_mk.rds
## ---------------------------------------------------------------------------
SC_RDS      <- Sys.getenv("MK_SC_RDS",      "data/mouse_mk.rds")
BULK_COUNTS <- Sys.getenv("MK_BULK_COUNTS", "data/counts_matrix.tsv")
DE_TABLE    <- Sys.getenv("MK_DE_TABLE",    "results/Fig4A_S5A_differential_expression/DE_genes_sig.tsv")
OUTDIR      <- Sys.getenv("MK_OUTDIR",      "results/Fig4C-E_S5D-F_sun_atlas_integration")

CLUSTER_COL   <- "RNA_snn_res.0.1"   # clustering that defines the 4 subsets
N_PSEUDOREP   <- 3                   # pseudo-replicates per subset
N_TOP_MARKERS <- 12                  # subset-defining markers per subset
N_HVG         <- 2000                # shared HVGs used for correlation
SEED          <- 13579               # seed used for the reported results
EXPORT_TIFF   <- FALSE               # TRUE also writes 300 dpi TIFFs

## Cluster-to-subset map, fixed from per-cluster marker expression
## (written to qc_cluster_marker_means.csv in the output folder on every run):
##   0  Top2a / Birc5 high, Cxcl12 absent     -> Active-Cycling
##   1  Cxcl12 / Mpl high                     -> HSC-Niche
##   2  H2-Ab1 / Cd74 / Irf7 high             -> Inflammatory-Immune
##   3  Tubb1 / Gp1ba / Ppbp high             -> Platelet-producing
## Cluster IDs can change if the reference object is rebuilt. Re-check the
## QC table and update this map before trusting any downstream output.
SUBSET_MAP <- c("0" = "Active-Cycling",
                "1" = "HSC-Niche",
                "2" = "Inflammatory-Immune",
                "3" = "Platelet-producing")
SUBSET_LEVELS <- c("Platelet-producing", "HSC-Niche",
                   "Inflammatory-Immune", "Active-Cycling")

## ---------------------------------------------------------------------------
## 0b. LIBRARIES AND HELPERS
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(edgeR)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

if (!file.exists(DE_TABLE)) {
  stop("DE table not found: ", DE_TABLE,
       "\nRun figures/Fig4A_S5A_differential_expression/Fig4A_S5A_differential_expression.R first,",
       "\nor set MK_DE_TABLE to an existing DE_genes_sig.tsv.")
}
for (f in c(SC_RDS, BULK_COUNTS)) {
  if (!file.exists(f)) stop("Input not found: ", f, "\nSee data/README.md.")
}
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

## Raw counts matrix, Seurat v4 or v5
get_counts <- function(obj, assay = "RNA") {
  m <- tryCatch(LayerData(obj, assay = assay, layer = "counts"),
                error = function(e) NULL)
  if (is.null(m)) {
    obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
    m <- tryCatch(LayerData(obj, assay = assay, layer = "counts"),
                  error = function(e) GetAssayData(obj, assay = assay, slot = "counts"))
  }
  as(m, "CsparseMatrix")
}

## "ENSMUSG..._Symbol" -> "Symbol"
sym_from_rowname <- function(x) sub(".*_(.*)$", "\\1", x)

## log2 counts-per-million (edgeR)
logcpm <- function(m) cpm(DGEList(counts = m), log = TRUE, prior.count = 1)

## Row z-score with zero-variance rows set to 0
zrow <- function(m) { z <- t(scale(t(m))); z[is.na(z)] <- 0; z }

## Output writers (PDF always, TIFF optional)
gg_out <- function(p, stem, width, height) {
  ggsave(file.path(OUTDIR, paste0(stem, ".pdf")), p, width = width, height = height)
  if (EXPORT_TIFF)
    ggsave(file.path(OUTDIR, paste0(stem, ".tiff")), p, width = width,
           height = height, dpi = 300, compression = "lzw")
}
heatmap_out <- function(mat, stem, width, height, ...) {
  pheatmap(mat, filename = file.path(OUTDIR, paste0(stem, ".pdf")),
           width = width, height = height, ...)
  if (EXPORT_TIFF)
    pheatmap(mat, filename = file.path(OUTDIR, paste0(stem, ".tiff")),
             width = width, height = height, ...)
}

ann_colors <- list(
  Subset   = c(`Platelet-producing` = "#FF7F00", `HSC-Niche` = "#FF00FF",
               `Inflammatory-Immune` = "#4DAF4A", `Active-Cycling` = "#E41A1C"),
  Group    = c(MHCpos = "#3182bd", MHCneg = "#f0a30a"),
  Category = c(`Subset signature` = "grey70", `MHC II / co-stim` = "black")
)

###############################################################################
## 1. REFERENCE: ANNOTATE THE FOUR MK SUBSETS
###############################################################################
sc <- readRDS(SC_RDS)
DefaultAssay(sc) <- "RNA"
if (length(VariableFeatures(sc)) == 0) {
  sc <- NormalizeData(sc, verbose = FALSE)
  sc <- FindVariableFeatures(sc, nfeatures = N_HVG, verbose = FALSE)
}

if (!CLUSTER_COL %in% colnames(sc@meta.data)) {
  stop("Cluster column '", CLUSTER_COL, "' not found. Available: ",
       paste(grep("snn_res|clusters", colnames(sc@meta.data), value = TRUE),
             collapse = ", "))
}
clusters <- as.character(sc@meta.data[[CLUSTER_COL]])
if (!setequal(unique(clusters), names(SUBSET_MAP))) {
  stop("Clusters in the object (", paste(sort(unique(clusters)), collapse = ", "),
       ") do not match SUBSET_MAP (", paste(names(SUBSET_MAP), collapse = ", "), ").")
}

mk_subset <- factor(unname(SUBSET_MAP[clusters]), levels = SUBSET_LEVELS)
names(mk_subset) <- colnames(sc)
sc <- AddMetaData(sc, metadata = mk_subset, col.name = "MK_subset")
Idents(sc) <- "MK_subset"
print(table(sc$MK_subset))

## Marker evidence for the map: mean expression of canonical markers per cluster
subset_markers <- list(
  `Platelet-producing`  = c("Pf4", "Ppbp", "Itga2b", "Gp1ba", "Gp9", "Tubb1", "Clec1b", "Vwf"),
  `HSC-Niche`           = c("Cxcl12", "Mpl", "Thbs1", "Igf1", "Sdc1", "Gata2", "Meg3", "Pdzk1ip1"),
  `Inflammatory-Immune` = c("H2-Aa", "H2-Ab1", "Cd74", "B2m", "Tap1", "Irf7", "Lgals3", "Cd53", "Lsp1"),
  `Active-Cycling`      = c("Mki67", "Top2a", "Cdk1", "Ccnb1", "Ccnb2", "Birc5", "Stmn1", "Pclaf"))
subset_markers <- lapply(subset_markers, function(g) intersect(g, rownames(sc)))

qc_means <- AverageExpression(sc, features = unique(unlist(subset_markers)),
                              group.by = CLUSTER_COL, assays = "RNA")$RNA
write.csv(round(as.matrix(qc_means), 3),
          file.path(OUTDIR, "qc_cluster_marker_means.csv"))

dp <- DotPlot(sc, features = unique(unlist(subset_markers)), group.by = "MK_subset") +
  RotatedAxis() + ggtitle("Subset marker check")
gg_out(dp, "qc_subset_marker_dotplot", width = 11, height = 4)

###############################################################################
## 2. BULK COUNTS AND DE SIGNATURES
###############################################################################
bulk_raw <- read.table(BULK_COUNTS, header = TRUE, sep = "\t", row.names = 1,
                       check.names = FALSE)
stopifnot(ncol(bulk_raw) == 8)   # columns 1-4 MHC II+, 5-8 MHC II-
colnames(bulk_raw) <- c(paste0("MHCpos_", 1:4), paste0("MHCneg_", 1:4))
bulk_group <- factor(c(rep("MHCpos", 4), rep("MHCneg", 4)),
                     levels = c("MHCpos", "MHCneg"))

## Gene symbols, immunoglobulin transcripts removed, duplicate symbols summed
bulk_sym <- sym_from_rowname(rownames(bulk_raw))
keep     <- !grepl("^Ighv|^Igha|^Igkv|^Iglv|^Igkc|^Iglc|^Jchain", bulk_sym)
bulk_mat <- rowsum(as.matrix(bulk_raw[keep, ]), group = bulk_sym[keep])
mode(bulk_mat) <- "integer"

de      <- read.delim(DE_TABLE, stringsAsFactors = FALSE)
sig_pos <- de$gene[de$sig == "Up_MHCpos"]
sig_neg <- de$gene[de$sig == "Up_MHCneg"]
message(sprintf("Signatures: MHC II+ %d genes, MHC II- %d genes",
                length(sig_pos), length(sig_neg)))

###############################################################################
## 3. PSEUDOBULK EACH REFERENCE SUBSET
###############################################################################
make_pseudobulk <- function(counts, cell_subset, n_rep, seed) {
  set.seed(seed)
  cells <- data.frame(cell = colnames(counts), subset = as.character(cell_subset),
                      stringsAsFactors = FALSE)
  cells <- cells[!is.na(cells$subset), ]
  cells$rep <- ave(seq_len(nrow(cells)), cells$subset,
                   FUN = function(i) sample(rep_len(seq_len(n_rep), length(i))))
  grp <- factor(paste(cells$subset, cells$rep, sep = "_rep"))
  pb  <- t(rowsum(t(as.matrix(counts[, cells$cell])), group = grp))
  attr(pb, "meta") <- data.frame(group  = colnames(pb),
                                 subset = sub("_rep[0-9]+$", "", colnames(pb)),
                                 stringsAsFactors = FALSE)
  pb
}
pb      <- make_pseudobulk(get_counts(sc), sc$MK_subset, N_PSEUDOREP, SEED)
pb_meta <- attr(pb, "meta")
write.csv(as.data.frame(pb), file.path(OUTDIR, "sun_pseudobulk_counts.csv"))

###############################################################################
## 4. UNIFIED GENE PANEL
###############################################################################
mhc_panel <- intersect(c("H2-Ab1", "H2-Aa", "H2-Eb1", "H2-DMa", "H2-DMb1", "Cd74",
                         "Ciita", "Cd80", "Cd86", "Cd40", "B2m", "Tap1", "Tap2",
                         "Psmb8", "Psmb9"), rownames(sc))

markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
write.csv(markers, file.path(OUTDIR, "sun_subset_markers.csv"), row.names = FALSE)
sig_genes <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = N_TOP_MARKERS) %>%
  pull(gene) %>%
  unique()

panel_genes <- intersect(unique(c(sig_genes, mhc_panel)),
                         intersect(rownames(pb), rownames(bulk_mat)))
gene_cat <- setNames(ifelse(panel_genes %in% mhc_panel, "MHC II / co-stim",
                            "Subset signature"), panel_genes)
write.csv(data.frame(gene = panel_genes, category = gene_cat),
          file.path(OUTDIR, "unified_panel_genes.csv"), row.names = FALSE)
message(sprintf("Unified panel: %d genes (%d MHC II / co-stim)",
                length(panel_genes), sum(gene_cat == "MHC II / co-stim")))

###############################################################################
## 5. PANEL HEATMAPS: REFERENCE PSEUDOBULKS AND BULK SAMPLES
###############################################################################
pb_lcpm   <- logcpm(as.matrix(pb))
bulk_lcpm <- logcpm(bulk_mat)

ann_pb   <- data.frame(Subset = pb_meta$subset, row.names = pb_meta$group)
ann_bulk <- data.frame(Group = bulk_group, row.names = colnames(bulk_mat))
ann_row  <- data.frame(Category = gene_cat, row.names = panel_genes)
div_pal  <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

heatmap_out(zrow(pb_lcpm[panel_genes, ])[, order(pb_meta$subset)],
            "Fig4C_pseudobulk_panel_heatmap", width = 7, height = 9,
            cluster_cols = FALSE, cluster_rows = TRUE,
            annotation_col = ann_pb, annotation_row = ann_row,
            annotation_colors = ann_colors, color = div_pal,
            fontsize_row = 7, border_color = NA,
            main = "Sun et al. pseudobulk, unified panel (row z-score)")

heatmap_out(zrow(bulk_lcpm[panel_genes, ])[, order(bulk_group)],
            "Fig4D_bulk_panel_heatmap", width = 5.5, height = 9,
            cluster_cols = FALSE, cluster_rows = TRUE,
            annotation_col = ann_bulk, annotation_row = ann_row,
            annotation_colors = ann_colors, color = div_pal,
            fontsize_row = 7, border_color = NA,
            main = "Bulk MHC II+/-, unified panel (row z-score)")

###############################################################################
## 6. SPEARMAN MAPPING OF BULK REPLICATES TO SUBSETS
###############################################################################
pb_subset      <- t(rowsum(t(as.matrix(pb)), group = pb_meta$subset))
pb_subset_lcpm <- logcpm(pb_subset)

hvg <- head(intersect(VariableFeatures(sc),
                      intersect(rownames(pb_subset_lcpm), rownames(bulk_lcpm))), N_HVG)
message(sprintf("Correlation over %d shared HVGs", length(hvg)))

cor_mat <- cor(bulk_lcpm[hvg, ], pb_subset_lcpm[hvg, ], method = "spearman")
write.csv(cor_mat, file.path(OUTDIR, "bulk_vs_subset_spearman.csv"))

heatmap_out(cor_mat[order(bulk_group), ], "Fig4E_bulk_vs_subset_correlation",
            width = 7.5, height = 5,
            cluster_rows = FALSE, cluster_cols = TRUE,
            annotation_row = ann_bulk, annotation_colors = ann_colors,
            display_numbers = TRUE, number_format = "%.2f",
            color = colorRampPalette(brewer.pal(9, "Blues"))(100),
            main = "Bulk vs Sun subsets (Spearman rho)")

nearest <- data.frame(sample      = rownames(cor_mat),
                      group       = bulk_group,
                      best_subset = colnames(cor_mat)[max.col(cor_mat)],
                      best_rho    = round(apply(cor_mat, 1, max), 3))
write.csv(nearest, file.path(OUTDIR, "bulk_nearest_subset.csv"), row.names = FALSE)
print(nearest)

###############################################################################
## 7. MHC II+ / MHC II- SIGNATURE SCORES ACROSS REFERENCE CELLS
###############################################################################
sc <- AddModuleScore(sc, features = list(intersect(sig_pos, rownames(sc))),
                     name = "MHCpos_sig", seed = SEED)
sc <- AddModuleScore(sc, features = list(intersect(sig_neg, rownames(sc))),
                     name = "MHCneg_sig", seed = SEED)

sig_labels <- c(MHCpos_sig1 = "MHC II+ signature", MHCneg_sig1 = "MHC II- signature")
score_df <- FetchData(sc, vars = c("MK_subset", "MHCpos_sig1", "MHCneg_sig1")) %>%
  pivot_longer(c("MHCpos_sig1", "MHCneg_sig1"),
               names_to = "signature", values_to = "score") %>%
  mutate(signature = sig_labels[signature])

vln <- ggplot(score_df, aes(MK_subset, score, fill = MK_subset)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.9) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  facet_wrap(~signature, scales = "free_y") +
  scale_fill_manual(values = ann_colors$Subset, guide = "none") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = "Module score",
       title = "Bulk MHC signatures scored across Sun et al. subsets")
gg_out(vln, "FigS5E_signature_scores_violin", width = 9, height = 4.5)

if ("umap" %in% Reductions(sc)) {
  fp <- FeaturePlot(sc, features = c("MHCpos_sig1", "MHCneg_sig1"), order = TRUE,
                    min.cutoff = "q05", max.cutoff = "q95") &
    scale_color_viridis_c()
  gg_out(fp, "FigS5F_signature_scores_umap", width = 10, height = 4.5)
}

###############################################################################
## 8. VALIDATION CONTROL
## Each subset pseudo-replicate is mapped with the same correlation step used
## for the bulk samples. A specific mapping returns every replicate to its
## own subset.
###############################################################################
ctrl_cor <- cor(pb_lcpm[hvg, ], pb_subset_lcpm[hvg, ], method = "spearman")
ann_ctrl <- data.frame(Subset = pb_meta$subset, row.names = rownames(ctrl_cor))

heatmap_out(ctrl_cor[order(pb_meta$subset), ], "FigS5D_mapping_validation_control",
            width = 6.5, height = 5.5,
            cluster_rows = FALSE, cluster_cols = TRUE,
            annotation_row = ann_ctrl, annotation_colors = ann_colors,
            display_numbers = TRUE, number_format = "%.2f",
            color = colorRampPalette(brewer.pal(9, "Greens"))(100),
            main = "Control: subset pseudo-replicates vs subsets")

ctrl_assign <- data.frame(pseudobulk  = rownames(ctrl_cor),
                          true_subset = pb_meta$subset,
                          assigned    = colnames(ctrl_cor)[max.col(ctrl_cor)])
ctrl_assign$correct <- ctrl_assign$true_subset == ctrl_assign$assigned
write.csv(ctrl_assign, file.path(OUTDIR, "mapping_validation_control.csv"),
          row.names = FALSE)
message(sprintf("Validation control: %d of %d pseudo-replicates map to their own subset",
                sum(ctrl_assign$correct), nrow(ctrl_assign)))

###############################################################################
## 9. SESSION INFO
###############################################################################
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
message("Done. Outputs in: ", normalizePath(OUTDIR))
