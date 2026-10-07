###############################################################################
## Fig 4B and Fig S5C: transcription factor activity, MHC II+ vs MHC II- MKs
##
##   Fig 4B   top 25 TFs ranked by absolute difference in mean activity
##            (delta activity = mean MHC II+ minus mean MHC II-)
##   Fig S5C  regulatory networks for the top 5 TFs in each direction
##
## TF activity is inferred per sample with VIPER (decoupleR::run_viper) on
## variance-stabilized expression (DESeq2 vst, blind = TRUE), using mouse
## DoRothEA regulons of confidence A to C (minimum regulon size 5).
##
## Run from the repository root:
##   Rscript figures/Fig4B_S5C_tf_activity/Fig4B_S5C_tf_activity.R
###############################################################################

source("R/bulk_rnaseq_helpers.R")
suppressPackageStartupMessages({
  library(decoupleR)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(igraph)
  library(ggraph)
  library(patchwork)
})

OUTDIR        <- results_dir("Fig4B_S5C_tf_activity")
N_TOP_TFS     <- 25    # bars in Fig 4B
N_NETWORK_TFS <- 5     # TFs per direction in Fig S5C
N_TARGETS     <- 10    # targets per TF network, chosen by adjusted P
SEED          <- 13579 # network layouts only

dds    <- build_bulk_dds()
res_df <- bulk_results(dds)

## ---------------------------------------------------------------------------
## TF activity per sample
## ---------------------------------------------------------------------------
expr_mat <- assay(vst(dds, blind = TRUE))

## dorothea is not attached: it has its own deprecated run_viper() that would
## mask decoupleR::run_viper(). The regulons are loaded without attaching it.
data(dorothea_mm, package = "dorothea")
regulon <- dorothea_mm %>% dplyr::filter(confidence %in% c("A", "B", "C"))
message("DoRothEA regulons (A-C): ", length(unique(regulon$tf)), " TFs")

tf_activity <- decoupleR::run_viper(mat = expr_mat, network = regulon,
                                    .source = "tf", .target = "target", .mor = "mor",
                                    minsize = 5)
write.table(tf_activity, file.path(OUTDIR, "TF_activity_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

groups <- data.frame(condition = colnames(dds), Group = as.character(colData(dds)$Group))
tf_summary <- tf_activity %>%
  left_join(groups, by = "condition") %>%
  group_by(source, Group) %>%
  summarise(mean_activity = mean(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = mean_activity) %>%
  mutate(diff = MHCpos - MHCneg,
         Up_in = ifelse(diff > 0, "MHCpos", "MHCneg")) %>%
  arrange(desc(abs(diff)))
write.table(tf_summary, file.path(OUTDIR, "TF_activity_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---------------------------------------------------------------------------
## Fig 4B: top TFs by absolute delta activity
## ---------------------------------------------------------------------------
top_tfs <- head(tf_summary, N_TOP_TFS)
bar <- ggplot(top_tfs, aes(x = reorder(source, diff), y = diff, fill = Up_in)) +
  geom_col(colour = "black", width = 0.7) +
  geom_text(aes(label = round(diff, 2), hjust = ifelse(diff > 0, -0.2, 1.2)), size = 3) +
  coord_flip() +
  scale_fill_manual(values = GROUP_COLORS, labels = GROUP_LABELS, name = "Higher in") +
  ## plotmath Delta renders on every PDF device (a unicode string fails on macOS pdf())
  labs(x = "Transcription factor", y = expression(Delta * " activity (MHC II+ minus MHC II-)")) +
  ## extra room on both sides so the value labels are not clipped
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(colour = "black"))
ggsave(file.path(OUTDIR, "Fig4B_TF_activity_barplot.pdf"), bar, width = 7, height = 7)

## ---------------------------------------------------------------------------
## Fig S5C: TF regulatory networks
## ---------------------------------------------------------------------------
top_pos <- tf_summary %>% dplyr::filter(diff > 0) %>% arrange(desc(diff)) %>% slice_head(n = N_NETWORK_TFS) %>% pull(source)
top_neg <- tf_summary %>% dplyr::filter(diff < 0) %>% arrange(diff)       %>% slice_head(n = N_NETWORK_TFS) %>% pull(source)
log2fc  <- deframe(res_df %>% dplyr::select(gene, log2FoldChange))

tf_network_plot <- function(tf_name, group_label) {
  targets <- regulon %>% dplyr::filter(tf == tf_name, target %in% rownames(expr_mat)) %>% dplyr::select(tf, target)
  keep <- res_df %>% dplyr::filter(gene %in% targets$target) %>% arrange(padj) %>%
    slice_head(n = N_TARGETS) %>% pull(gene)
  targets <- targets %>% dplyr::filter(target %in% keep)
  if (nrow(targets) == 0) return(NULL)

  g <- graph_from_data_frame(targets, directed = TRUE)
  V(g)$type   <- ifelse(V(g)$name == tf_name, "TF", "Target")
  V(g)$log2FC <- ifelse(V(g)$type == "Target", log2fc[V(g)$name], 0)
  V(g)$size   <- ifelse(V(g)$type == "TF", 8, 5)

  set.seed(SEED)
  ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.4, colour = "grey70") +
    geom_node_point(aes(fill = log2FC, size = size), shape = 21, colour = "black") +
    scale_fill_gradient2(low = "orange", mid = "white", high = "dodgerblue",
                         midpoint = 0, name = "log2FC") +
    geom_node_text(aes(label = name), repel = TRUE, size = 4, fontface = "bold") +
    labs(title = paste0(tf_name, " (", group_label, ")")) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

save_networks <- function(tfs, group_label, file) {
  plots <- Filter(Negate(is.null), lapply(tfs, tf_network_plot, group_label = group_label))
  ggsave(file.path(OUTDIR, file), wrap_plots(plots, ncol = 5), width = 22, height = 5)
}
save_networks(top_pos, "MHC II+", "FigS5C_TF_networks_MHCpos.pdf")
save_networks(top_neg, "MHC II-", "FigS5C_TF_networks_MHCneg.pdf")
message("Network TFs, MHC II+: ", paste(top_pos, collapse = ", "))
message("Network TFs, MHC II-: ", paste(top_neg, collapse = ", "))

writeLines(c(capture.output(sessionInfo()), "",
             paste("DoRothEA regulons: dorothea", packageVersion("dorothea"))),
           file.path(OUTDIR, "sessionInfo.txt"))
message("Done. Outputs in: ", normalizePath(OUTDIR))
