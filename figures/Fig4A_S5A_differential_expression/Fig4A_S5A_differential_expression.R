###############################################################################
## Fig 4A and Fig S5A: differential expression, MHC II+ vs MHC II- MKs
##
##   Fig 4A   volcano plot (DESeq2 Wald test, |log2FC| > 1, adjusted P < 0.05)
##   Fig S5A  heatmap of the top 50 genes by adjusted P (rlog, row z-score)
##
## Also writes the DE tables used by the other analyses, including
## DE_genes_sig.tsv, the input for the MHC II+/- signatures in Fig 4C-E/S5D-F.
##
## Run from the repository root:
##   Rscript figures/Fig4A_S5A_differential_expression/Fig4A_S5A_differential_expression.R
###############################################################################

source("R/bulk_rnaseq_helpers.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(viridisLite)
})

OUTDIR      <- results_dir("Fig4A_S5A_differential_expression")
N_LABEL     <- 10   # top genes labelled per direction on the volcano
N_HEATMAP   <- 50   # genes in the Fig S5A heatmap

dds <- build_bulk_dds()
res <- bulk_results(dds)

## ---------------------------------------------------------------------------
## DE tables
## ---------------------------------------------------------------------------
write.table(res, file.path(OUTDIR, "DE_all_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sig <- res[res$sig != "NS", c("baseMean", "log2FoldChange", "lfcSE", "stat",
                              "pvalue", "padj", "gene", "sig")]
write.table(sig, file.path(OUTDIR, "DE_genes_sig.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("DE genes: %d up in MHC II+, %d up in MHC II-",
                sum(sig$sig == "Up_MHCpos"), sum(sig$sig == "Up_MHCneg")))

## ---------------------------------------------------------------------------
## Fig 4A: volcano plot
## ---------------------------------------------------------------------------
plot_df <- res[!is.na(res$padj), ]
plot_df$sig <- factor(plot_df$sig, levels = c("NS", "Up_MHCneg", "Up_MHCpos"))
labels <- plot_df %>%
  filter(sig != "NS") %>%
  group_by(sig) %>%
  slice_max(order_by = -log10(padj), n = N_LABEL, with_ties = FALSE) %>%
  ungroup()

volcano <- ggplot(plot_df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(fill = sig), shape = 21, colour = "black", stroke = 0.3,
             size = 2.5, alpha = 0.9) +
  geom_text_repel(data = labels, aes(label = gene, colour = sig), size = 3,
                  max.overlaps = 20, show.legend = FALSE) +
  geom_vline(xintercept = c(-DE_LOG2FC, DE_LOG2FC), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(DE_PADJ), linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = c(NS = "grey80", Up_MHCneg = "orange", Up_MHCpos = "dodgerblue"),
                    labels = c(NS = "Not significant",
                               Up_MHCneg = "Upregulated in MHC II (-)",
                               Up_MHCpos = "Upregulated in MHC II (+)"),
                    name = NULL) +
  scale_colour_manual(values = c(Up_MHCneg = "darkorange3", Up_MHCpos = "dodgerblue4")) +
  labs(title = "MHC II- vs MHC II+ megakaryocytes",
       x = "log2 fold change (MHC II+ / MHC II-)", y = "-log10 (FDR)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")
ggsave(file.path(OUTDIR, "Fig4A_volcano.pdf"), volcano, width = 7, height = 6)

## ---------------------------------------------------------------------------
## Fig S5A: top 50 genes by adjusted P
## ---------------------------------------------------------------------------
rld  <- rlog(dds, blind = FALSE)
top  <- head(res$gene, N_HEATMAP)
mat  <- assay(rld)[top, ]
mat  <- t(scale(t(mat)))

ann  <- data.frame(Group = colData(dds)$Group, row.names = colnames(mat))
ann_colors <- list(Group = c(MHCneg = "#00BFC4", MHCpos = "#F8766D"))

pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = ann, annotation_colors = ann_colors,
         color = viridis(100), border_color = "black",
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 8, fontsize_col = 10,
         main = "Top 50 differentially expressed genes",
         filename = file.path(OUTDIR, "FigS5A_top50_DE_heatmap.pdf"),
         width = 6, height = 9)

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
message("Done. Outputs in: ", normalizePath(OUTDIR))
