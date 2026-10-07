###############################################################################
## Fig S5B: GO biological process enrichment of MHC II+ and MHC II- genes
##
## Genes up in MHC II+ and up in MHC II- (DESeq2, |log2FC| > 1, adjusted
## P < 0.05) are tested separately with clusterProfiler::enrichGO
## (org.Mm.eg.db, ontology BP, Benjamini-Hochberg, P and q < 0.05). The
## background is the clusterProfiler default (all annotated mouse genes).
## The top 20 terms per group are shown as dot plots side by side.
##
## Run from the repository root:
##   Rscript figures/FigS5B_go_enrichment/FigS5B_go_enrichment.R
###############################################################################

source("R/bulk_rnaseq_helpers.R")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(gridExtra)
})

OUTDIR  <- results_dir("FigS5B_go_enrichment")
N_TERMS <- 20

res <- bulk_results(build_bulk_dds())
up_genes   <- res$gene[res$sig == "Up_MHCpos"]
down_genes <- res$gene[res$sig == "Up_MHCneg"]
message(sprintf("Testing %d MHC II+ genes and %d MHC II- genes", length(up_genes), length(down_genes)))

run_go <- function(symbols) {
  ids <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  enrichGO(gene = ids$ENTREZID, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
           ont = "BP", pAdjustMethod = "BH",
           pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
}
go_up   <- run_go(up_genes)
go_down <- run_go(down_genes)

write.table(as.data.frame(go_up),   file.path(OUTDIR, "GO_BP_MHCpos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(go_down), file.path(OUTDIR, "GO_BP_MHCneg.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

style <- theme(axis.text.y = element_text(size = 8),
               axis.text.x = element_text(size = 10),
               plot.title  = element_text(size = 12, face = "bold"))
dp_pos <- dotplot(go_up,   showCategory = N_TERMS, title = "MHC (+) enriched genes") + style
dp_neg <- dotplot(go_down, showCategory = N_TERMS, title = "MHC (-) enriched genes") + style

pdf(file.path(OUTDIR, "FigS5B_GO_BP_dotplot.pdf"), width = 14, height = 8)
grid.arrange(dp_pos, dp_neg, ncol = 2)
invisible(dev.off())

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
message("Done. Outputs in: ", normalizePath(OUTDIR))
