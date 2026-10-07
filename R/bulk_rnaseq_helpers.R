###############################################################################
## Shared helpers for the bulk RNA-seq analyses (MHC II+ vs MHC II- MKs)
##
## build_bulk_dds() reproduces the DESeq2 model behind Fig 4A, Fig 4B and
## Fig S5A to S5C:
##   - gene symbols taken from "ENSMUSG..._Symbol" row names
##   - immunoglobulin variable/alpha genes removed (^Ighv, ^Igha, ^Igkv)
##   - duplicate symbols summed, counts converted to integer
##   - design ~ Group, reference level MHCneg (log2FC > 0 = higher in MHC II+)
##   - DESeq2 Wald test, no log2 fold-change shrinkage
## Sourced by the figure scripts. Run every script from the repository root.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
})

BULK_COUNTS_DEFAULT <- "data/counts_matrix.tsv"
DE_PADJ   <- 0.05   # adjusted P cut-off for differential expression
DE_LOG2FC <- 1      # absolute log2 fold-change cut-off

GROUP_COLORS <- c(MHCpos = "dodgerblue", MHCneg = "orange")
GROUP_LABELS <- c(MHCpos = "MHC II (+)", MHCneg = "MHC II (-)")

## Read the RSEM expected-counts matrix and return integer counts by symbol.
## Columns 1 to 4 are MHC II+ and columns 5 to 8 are MHC II-.
read_bulk_counts <- function(path = Sys.getenv("MK_BULK_COUNTS", BULK_COUNTS_DEFAULT)) {
  if (!file.exists(path)) stop("Input not found: ", path, "\nSee data/README.md.")
  counts <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  stopifnot(ncol(counts) == 8)
  colnames(counts) <- c(paste0("MHCpos_", 1:4), paste0("MHCneg_", 1:4))

  symbols <- sub(".*_(.*)$", "\\1", rownames(counts))
  keep    <- !grepl("^Ighv|^Igha|^Igkv", symbols)
  mat     <- rowsum(as.matrix(counts[keep, ]), group = symbols[keep])
  mode(mat) <- "integer"
  mat
}

## Build and fit the DESeq2 model.
build_bulk_dds <- function(counts = read_bulk_counts()) {
  coldata <- data.frame(
    Group = factor(ifelse(grepl("^MHCpos", colnames(counts)), "MHCpos", "MHCneg"),
                   levels = c("MHCneg", "MHCpos")),
    row.names = colnames(counts))
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ Group)
  DESeq(dds, quiet = TRUE)
}

## DESeq2 results as a data frame, ordered by adjusted P, with a DE label.
bulk_results <- function(dds) {
  res <- results(dds)
  res <- res[order(res$padj), ]
  as.data.frame(res) %>%
    mutate(gene = rownames(res),
           sig  = case_when(
             padj < DE_PADJ & log2FoldChange >  DE_LOG2FC ~ "Up_MHCpos",
             padj < DE_PADJ & log2FoldChange < -DE_LOG2FC ~ "Up_MHCneg",
             TRUE ~ "NS"))
}

## Create (if needed) and return a figure's results folder.
results_dir <- function(name) {
  d <- Sys.getenv("MK_RESULTS_ROOT", "results")
  d <- file.path(d, name)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}
