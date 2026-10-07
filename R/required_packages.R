###############################################################################
## Every package the analyses in this repository need, in one place.
## Used by setup.R and by the pre-flight check in run_all.R.
###############################################################################

REQUIRED_CRAN <- c("Seurat", "SeuratObject", "dplyr", "tidyr", "tibble", "ggplot2",
                   "ggrepel", "pheatmap", "viridisLite", "RColorBrewer", "gridExtra",
                   "igraph", "ggraph", "patchwork")

REQUIRED_BIOC <- c("DESeq2", "edgeR", "clusterProfiler", "org.Mm.eg.db", "enrichplot",
                   "decoupleR", "dorothea", "viper")

## Names of required packages that are not installed.
missing_packages <- function(pkgs = c(REQUIRED_CRAN, REQUIRED_BIOC)) {
  pkgs[!vapply(pkgs, function(p) suppressWarnings(requireNamespace(p, quietly = TRUE)),
               logical(1))]
}

## Install whatever is missing: CRAN packages first, then Bioconductor.
install_missing_packages <- function() {
  repos <- getOption("repos")
  if (is.null(repos) || is.na(repos["CRAN"]) || repos["CRAN"] == "@CRAN@") {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }

  cran <- intersect(missing_packages(), REQUIRED_CRAN)
  if (length(cran)) install.packages(cran)

  bioc <- intersect(missing_packages(), REQUIRED_BIOC)
  if (length(bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(bioc, update = FALSE, ask = FALSE)
  }

  still <- missing_packages()
  if (length(still)) {
    stop("These packages could not be installed: ", paste(still, collapse = ", "),
         "\nInstall them manually, then run this again.")
  }
  invisible(TRUE)
}
