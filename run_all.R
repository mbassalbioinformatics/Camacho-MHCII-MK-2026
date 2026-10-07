###############################################################################
## Run every analysis in this repository, in dependency order.
##
## From the repository root:  Rscript run_all.R
##
## 1. Any missing R packages are installed first (see R/required_packages.R).
## 2. The differential expression script runs first, because it writes the DE
##    table used by the Fig 4C-E / S5D-F integration.
## Each script runs in its own R process, exactly as it would on its own.
###############################################################################

source("R/required_packages.R")
todo <- missing_packages()
if (length(todo)) {
  message("Installing missing packages first: ", paste(todo, collapse = ", "))
  install_missing_packages()
}

scripts <- c(
  "figures/Fig4A_S5A_differential_expression/Fig4A_S5A_differential_expression.R",
  "figures/FigS5B_go_enrichment/FigS5B_go_enrichment.R",
  "figures/Fig4B_S5C_tf_activity/Fig4B_S5C_tf_activity.R",
  "figures/Fig4C-E_S5D-F_sun_atlas_integration/Fig4C-E_S5D-F_sun_atlas_integration.R"
)

rscript <- file.path(R.home("bin"), "Rscript")
for (s in scripts) {
  message("\n==== ", s, " ====")
  status <- system2(rscript, shQuote(s))
  if (status != 0) stop("Stopped: ", s, " exited with status ", status)
}
message("\nAll analyses finished. Outputs are in results/.")
