###############################################################################
## Install every package the analyses need. Packages already installed are
## skipped. From the repository root:  Rscript setup.R
##
## run_all.R does this automatically, so running setup.R first is optional.
###############################################################################

source("R/required_packages.R")

todo <- missing_packages()
if (length(todo) == 0) {
  message("All required packages are already installed.")
} else {
  message("Installing: ", paste(todo, collapse = ", "))
  install_missing_packages()
  message("Done. All required packages are installed.")
}
