# =============================================================================
# Run once (from the RStudio project "EPA with Unknown Clusters"):
# installs the CRAN dependencies and the two project packages.
# Run it in a fresh R session (Session > Restart R): reinstalling a package that is loaded in
# the session corrupts its lazy-load database ("lazy-load database ... is corrupt").
# =============================================================================
# if (any(c("PanelKmeansInference", "clusteredEPA") %in% loadedNamespaces()))
#   stop("PanelKmeansInference or clusteredEPA is loaded in this session. ",
#        "Restart R (Ctrl+Shift+F10 in RStudio) and run this script again before loading them.")
# RStudio can load the two packages by itself after a restart (e.g. for autocompletion in open
# scripts). Unload them first (clusteredEPA before PanelKmeansInference, which it imports).
for (p in c("clusteredEPA", "PanelKmeansInference")) {
  if (p %in% loadedNamespaces()) {
    ok <- tryCatch({ unloadNamespace(p); TRUE }, error = function(e) FALSE)
    if (!ok) stop(p, " is loaded and cannot be unloaded. Close the open scripts that use it, ",
                  "restart R, or run in a terminal: Rscript Simulation/00_setup.R")
  }
}
cran <- c("devtools", "roxygen2", "Rcpp", "here", "dplyr", "tidyr", "openxlsx", "mclust", "ggplot2",
          "intervals", "zoo", "xts", "readxl", "glmnet", "xgboost", "foreach", "pbapply",
          "patchwork", "ggrepel",               # figures of the application (R/prepare_data.R)
          "plm", "lmtest", "sandwich",          # cluster regressions (step 6 of run_application.R)
          "e1071", "randomForest")              # forecasts (only if regenerate_forecasts = TRUE)
miss <- setdiff(cran, rownames(installed.packages()))
if (length(miss)) install.packages(miss)
# fbi (FRED-MD tools) is on GitHub; needed only to regenerate the forecasts of the application
if (!requireNamespace("fbi", quietly = TRUE)) try(devtools::install_github("cykbennie/fbi"))

# library(here)
# devtools::document(here("R", "PanelKmeansInference"))
# devtools::install(here("R", "PanelKmeansInference"), upgrade = FALSE)
# devtools::document(here("R", "clusteredEPA"))
# devtools::install(here("R", "clusteredEPA"), upgrade = FALSE)

# PanelKmeansInference is installed first because clusteredEPA depends on it. If a local copy
# of a package is in R/<package> (development: its own git repository, ignored by this one),
# it is documented and installed from there; otherwise it is installed from GitHub at the
# version used in the paper.
library(here)
pkg_versions <- c(PanelKmeansInference = "v0.2.4", clusteredEPA = "v0.2.4")
for (p in names(pkg_versions)) {
  local_dir <- here("R", p)
  if (file.exists(file.path(local_dir, "DESCRIPTION"))) {
    message("Installing ", p, " from the local copy in R/", p)
    devtools::document(local_dir)
    devtools::install(local_dir, upgrade = FALSE)
  } else {
    message("Installing ", p, " ", pkg_versions[[p]], " from GitHub")
    remotes::install_github(paste0("akoguzhan/", p, "@", pkg_versions[[p]]),
                            upgrade = "never", force = TRUE)
  }
}

# Compile (once) the C++ walk used by the selective tests. It needs a C++ compiler (Rtools on
# Windows, Xcode command line tools on macOS). Without it the packages fall back to an R version
# that gives identical results but is about 10 times slower.
library(PanelKmeansInference)
ok <- !is.null(PanelKmeansInference:::.get_cpp_walk())
message(if (ok) "C++ walk compiled and cached." else "C++ walk NOT available: R fallback will be used.")
ok2 <- !is.null(clusteredEPA:::.get_cpp_gen("rolling_errors_cpp"))
message(if (ok2) "C++ simulation helpers compiled and cached." else
          "C++ simulation helpers NOT available: R fallback will be used.")
