# =============================================================================
# Run once (from the RStudio project "EPA with Unknown Clusters"):
# installs the CRAN dependencies and the two project packages.
# Run it in a fresh R session (Session > Restart R): reinstalling a package that is loaded in
# the session corrupts its lazy-load database ("lazy-load database ... is corrupt").
# =============================================================================
if (any(c("PanelKmeansInference", "clusteredEPA") %in% loadedNamespaces()))
  stop("PanelKmeansInference or clusteredEPA is loaded in this session. ",
       "Restart R (Ctrl+Shift+F10 in RStudio) and run this script again before loading them.")
cran <- c("devtools", "roxygen2", "Rcpp", "here", "dplyr", "tidyr", "openxlsx", "mclust", "ggplot2",
          "intervals", "zoo", "xts", "readxl", "glmnet", "xgboost", "foreach", "pbapply",
          "patchwork", "ggrepel",               # figures of the application (R/prepare_data.R)
          "plm", "lmtest", "sandwich",          # cluster regressions (step 6 of run_application.R)
          "e1071", "randomForest")              # forecasts (only if regenerate_forecasts = TRUE)
miss <- setdiff(cran, rownames(installed.packages()))
if (length(miss)) install.packages(miss)
# fbi (FRED-MD tools) is on GitHub; needed only to regenerate the forecasts of the application
if (!requireNamespace("fbi", quietly = TRUE)) try(devtools::install_github("cykbennie/fbi"))

library(here)
devtools::document(here("R", "PanelKmeansInference"))
devtools::install(here("R", "PanelKmeansInference"), upgrade = FALSE)
devtools::document(here("R", "clusteredEPA"))
devtools::install(here("R", "clusteredEPA"), upgrade = FALSE)

# Compile (once) the C++ walk used by the selective tests. It needs a C++ compiler (Rtools on
# Windows, Xcode command line tools on macOS). Without it the packages fall back to an R version
# that gives identical results but is about 10 times slower.
library(PanelKmeansInference)
ok <- !is.null(PanelKmeansInference:::.get_cpp_walk())
message(if (ok) "C++ walk compiled and cached." else "C++ walk NOT available: R fallback will be used.")
ok2 <- !is.null(clusteredEPA:::.get_cpp_gen("rolling_errors_cpp"))
message(if (ok2) "C++ simulation helpers compiled and cached." else
          "C++ simulation helpers NOT available: R fallback will be used.")
