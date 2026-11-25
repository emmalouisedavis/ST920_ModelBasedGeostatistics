# ST920 | Childhood malaria in The Gambia
# 01_dependencies.R

# One-off installation (uncomment only if needed):
# install.packages(c("geoR", "dplyr", "ggplot2", "patchwork", "sf",
#                    "scales", "posterior", "bayesplot"))
# install.packages("cmdstanr", repos = c(
#   "https://stan-dev.r-universe.dev", getOption("repos")))
# cmdstanr::check_cmdstan_toolchain()
# cmdstanr::install_cmdstan(cores = 1)

# For Windows:
# install.packages("INLA", repos = c(getOption("repos"),
#   INLA = "https://inla.r-inla-download.org/R/stable"), dependencies = TRUE)

# For Mac OS:
# install.packages(
#  "INLA",
#  repos = c(
#    CRAN = "https://cloud.r-project.org",
#    INLA = "https://inla.r-inla-download.org/R/stable"
#  ),
#  type = "binary",
#  dependencies = NA
#)

library(geoR)
library(dplyr)
library(ggplot2)
library(patchwork)
library(sf)
library(scales)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(INLA)
