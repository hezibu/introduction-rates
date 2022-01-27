# Function for installing pacman package if it is uninstalled
is_inst <- function(pkg) {
  nzchar(system.file(package = pkg))
}

if (!is_inst("pacman")) {install.packages("pacman")}
library("pacman")

# Load libraries
p_load(tidyverse,
       pbapply,
       optimx,
       purrr,
       matrixStats,
       pracma,
       magrittr,
       doParallel,
       RcppArmadillo,
       GenSA,
       cowplot,
       ucminf,
       HelpersMG)
# sacII


# Load s&c package:
install.packages("snc2022_1.0.tar.gz", repos = NULL, type="source")
#devtools::install_github("zdk123/SpiecEasi")