# Load libraries

# When working in Linux some packages needed
# some linux libraries. In the terminal type:
# sudo apt-get install build-essential
# 
# sudo apt-get install libcurl4-openssl-dev
# sudo apt-get install libssl-dev
# sudo apt-get install libxml2-dev
# sudo apt install cmake
# To install a package for all users, from the root terminal (puTTY), type:
# sudo R
# and then install.packages(...)

library(tidyverse)
library(pbapply)
library(optimx)
library(purrr)
library(matrixStats)
library(pracma)
library(magrittr)
library(doParallel)
library(RcppArmadillo)
library(GenSA)
library(cowplot)
library(ucminf)
library(HelpersMG)
library(future.apply)
library(progressr)
library(snc2022) # my package


# Load s&c package:
# install.packages("snc2022_1.0.tar.gz", repos = NULL, type="source")