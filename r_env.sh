#!/usr/bin/bash

R_ENV="/home/workspace/environment/mm_r_env"

conda create -y -p $R_ENV -c conda-forge \
    r-base=4.3 \
    r-essentials \
    r-biocmanager \
    r-irkernel \
    r-matrix \
    r-hdf5r \
    r-ggplot2 \
    r-tidyverse	\
    r-patchwork \
    r-dplyr \
    r-readxl \
    r-seurat \
    r-seuratobject \
    r-remotes \
    r-harmony \
    r-clustree \
    r-rcurl \

conda activate $R_ENV

Rscript -e '
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

remotes::install_github("satijalab/azimuth")

remotes::install_cran("qs", type = "source", configure.args = "--with-simd=AVX2")
'

R -e "IRkernel::installspec(name = 'mm_r_env', displayname = 'R (mm_r_env)')"