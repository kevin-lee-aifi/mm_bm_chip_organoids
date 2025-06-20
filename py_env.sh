#!/usr/bin/bash

PY_ENV="/home/workspace/environment/mm_py_env"

# Create conda environment with Python and R packages
conda create -y -p $PY_ENV -c conda-forge \
    python=3.9 \
    ipykernel \
    h5py \
    scanpy \
    numpy==1.24.4 \
    pandas anndata \
    matplotlib \
    seaborn \
    dill \
    r-base r-essentials rpy2

# Activate the environment
conda activate $PY_ENV

# Install Python packages via pip
pip install scrublet \
    igraph \
    fa2-modified \
    harmonypy \
    scib[main] \
    celltypist \
    decoupler \
    anndata2ri \
    gprofiler-official \
    gseapy \
    lxml

# Install Bioconductor packages via R
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("SingleCellExperiment")'

# Install the Jupyter kernel
python -m ipykernel install --user --name=tissdiss_py_env --display-name="Python (mm_py_env)"