#!/usr/bin/bash

PY_ENV="/home/workspace/environment/mm_py_env"

conda create -y -p $PY_ENV -c conda-forge \
    python=3.9 \
    ipykernel \
    h5py \
    scanpy \
    numpy==1.24.4 \
    pandas anndata \
    matplotlib \
    seaborn \
    dill

conda activate $PY_ENV

pip install scrublet \
    igraph \
    fa2-modified \
    harmonypy \
    scib[main] \
    celltypist

python -m ipykernel install --user --name=tissdiss_py_env --display-name="Python (mm_py_env)"