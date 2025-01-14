import numpy as np
import pandas as pd
import scanpy as sc
import h5py as h5
import anndata as ad
import scrublet
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from joblib import dump, load
from pyclustree import clustree
import harmonypy as hm
import os
import scib
import dill
import celltypist
from celltypist import models

def ingest_mplex_gex_h5_list(cr_outs_path):
    """
    Process all 10x Genomics single-cell RNA sequencing H5 files within a Cell Ranger output directory,
    identify doublets using Scrublet, combine all samples, and filter out predicted doublets.
    
    Parameters
    ----------
    cr_outs_path : str
       Path to the Cell Ranger output directory
       
    Returns
    -------
    concat : AnnData
       Combined and filtered AnnData object containing all samples with doublets removed
       Contains:
       - Raw count matrix (adata.X)
       - Cell metadata (adata.obs) including:
           - sample: Sample identifier
           - Doublet_Score: Scrublet doublet scores
           - Predicted_Doublet: Boolean indicating doublet prediction
       - Gene metadata (adata.var)
    
    Process
    -------
    1. Find all H5 files recursively in the Cell Ranger output directory
    2. For each H5 file:
       - Extract sample name from path
       - Load data into AnnData object
       - Ensure gene names are unique
       - Run Scrublet for doublet detection
       - Add doublet scores to cell metadata
    3. Combine all samples into one AnnData object
    """

    # Dictionary mapping user-friendly sample names to their corresponding IDs
    sample_dict = {
       'week2': "OR07965-01",     # Maps time point labels to sample IDs
       'week3': "OR07965-02", 
       'week4': "OR00001",
       'bm': "BMC07965-007",      # Bone marrow sample
       'msc': "CELL00911"         # Mesenchymal stem cell sample
    }
    
    # Create reverse mapping from sample IDs to their user-friendly names
    id_to_sample = {v: k for k, v in sample_dict.items()}
    
    # Find all filtered_feature_bc_matrix.h5 files in the directory structure
    h5_paths = [os.path.join(root, 'sample_filtered_feature_bc_matrix.h5') 
               for root, _, files in os.walk(cr_outs_path) 
               if 'sample_filtered_feature_bc_matrix.h5' in files]
    
    # Dictionary to store AnnData objects for each sample
    adatas = {}
    
    # Process each H5 file
    for path in h5_paths:
        # Extract sample name from path (e.g., 'BMC07965-007_3')
        name = path.split('per_sample_outs/')[1].split('/')[0]
        
        # Read the H5 file and create AnnData object
        adata = sc.read_10x_h5(path)
        adata.var_names_make_unique()
        adatas[name] = adata
        
        # Run Scrublet for doublet detection
        scrub = scrublet.Scrublet(adata.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)

        # Add Scrublet results to cell metadata
        adata.obs['Doublet_Score'] = doublet_scores
        adata.obs['Predicted_Doublet'] = predicted_doublets

    # Combine all samples into one AnnData object
    concat = ad.concat(adatas, label='sample', join='outer', merge='same')

    # Add metadata column for batched replicates
    concat.obs['base_sample'] = concat.obs['sample'].str.replace(r'_\d+$', '', regex=True)

    # Add sample names (week2, week3, etc.)
    concat.obs['sample_type'] = concat.obs['base_sample'].replace(id_to_sample)
   
    return concat



def process_adata(adata):
    """
    Process an AnnData object through a standard single-cell RNA sequencing workflow,
    including QC metric calculation, normalization, dimensionality reduction, and clustering.
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData object containing raw count data
        Must have:
        - Raw counts in .X
        - 'sample' column in .obs for batch correction
       
    Returns
    -------
    adata : AnnData
        Processed AnnData object containing:
        - QC metrics in .obs (mt/ribo/hb percentages)
        - Raw counts in .layers['counts']
        - Normalized log counts in .X
        - PCA in .obsm['X_pca']
        - UMAP in .obsm['X_umap']
        - Leiden clusters in .obs['leiden']
        - Differentially expressed genes in .uns['rank_genes_groups']
    
    Process
    -------
    1. Calculate QC metrics:
        - Mitochondrial gene percentage (MT-)
        - Ribosomal gene percentage (RPS/RPL)
        - Hemoglobin gene percentage (HB)
    2. Store raw counts and normalize:
        - Save raw counts to .layers['counts']
        - Normalize to library size
        - Log transform
    3. Dimensionality reduction:
        - Select highly variable genes (top 2000)
        - Perform PCA
        - Calculate nearest neighbors
        - Generate UMAP
    4. Clustering and differential expression:
        - Perform Leiden clustering
        - Calculate marker genes for each cluster
    """
    
    # Calculate gene QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-")  # Mitochondrial genes
    adata.var["ig"] = adata.var_names.str.contains("^IG")  # IG genes

    adata = adata[:, (adata.var['mt'] == False) & (adata.var['ig'] == False)]  # Filtering out ig and mt genes
    
    # Calculate cell QC metrics
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=["mt", "ig"],
        inplace=True,
        log1p=True
    )
    
    # Store raw counts and normalize data
    adata.layers["counts"] = adata.X.copy()  # Save raw counts
    sc.pp.normalize_total(adata)  # Normalize to library size
    sc.pp.log1p(adata)  # Log transform
    
    # Select variable genes and reduce dimensions
    sc.pp.highly_variable_genes(
        adata, 
        n_top_genes=2000, 
        batch_key="sample"  # Account for batch effects
    )
    sc.tl.pca(adata)  # Perform PCA
    sc.pp.neighbors(adata)  # Calculate nearest neighbors
    sc.tl.draw_graph(adata)
    sc.tl.umap(adata)  # Generate UMAP embedding
    
    # Cluster cells and find markers
    sc.tl.leiden(
        adata, 
        flavor="igraph", 
        n_iterations=2,
        resolution=1.1
    )

    sc.tl.paga(adata, groups='leiden')
    
    sc.tl.rank_genes_groups(
        adata, 
        groupby="leiden", 
        method="wilcoxon"
    )

    sc.tl.diffmap(adata)
    
    return adata