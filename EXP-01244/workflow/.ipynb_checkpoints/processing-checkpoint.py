hdir = "/home/workspace/"
wdir = hdir + "/TissDiss/EXP-01244"
degdir = wdir + "/EXP_01244_DEGs/"
srldir = wdir + "/EXP-01244_serial_obj/"

import sys
sys.path.insert(0, hdir + "TissDiss")

from tissdiss_py_util import *

def ingest_mplex_gex_h5_list(cr_outs_path):

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

    # Downsampling to bone marrow cell count, excluding the week 3 sample
    samples = adata.obs['base_sample'].values.unique().tolist()
    
    ds_adata = {}
    
    for sample in samples:
        if sample == 'OR07965-02':
            mask = adata.obs['base_sample'] == sample
            ds_adata[sample] = adata[mask]
            continue
    
        mask = adata.obs['base_sample'] == sample
        downsampled = sc.pp.subsample(
            adata[mask].copy(),
            n_obs = sum(adata.obs['base_sample'] == 'BMC07965-007'),
            random_state = 0,
            copy = True
        )
    
        ds_adata[sample] = downsampled
    
    adata = ad.concat(ds_adata, label = "base_sample", join = 'outer', merge = 'same')
    
    # Filtering out mitochondrial and Ig genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ig"] = adata.var_names.str.contains("^IG")
    
    adata = adata[:, (adata.var['mt'] == False) & (adata.var['ig'] == False)]

    # Applying scanpy processing pipeline
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=["mt", "ig"],
        inplace=True,
        log1p=True
    )

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata, 
        n_top_genes=2000, 
        batch_key="sample"
    )
    sc.tl.pca(adata)

    hm_adata = adata.copy()

    X = hm_adata.obsm['X_pca'].astype(np.float64)

    harmony_out = hm.run_harmony(X, hm_adata.obs, 'sample_type')

    hm_adata.obsm['X_pca'] = harmony_out.Z_corr.T

    adatas = {
        'adata': adata,
        'hm_adata': hm_adata
    }

    adatas['adata'].obs['harmony'] = 'no harmony'
    adatas['hm_adata'].obs['harmony'] = 'harmony'

    full_adata = ad.concat(adatas, join = 'outer', merge = 'same')

    sc.pp.neighbors(full_adata)
    sc.tl.draw_graph(full_adata)
    sc.tl.umap(full_adata)
    sc.tl.leiden(
        full_adata, 
        flavor="igraph", 
        n_iterations=2,
        resolution=1.1
    )
    sc.tl.paga(full_adata, groups='leiden')
    sc.tl.rank_genes_groups(
        full_adata, 
        groupby="leiden", 
        method="wilcoxon"
    )

    sc.tl.diffmap(full_adata)

    return full_adata[full_adata.obs['harmony'] == 'no harmony']

EXP01244_cr_outs_path = os.path.join(hdir, "TissDiss/EXP-01244/EXP-01244_cr_outs")

adata = ingest_mplex_gex_h5_list(EXP01244_cr_outs_path)

adata = process_adata(adata)

adata.write(srldir + 'processed_adata.h5ad', compression='gzip')