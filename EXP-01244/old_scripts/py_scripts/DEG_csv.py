hdir = "/home/workspace/"
wdir = hdir + "/TissDiss/EXP-01244"
degdir = wdir + "/EXP_01244_DEGs/"
srldir = wdir + "/EXP-01244_serial_obj/"

import sys
sys.path.insert(0, hdir + "TissDiss")

from tissdiss_py_util import *

adata = sc.read_h5ad(srldir + "processed_adata.h5ad")
adata = adata[adata.obs['harmony'] == 'no harmony']

sample_types = adata.obs['sample_type'].unique()

for tissue_type in sample_types:
    mask = adata.obs['sample_type'] == tissue_type
    ranking = adata[mask].uns['rank_genes_groups']
    gene_names = ranking['names']
    gene_scores = ranking['scores']
    gene_pvals_adj = ranking['pvals_adj']
    clusters = gene_names.dtype.names
    
    de_genes = []
    for cluster in clusters:
        names = gene_names[cluster]
        scores = gene_scores[cluster]
        pvals_adj = gene_pvals_adj[cluster]
        for name, score, pval_adj in zip(names, scores, pvals_adj):
            if pval_adj < 0.05:
                de_genes.append((cluster, name, score, pval_adj))
    
    df_de_genes = pd.DataFrame(de_genes, columns=['Cluster', 'Gene', 'Score', 'Adjusted p-value'])
    df_de_genes.to_csv(degdir + f'{tissue_type}_DEGs.csv', index=False)