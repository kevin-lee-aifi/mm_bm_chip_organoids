import sys
sys.path.insert(0, '/home/jupyter/TissDiss')

from tissdiss_py_util import *

hdir = '/home/jupyter'
wdir = hdir + "/TissDiss/EXP-01244"
pltdir = wdir + "/EXP-01244_plots/"
srldir = wdir + "/EXP-01244_serial_obj/"

adata = sc.read_h5ad(srldir + 'raw_adata.h5ad')

def functions(adata):
    sc.pp.neighbors(adata)
    sc.tl.draw_graph(adata)
    sc.tl.umap(adata)
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

    return adata

def process_adata(adata):

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ig"] = adata.var_names.str.contains("^IG")

    adata = adata[:, (adata.var['mt'] == False) & (adata.var['ig'] == False)]

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

    hm_adata.obsm['hm_X_pca'] = harmony_out.Z_corr.T

    adatas = {
        'adata': functions(adata),
        'hm_adata': functions(hm_adata)
    }

    adatas['adata'].obs['harmony'] = 'no harmony'
    adatas['hm_adata'].obs['harmony'] = 'harmony'

    full_adata = ad.concat(adatas, join = 'outer', merge = 'same')

    return full_adata

adata = process_adata(adata)

sc.pl.umap(adata, color='harmony', legend_loc='on data')

plt.savefig(pltdir + 'harmony_test_umap.pdf', format='pdf', bbox_inches='tight')

print('Done')