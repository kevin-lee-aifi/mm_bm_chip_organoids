hdir = "/home/workspace/"
wdir = hdir + "/TissDiss/EXP-01244"
degdir = wdir + "/EXP_01244_DEGs/"
pltdir = wdir + "/EXP-01244_plots/"
srldir = wdir + "/EXP-01244_serial_obj/"

import sys
sys.path.insert(0, hdir + "TissDiss")

from tissdiss_py_util import *

adata = sc.read_h5ad(srldir + "processed_adata.h5ad")
adata = adata[adata.obs['harmony'] == 'no harmony']

fig, ax = plt.subplots(1, 2, figsize=(18, 8))

sc.pl.umap(adata, color="sample_type", size=1.5, ax=ax[0], show=False, legend_loc='best')

sc.pl.umap(adata, color="leiden", size=1.5, ax=ax[1], show=False, legend_loc='on data')

plt.savefig(pltdir + 'UMAP_plots.pdf', format='pdf', bbox_inches='tight')