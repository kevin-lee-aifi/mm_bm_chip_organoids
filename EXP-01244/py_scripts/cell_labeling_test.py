hdir = "/home/workspace/"
wdir = hdir + "/TissDiss/EXP-01244"
degdir = wdir + "/EXP_01244_DEGs/"
srldir = wdir + "/EXP-01244_serial_obj/"

import sys
sys.path.insert(0, hdir + "TissDiss")

from tissdiss_py_util import *

adata = sc.read_h5ad(srldir + 'processed_adata.h5ad')
adata = adata[adata.obs['harmony'] == 'no harmony']

