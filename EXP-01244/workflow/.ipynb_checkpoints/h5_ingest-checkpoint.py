import sys
sys.path.insert(0, '/home/jupyter/TissDiss')

from tissdiss_py_util import *

hdir = '/home/jupyter'
wdir = hdir + "/TissDiss/EXP-01244"
srldir = wdir + "/EXP-01244_serial_obj/"

EXP01244_cr_outs_path = os.path.join(hdir, "TissDiss/EXP-01244/EXP-01244_cr_outs")

adata = ingest_mplex_gex_h5_list(EXP01244_cr_outs_path)

adata.write(srldir + 'raw_adata.h5ad', compression='gzip')

print('adata saved')