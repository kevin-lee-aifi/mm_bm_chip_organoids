import sys
sys.path.insert(0, '/home/jupyter/TissDiss')

from tissdiss_py_util import *

hdir = '/home/jupyter'
wdir = hdir + "/TissDiss/EXP-01244"
srldir = wdir + "/EXP-01244_serial_obj/"

EXP01244_cr_outs_path = os.path.join(hdir, "TissDiss/EXP-01244/EXP-01244_cr_outs")

adata = ingest_mplex_gex_h5_list(EXP01244_cr_outs_path)

sample_dict = {
    'week2': 'Week 2',
    'week3': 'Week 3',
    'week4': 'Week 4',
    'bm': 'BMMC Start Sample',
    'msc': 'MSC Start Sample'
}

adata.obs['names'] = adata.obs['sample_type'].replace(sample_dict)

adata.write(srldir + 'raw_adata.h5ad', compression='gzip')