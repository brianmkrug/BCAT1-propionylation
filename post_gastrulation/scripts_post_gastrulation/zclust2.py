import numpy as np,pandas as pd,scanpy as sc,anndata as ad
H=np.load('work/zeng_harmony.npy'); c=pd.read_parquet('work/zeng_cells_qc.parquet')
A=ad.AnnData(obs=pd.DataFrame(index=np.arange(len(c)).astype(str))); A.obsm['X_harmony']=H[:,:30]
sc.pp.neighbors(A,n_neighbors=15,use_rep='X_harmony',random_state=0)
print('nn done',flush=True)
sc.tl.leiden(A,resolution=1.0,flavor='igraph',n_iterations=2,random_state=0,key_added='leiden')
pd.DataFrame({'leiden':A.obs.leiden.values}).to_parquet('work/zeng_leiden.parquet'); print('clusters',A.obs.leiden.nunique(),flush=True)
import scipy.sparse as sp; sp.save_npz('work/zeng_knn.npz',A.obsp['connectivities'])
