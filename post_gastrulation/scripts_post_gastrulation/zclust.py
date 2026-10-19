import numpy as np,pandas as pd,scipy.sparse as sp,scanpy as sc,anndata as ad
from sklearn.decomposition import PCA
Z=sp.load_npz('work/zeng_sub.npz'); n=Z.shape[0]
mu=np.asarray(Z.mean(0)).ravel(); sd=np.sqrt(np.asarray(Z.multiply(Z).mean(0)).ravel()-mu**2)+1e-6
rng=np.random.default_rng(0); sub=rng.choice(n,120000,replace=False)
def scale(M): return np.clip((M.toarray()-mu)/sd,-10,10).astype(np.float32)
pca=PCA(50,random_state=0).fit(scale(Z[sub]))
P=np.vstack([pca.transform(scale(Z[a:a+30000])) for a in range(0,n,30000)]).astype(np.float32)
np.save('work/zeng_pca.npy',P); np.save('work/zeng_scale.npy',np.vstack([mu,sd]))
c=pd.read_parquet('work/zeng_cells_qc.parquet')
A=ad.AnnData(obs=c.set_index(c.barcode+'_'+c['sample'])); A.obsm['X_pca']=P
import harmonypy
h=harmonypy.run_harmony(P,c,'sample',max_iter_harmony=10); H=np.asarray(h.Z_corr); H=H.T if H.shape[0]!=n else H
A.obsm['X_harmony']=H.astype(np.float32); np.save('work/zeng_harmony.npy',A.obsm['X_harmony'])
sc.pp.neighbors(A,n_neighbors=15,use_rep='X_harmony',random_state=0)
sc.tl.leiden(A,resolution=1.0,flavor='igraph',n_iterations=2,random_state=0,key_added='leiden')
sc.tl.umap(A,random_state=0)
A.obs[['leiden']].assign(umap1=A.obsm['X_umap'][:,0],umap2=A.obsm['X_umap'][:,1]).to_parquet('work/zeng_leiden.parquet')
print('clusters',A.obs.leiden.nunique())
