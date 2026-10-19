import numpy as np,pandas as pd,scipy.sparse as sp,gc
from sklearn.decomposition import PCA
def stats(M):
    mu=np.asarray(M.mean(0)).ravel(); sd=np.sqrt(np.asarray(M.multiply(M).mean(0)).ravel()-mu**2)+1e-6; return mu,sd
def proj(M,mu,sd,pca): return np.vstack([pca.transform(np.clip((M[a:a+25000].toarray()-mu)/sd,-10,10)) for a in range(0,M.shape[0],25000)]).astype(np.float32)
X=sp.load_npz('work/xu_sub.npz'); xm,xs=stats(X)
rng=np.random.default_rng(1); sub=rng.choice(X.shape[0],80000,replace=False)
pca=PCA(50,random_state=0).fit(np.clip((X[sub].toarray()-xm)/xs,-10,10))
PX=proj(X,xm,xs,pca); del X; gc.collect()
Z=sp.load_npz('work/zeng_sub.npz'); zm,zs=stats(Z); PZ=proj(Z,zm,zs,pca); del Z; gc.collect()
np.save('work/joint_pca_xu.npy',PX); np.save('work/joint_pca_zeng.npy',PZ)
import harmonypy
meta=pd.DataFrame({'ds':['xu']*len(PX)+['zeng']*len(PZ)})
h=harmonypy.run_harmony(np.vstack([PX,PZ]),meta,'ds',max_iter_harmony=10); H=np.asarray(h.Z_corr); H=H.T if H.shape[0]!=len(meta) else H
HX,HZ=H[:len(PX)].astype(np.float32),H[len(PX):].astype(np.float32); np.save('work/joint_h_xu.npy',HX); np.save('work/joint_h_zeng.npy',HZ)
print('harmony done',flush=True)
import pynndescent
idx=pynndescent.NNDescent(HX,n_neighbors=30,metric='euclidean',random_state=0,n_jobs=2); idx.prepare()
nn=np.vstack([idx.query(HZ[a:a+50000],k=30)[0] for a in range(0,len(HZ),50000)])
np.save('work/zeng_xu_nn.npy',nn.astype(np.int32))
c=pd.read_csv('xu_cells_BCAT.csv.gz')
out={}
for col in ['developmental system','final_annotation']:
    lab=c[col].astype('category'); codes=lab.cat.codes.values[nn]; K=len(lab.cat.categories)
    cnt=np.zeros((len(HZ),K),np.int16)
    for j in range(30): np.add.at(cnt,(np.arange(len(HZ)),codes[:,j]),1)
    best=cnt.argmax(1); out['xu_'+col.replace(' ','_')]=lab.cat.categories[best].values; out['conf_'+col.replace(' ','_')]=cnt.max(1)/30
# self-consistency: leave-one-out on Xu subset
q=rng.choice(len(HX),20000,replace=False); nnq=idx.query(HX[q],k=31)[0][:,1:]
lab=c['developmental system'].values; acc=np.mean([pd.Series(lab[r]).value_counts().index[0]==lab[i] for r,i in zip(nnq,q)]); print('Xu LOO system acc',round(acc,3))
pd.DataFrame(out).to_parquet('work/zeng_transfer.parquet'); print('done')
