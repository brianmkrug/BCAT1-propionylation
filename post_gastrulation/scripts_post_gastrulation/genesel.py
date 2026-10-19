import pandas as pd,numpy as np
z=pd.read_parquet('work/zeng_genes.parquet'); x=pd.read_parquet('work/xu_genes.parquet')
shared=set(z.ensembl)&set(x.gene_id)
# Zeng HVG: normalized dispersion within 20 mean bins
zz=z[(z['mean']>0.01)&z.ensembl.isin(shared)&~z.symbol.str.match(r'^(MT-|RPS|RPL|HB[ABDEGZQM]\d?$|MALAT1|XIST)')].copy()
zz['disp']=np.log(zz['var']/zz['mean']); zz['bin']=pd.qcut(zz['mean'],20,labels=False)
zz['nd']=zz.groupby('bin')['disp'].transform(lambda d:(d-d.median())/(d.quantile(.75)-d.quantile(.25)+1e-9))
hvg=set(zz.nlargest(2500,'nd').ensembl)
# Xu markers per type
pb=np.load('work/xu_pb_type.npy'); nz=np.load('work/xu_nz_type.npy'); types=pd.read_csv('work/xu_types.csv').iloc[:,0].values
c=pd.read_csv('xu_cells_BCAT.csv.gz'); ncell=(c['developmental system']+'|'+c['final_annotation']).value_counts().reindex(types).values
cpm=pb/pb.sum(0)*1e6; frac=nz/ncell; lc=np.log2(cpm+1)
mk=set()
for j in np.where(ncell>=30)[0]:
    other=np.delete(lc,j,axis=1).mean(1); fc=lc[:,j]-other; fc[frac[:,j]<0.25]=-99
    mk|=set(x.gene_id.values[np.argsort(-fc)[:10]])
mk&=shared
markers='SOX2 PAX6 NES VIM HES5 FABP7 SOX1 SLC1A3 EOMES NEUROG2 NHLH1 STMN2 DCX ELAVL4 TUBB3 GAP43 SYT1 SOX10 FOXD3 TFAP2A PAX3 PDGFRA OLIG2 OLIG1 AQP4 GFAP MPZ PLP1 COL1A1 COL3A1 PDGFRB TWIST1 PRRX1 MEOX1 TCF15 PAX1 MYOD1 MYOG MYF5 DES MYL7 TNNT2 NKX2-5 CDH5 PECAM1 KDR CLDN5 HBA1 HBG1 HBZ GYPA PTPRC LYZ SPI1 C1QA P2RY12 CX3CR1 EPCAM FOXA2 SOX17 CDH1 AFP ALB APOA2 KRT5 TP63 KRT15 KRT8 KRT18 RGS5 TTR CRYAA VSX2 TBXT MIXL1 MESP1 HAND1 FOXF1 OSR1 NANOS3 DPPA3 BCAT1 BCAT2 MYC GAD2 SLC17A6 RELN HOXB4 HOXA2 OTX2 FOXG1 EMX2 EN1 GBX2 NKX2-1 DLX2 HES1 MKI67 TOP2A'.split()
mm=set(z[z.symbol.isin(markers)&z.ensembl.isin(shared)].ensembl)
G=sorted(hvg|mk|mm); print('hvg',len(hvg),'xu markers',len(mk),'manual',len(mm),'union',len(G))
pd.Series(G).to_csv('work/geneset_ensembl.csv',index=False)
