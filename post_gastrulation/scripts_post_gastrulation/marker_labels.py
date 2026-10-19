import pandas as pd
fine={0:'Mesenchyme (cycling)',1:'Lateral plate / splanchnic mesoderm',2:'Mesenchyme (craniofacial, CRABP1+)',3:'Mesenchyme (cycling)',4:'Mesenchyme / fibroblast (DCN+)',5:'Endothelium',
6:'Chondrogenic / sclerotome mesenchyme',7:'W3 posterior mesoderm / primitive streak (CDX2+)',8:'Mixed neural-mesenchyme (likely doublets)',9:'W3 early mesoderm (HAND1+, cardiac-LPM)',
10:'Neural progenitor (cycling)',11:'Neural progenitor (cycling)',12:'Neural progenitor',13:'Somite / dermomyotome (MEOX1+)',14:'Dorsal neural tube / roof plate (ZIC1+, MSX1+)',
15:'Late radial glia / glial progenitor (BCAN+)',16:'Neural progenitor (HES5+)',17:'Ventral forebrain progenitor (NKX2-1+, OLIG2+)',18:'GABAergic neurons (DLX+)',19:'Neurogenic progenitor (NEUROG1+, HES6+)',
20:'Cortical neurons + RG mix (NEUROD6+)',21:'Intermediate progenitor (EOMES+, NEUROD1+)',22:'Immune (macrophage/microglia) + endothelial',23:'Cortical radial glia (PAX6+, EMX2+, FOXG1+)',
24:'Eye (RPE / optic, PMEL+)',25:'Neurons (NEFM+, early)',26:'Epithelium (endoderm / surface ectoderm)',27:'Neural crest / Schwann precursor (SOX10+)',28:'Placodal / olfactory neurons',
29:'Cortical excitatory neurons (TBR1+)',30:'Cranial / sensory ganglion neurons (ISL1+, SIX1+)',31:'W3 epithelium (ectoderm/endoderm)',32:'Myoblast / muscle (MYF5+)',33:'Immune (macrophage)',
34:'Peripheral sensory neurons (NEFM+, ISL1+)',35:'Erythroid',36:'Doublet (erythroid-neuron)',37:'Neurons (GAD2+, LHX5+)',38:'Neurons (TFAP2A+)',39:'Diencephalic neurons (TCF7L2+, LHX9+)'}
coarse={0:'Mesenchyme',1:'Lateral plate mesoderm',2:'Mesenchyme',3:'Mesenchyme',4:'Mesenchyme',5:'Endothelium',6:'Mesenchyme',7:'Early mesoderm (W3)',8:'Doublet/mixed',9:'Early mesoderm (W3)',
10:'Neural progenitor',11:'Neural progenitor',12:'Neural progenitor',13:'Somite / muscle',14:'Neural progenitor',15:'Neural progenitor',16:'Neural progenitor',17:'Neural progenitor',18:'Neuron',19:'Neural progenitor',
20:'Neuron',21:'Neural progenitor',22:'Immune',23:'Neural progenitor',24:'Eye',25:'Neuron',26:'Epithelium',27:'Neural crest / Schwann',28:'Neuron',29:'Neuron',30:'Neuron',31:'Epithelium',32:'Somite / muscle',33:'Immune',
34:'Neuron',35:'Erythroid',36:'Doublet/mixed',37:'Neuron',38:'Neuron',39:'Neuron'}
cl=pd.read_parquet('work/zeng_leiden.parquet').leiden.astype(int)
pd.DataFrame({'leiden':cl,'marker_fine':cl.map(fine),'marker_coarse':cl.map(coarse)}).to_parquet('work/zeng_marker_labels.parquet')
