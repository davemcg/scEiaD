import sys
import scIB
import scanpy as sc
#adata = sc.read_h5ad('anndata/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__onlyDROPLET__batch__scVI__dims50__preFilter__mindist0.1__nneighbors500__knn10.h5ad')
args = sys.argv

adata = sc.read_h5ad(args[1])
adata
adata_filter  = adata[adata.obs['CellType'] != 'NA']
adata_filterSUB = adata[adata.obs['SubCellType'] != 'NA']
#silh = scIB.metrics.silhouette(adata, 'cluster', embed = 'X_scvi') # 3 hours

#lisi = scIB.metrics.lisi(adata = adata, batch_key = 'batch', label_key = 'CellType') # 8 hours
ari = scIB.metrics.ari(adata_filter, 'CellType', 'cluster')
ari_sub = scIB.metrics.ari(adata_filterSUB, 'SubCellType', 'cluster')

nmi = scIB.metrics.nmi(adata_filter, 'CellType', 'cluster')
nmi_sub = scIB.metrics.nmi(adata_filterSUB, 'SubCellType', 'cluster')

sc.pp.pca(adata, n_comps = int(args[2]))
before = scIB.metrics.pcr(adata, 'batch', embed = 'X_pca', recompute_pca = False)
after = scIB.metrics.pcr(adata, 'batch', embed = 'X_scvi', recompute_pca = False)
pcr = (before-after)/before


f = open(args[3], 'w')

f.write("pcr," + str(pcr) + '\n')
f.write("nmi," + str(nmi) + '\n')
f.write("nmi_sub," + str(nmi_sub) + '\n')

f.write("ari," + str(ari) + '\n')
f.write("ari_sub," + str(ari_sub) + '\n')

f.close()
