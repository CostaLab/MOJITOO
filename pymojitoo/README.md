# MOJITOO<img src="../inst/figures/LOGO.png" align="right" width="120" />
Mingbo Cheng<sup>1</sup>,
Zhijian Li<sup>1</sup>,
Ivan G. Costa<sup>1*</sup>


<sup>1</sup>Institute for Computational Genomics, Faculty of Medicine, RWTH Aachen University, Aachen, 52074 Germany

### Installation
```{shell}
git clone https://github.com/CostaLab/MOJITOO.git
cd pymojitoo
pip install .
```

### Usage

##### dummy adata

```{python}
## not for real, please use multimodal data instead
import mojitoo
import scanpy
import scanpy as sc
import episcanpy as epi
ds = sc.datasets
adata = ds.pbmc68k_reduced()
epi.tl.tfidf(adata)
epi.tl.lsi(adata)
mojitoo.mojitoo(adata, ["X_pca", "X_lsi"])
sc.pp.neighbors(adata, use_rep='X_mojitoo')
sc.pl.umap(adata, color='louvain')
```


##### mudata

```{python}
import muon
import mojitoo
import mudatasets as mds
import scanpy as sc
import episcanpy as epi

mdata = mds.load('pbmc3k_multiome',files=["filtered_feature_bc_matrix.h5"])

atac = mdata.mod['atac']
rna = mdata.mod['rna']


epi.pp.cal_var(atac)
epi.pp.select_var_feature(atac, nb_features=5000)
epi.tl.tfidf(atac)
epi.tl.lsi(atac, n_components=50)


sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.scale(rna, max_value=10)
sc.tl.pca(rna, svd_solver='arpack')


mdata.obsm["pca"] = rna.obsm["X_pca"]
mdata.obsm["lsi"] = atac.obsm["X_lsi"]


mojitoo.mojitoo(mdata, reduction_list=["pca", "lsi"],  dims_list=(range(50), range(1,50)),reduction_name='mojitoo', overwrite=True)

sc.pp.neighbors(mdata, use_rep='mojitoo')
sc.tl.louvain(mdata, resolution=0.5)
sc.tl.umap(mdata)
sc.pl.embedding(mdata, color='louvain', basis='umap')
```


### References:
MOJITOO: a fast and universal method for integration of multimodal single cell data [link](https://doi.org/10.1093/bioinformatics/btac220)
