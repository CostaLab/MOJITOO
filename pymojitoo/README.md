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

### dummy example
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

### References:
MOJITOO: a fast and universal method for integration of multimodal single cell data [link](https://doi.org/10.1093/bioinformatics/btac220)
