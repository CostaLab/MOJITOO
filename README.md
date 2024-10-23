# MOJITOO<img src="man/figures/logo.png" align="right" width="120" />
Mingbo Cheng<sup>1</sup>,
Zhijian Li<sup>1</sup>,
Ivan G. Costa<sup>1*</sup>


<sup>1</sup>Institute for Computational Genomics, Faculty of Medicine, RWTH Aachen University, Aachen, 52074 Germany

### Installation

##### Default R version
```{r}
install.packages("devtools")
devtools::install_github("CostaLab/MOJITOO", build_vignettes = TRUE)
```

##### Python version alternative

We have implemented a simple python version, please check out the installation detail at [pymojitoo](https://github.com/CostaLab/MOJITOO/tree/main/pymojitoo).

### Howto
The MOJITOO integration using `Seurat` can be found at [seurat](https://costalab.github.io/MOJITOO/articles/SeuratObject_integration.html). `ArchR` example of integration is [ArchR](https://costalab.github.io/MOJITOO/articles/ArchRObject_integration.html). Integration using `R matrix` example can be found [Matrix](https://costalab.github.io/MOJITOO/articles/Matrix_integration.html).


### References:
MOJITOO: a fast and universal method for integration of multimodal single cell data [link](https://doi.org/10.1093/bioinformatics/btac220)

citation:

```bibtex
@article{cheng2022mojitoo,
  title={MOJITOO: a fast and universal method for integration of multimodal single-cell data},
  author={Cheng, Mingbo and Li, Zhijian and Costa, Ivan G},
  journal={Bioinformatics},
  volume={38},
  number={Supplement\_1},
  pages={i282--i289},
  year={2022},
  publisher={Oxford University Press}
}
```


### Code from MOJITOO Benchmarking

https://github.com/CostaLab/MOJITOO-reproducibility.git

### Data Set

This is deposited in zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6348128.svg)](https://doi.org/10.5281/zenodo.6348128)
