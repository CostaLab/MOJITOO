from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import StandardScaler
from scipy import stats
from statsmodels.stats.multitest import multipletests
from typing import Iterable, List, Union, Optional, Set, Tuple, TypeVar
from collections import Counter
from anndata import AnnData
from mudata import MuData

def mojitoo(adata: Union[AnnData, MuData]=None,
            reduction_list:list = [],
            dims_list:list = [],
            reduction_name:str = "X_mojitoo",
            is_reduction_center:str = False,
            is_reduction_scale:str = False,
            fdr_method:str = "fdr_bh",
            corr_pval:float = 0.05,
            overwrite:bool = False,
            iscopy:bool = False,
            **kargs
        ):
    """
    MOJITOO multimodal integration

    Parameters
    ----
    adata: AnnData or MuData
        anndata.Anndata or mudata.MuData object
    reduction_list: list
        reductions in adata.obsm
    dims_list : list
        dims for each dimension reduction to use, e.g. range(1,30), empty list indicate all of the dimensions
    reduction_name:str
       save_name, default: X_mojitoo 
    is_reduction_center: bool
        if center before cca, default: False
    is_reduction_scale: bool
        if scale before cca, default: False
    fdr_method: str
        fdr methods to use, candidates: bonferroni, sidak, holm-sidak, holm, simes-hochberg, hommel, fdr_bh, fdr_by, fdr_tsbh, fdr_tsbky
    corr_pval: float
        corrected pval threshold, default: 0.05
    overwrite: bool
        if overwrite the cca name, defalt: False
    iscopy: bool
        if copy adata, default: False
    """
    if adata is None:
        raise ValueError("Please provide an AnnData or MuData object.")

    if not overwrite and reduction_name in adata.obsm.keys():
        raise ValueError("reduction name exists, please enable parameter overwrite and re-try")

    if len(reduction_name) < 2:
        raise ValueError("At least 2 dimension names in reduction_list")
    if dims_list is None:
        dims_list = []
    #print(len(reduction_list), len(dims_list))
    if len(dims_list) !=0 and len(dims_list) !=len(reduction_list):
        raise ValueError("dims_list should be consistent with reduction_list")

    assert(fdr_method in ["bonferroni",
                          "sidak",
                          "holm-sidak",
                          "holm",
                          "simes-hochberg",
                          "hommel",
                          "fdr_bh",
                          "fdr_by",
                          "fdr_tsbh",
                          "fdr_tsbky"])

    assert(corr_pval >0 and corr_pval < 1)
 
    for redu in reduction_list:
        if redu not in adata.obsm.keys():
            raise ValueError("%s not in adata.obsm" % redu)


    adata = adata.copy() if iscopy else adata


    if len(dims_list) == 0:
        for name in reduction_list:
            dims_list.append(range(adata.obsm[name].shape[1]))

    a_redu = None
    for i in range(len(reduction_list) -1):
        if i == 0:
            a_name = reduction_list[i]
            a_redu = adata.obsm[a_name]
            a_dims = list(dims_list[i])
            a_redu = a_redu[:, a_dims]
        b_name = reduction_list[i+1] 
        b_redu = adata.obsm[b_name]
        b_dims = list(dims_list[i+1])
        b_redu = b_redu[:, b_dims]

        ## center or scale the data
        if is_reduction_center or is_reduction_scale:
            scaler = StandardScaler(with_mean=is_reduction_center,
                                    with_std=is_reduction_scale) 
            a_redu = scaler.fit_transform(a_redu) #scale data
            b_redu = scaler.fit_transform(b_redu)

        ## CCA transformation
        cca = CCA(n_components=min(len(a_dims), len(b_dims)), **kargs)
        cca.fit(a_redu, b_redu)
        a, b = cca.transform(a_redu, b_redu) 

        ## correlation test
        correlation_test = [stats.pearsonr(a[:, i], b[:, i])[1] for i in range(a.shape[1])]
        fdr_bool = multipletests(correlation_test,
                                alpha=corr_pval,
                                method=fdr_method)[0]

        print(f"{i+1} round cc", Counter(fdr_bool)[True])
        cca_add = a[:, fdr_bool] + b[:, fdr_bool]

        a_redu = cca_add
    adata.obsm[reduction_name] = cca_add

    return adata if iscopy else None

