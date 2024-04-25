import scanpy as sc 
import numpy as np 
import pandas as pd 
import anndata as ad 

def mergeSingleAdata(adata, col="mapped"):
    df = pd.DataFrame(adata.X, columns=adata.var_names) 
    df["idx"] = list(adata.obs[col])
    df = df.groupby("idx").sum() 
    newadata = sc.AnnData(df, var=adata.var, obs=pd.DataFrame({"mapped":df.index}, index = df.index)) 
    return newadata

def merge(adatadict, col="mapped"):
    """add up all counts for a subregion. 

    Args:
        adatadict (_type_): raw data, make sure data is not normalized. 
        col (str, optional): subregion label column. Defaults to "mapped".
    """

    mergedadata = {}

    for i in adatadict:
        adata = adatadict[i] 
        newadata = mergeSingleAdata(adata, col=col)
        newadata.obs["slice"] = i 
        mergedadata[i] = newadata

    merged = ad.concat(mergedadata, label="slice")
    return merged



def Jaccard(list1, list2, k1 = 100, k2 = 200):
    u1 = np.union1d(list1[:k1], list2[:k1]) 
    u2 = np.intersect1d(list1[:k2], list2[:k2])

    t1 = np.intersect1d(u1, u2)

    return len(t1) / len(u1)

