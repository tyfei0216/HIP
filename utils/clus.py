import scanpy as sc 
import pandas as pd 
import os 
import numpy as np
import scipy as sp
import anndata as ad
from scipy.sparse import spmatrix
Logs = None 
LogPath = None

fields = ["output file", "smooth", "high variable", "neighbor", "npc", "resolution"]

def buildEmptyRecord():
    global Logs
    global fields
    a = pd.DataFrame(None, columns=fields)
    a = a.set_index("output file")
    Logs = a 
    return a 

def readTable():
    global LogPath
    global fields
    global Logs
    if LogPath is None:
        raise ValueError("log path not initialized")
    
    if os.path.exists(os.path.join(LogPath, "log.csv")):
        df = pd.read_csv(os.path.join(LogPath, "log.csv"), index_col=None)
        if len(np.intersect1d(fields, df.columns)) != len(fields):
            raise ValueError("log file doesn't meet requirements")
        df = df.set_index("output file")
    else:
        df = buildEmptyRecord()
    Logs = df 

def saveLogs():
    global LogPath 
    global Logs 
    if LogPath is None:
        raise ValueError("log path not initialized")

    if Logs is None:
        raise ValueError("Log table is empty")

    Logs.to_csv(os.path.join(LogPath, "log.csv"))
    
def initLogs(path):
    global LogPath
    
    if not os.path.exists(path):
        os.makedirs(path)

    LogPath = path 

    readTable()

def gau(mat, s = 0.8, truc=4):
    V = mat.copy()
    V[np.isnan(V)] = 0 
    VV = sp.ndimage.gaussian_filter(V, sigma=s, truncate=truc)

    W = np.ones_like(mat)
    W[np.isnan(mat)] = 0
    WW = sp.ndimage.gaussian_filter(W, sigma=s, truncate=truc)

    return VV/(WW+1e-4)

def newclass(g, t, s = 0.8, trunc=4):
    # t = adata.obsm["spatial"]
    t[:, 0] -= np.min(t[:, 0])
    t[:, 1] -= np.min(t[:, 1])
    # adata.obsm["spatial"] = t
    x = np.max(t[:, 0])+1 
    y = np.max(t[:, 1])+1
    mask = np.empty((x, y))
    mask.fill(np.nan)
    for i in t:
        mask[i[0], i[1]] = 0
    ss, ll = g.shape
    mat = np.empty((ll, x, y))
    for i in range(ss):
        xx = t[i, 0]
        yy = t[i, 1]
        mat[:, xx, yy] = g[i]
    for i in range(ll):
        mat[i] = gau(mat[i], s=s, truc=trunc)
    new = []
    for i in range(ss):
        xx = t[i, 0]
        yy = t[i, 1]
        # g.iloc[i] = mat[:, xx, yy] 
        new.append(mat[:, xx, yy])
    new = np.array(new)

    df = np.array(new)
    return df

def smoothadata(adata, s = 0.8, trunc=4):
    t = adata.obsm["spatial"].copy()
    t[:, 0] -= np.min(t[:, 0])
    t[:, 1] -= np.min(t[:, 1])
    # adata.obsm["spatial"] = t
    x = np.max(t[:, 0])+1 
    y = np.max(t[:, 1])+1
    # mask = np.empty((x, y))
    # mask.fill(np.nan)
    # for i in t:
    #     mask[i[0], i[1]] = 0
    _, ll = adata.shape
    mat = np.empty((ll, x, y))
    mat.fill(np.nan)
    for i in range(len(adata)):
        xx = t[i, 0]
        yy = t[i, 1]
        mat[:, xx, yy] = adata.X[i]
    for i in range(ll):
        mat[i] = gau(mat[i], s=s, truc=trunc)
    new = []
    for i in range(len(adata)):
        xx = t[i, 0]
        yy = t[i, 1]
        # g.iloc[i] = mat[:, xx, yy] 
        new.append(mat[:, xx, yy])
    new = np.array(new)
    adata.raw = adata 
    # adata.layers["beforeSmooth"] = adata
    adata.X = new
    return adata

def getSmooth(T25, s=0.8):

    print("smooth adata with s=", s)
    
    T25.obsm["spatial"] = T25.obsm["spatial"].astype(np.int32)
    
    T25_1 = T25[T25.obs["gene_area"]!="DG", ]
    df = newclass(T25_1.X.copy(), T25_1.obsm["spatial"].copy(), s=s)
    newdf = pd.DataFrame(df, index=T25_1.obs.index, columns=T25_1.var.index)
    newT25_1 = sc.AnnData(newdf, dtype=float)
    newT25_1.obs["gene_area"] = T25_1.obs["gene_area"]
    newT25_1.obsm["spatial"] = T25_1.obsm["spatial"]
    newT25_1 = newT25_1[newT25_1.obs["gene_area"]!="CA4",]
    # sc.pl.spatial(newT25_1, color="gene_area", spot_size=1)
    T25_2 = T25[(T25.obs["gene_area"]=="DG")|(T25.obs["gene_area"]=="CA4")|(T25.obs["gene_area"]=="CA3"), ]
    df = newclass(T25_2.X.copy(), T25_2.obsm["spatial"].copy(),s=s)
    newdf = pd.DataFrame(df, index=T25_2.obs.index, columns=T25_2.var.index)
    newT25_2 = sc.AnnData(newdf, dtype=float)
    newT25_2.obs["gene_area"] = T25_2.obs["gene_area"]
    newT25_2.obsm["spatial"] = T25_2.obsm["spatial"]
    newT25_2 = newT25_2[newT25_2.obs["gene_area"]!="CA3",]
    # sc.pl.spatial(newT25_2, color="gene_area", spot_size=1)
    adatas = [newT25_1, newT25_2]
    newA = ad.concat(adatas)
    newA.raw = T25
    print("finish smoothing")
    return newA

def process(adata:sc.AnnData, s=0.4, high_variables=-1, neighbor=20, npc=30, resolution=0.8, smooth="smoothadata", add_key="leiden") -> sc.AnnData:
    """ 
    The clustering process: do smoothing for all genes then carry out the standard single cell clustering pipeline
    

    Args:
        adata (sc.AnnData): data Anndata object
        s (float, optional): The variance for smoothing. Defaults to 0.4.
        high_variables (int, optional): The number of high variable genes genes to calculate pca. -1 means use all genes. Defaults to -1.
        neighbor (int, optional): The number of neighbors for leiden clustering. Defaults to 20.
        npc (int, optional): The number of pc for leiden clustering. Defaults to 30.
        resolution (float, optional): The resolution for leiden clustering. Defaults to 0.8.
        smooth (str, optional): smoothing method, just use the smoothadata function. getsmooth trys to do smoothing by separating DG 
        and CA. Defaults to "smoothadata".
        add_key (str, optional): The column name for saving clustering result in the obs. Defaults to "leiden".


    Returns:
        sc.AnnData: The clustering result
    """
    if "process" in adata.uns:
        # print(adata.uns["process"]["smooth"])
        if np.abs(adata.uns["process"]["smooth"] - s) > 0.001:
            adata.uns.pop("process")
            X = adata.X 
            if isinstance(X, spmatrix):
                X = X.toarray()
            adata = sc.AnnData(X, obs = adata.obs, var = adata.var, obsm = adata.obsm, uns = adata.uns)
            if adata.raw is None:
                raise ValueError("cannot find raw data from anndata")

            adata = adata.raw.to_adata()
            if smooth == "getsmooth":
                adata = getSmooth(adata, s)
            elif smooth == "smoothadata":
                adata = smoothadata(adata, s)
            else:
                raise NotImplementedError

            sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
            sc.pp.normalize_total(adata, target_sum=1e6)
            sc.pp.log1p(adata)
    else:
        X = adata.X 
        if isinstance(X, spmatrix):
            X = X.toarray()
        adata = sc.AnnData(X, obs = adata.obs, var = adata.var, obsm = adata.obsm, uns = adata.uns)
        if smooth == "getSmooth":
            adata = getSmooth(adata, s)
        elif smooth == "smoothadata":
            adata = smoothadata(adata, s)
        else:
            raise NotImplementedError
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)

    if "process" in adata.uns:
        if adata.uns["process"]["high variable"] != high_variables:            
            adata.uns.pop("process")
            print("find high variable n=", high_variables)
            
            if high_variables > 0:
                sc.pp.highly_variable_genes(adata, n_top_genes=high_variables)
            else:
                adata.var["highly_variable"] = True
                # sc.pp.regress_out(adata, ['total_counts'])
            
            # sc.pp.scale(adata, max_value=10)
    else:
        print("find high variable n=", high_variables)
        if high_variables > 0:
            sc.pp.highly_variable_genes(adata, n_top_genes=high_variables)
        else:
            adata.var["highly_variable"] = True
            # adata = adata[:, adata.var.highly_variable]
            # sc.pp.regress_out(adata, ['total_counts'])
            
        # sc.pp.scale(adata, max_value=10)

    if "process" not in adata.uns:
        print("finding neighbors n=", neighbor)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata, n_neighbors=neighbor, n_pcs=npc)
    else: 
        
        if adata.uns["process"]["neighbor"] != neighbor or adata.uns["process"]["npc"] != npc:
            print("finding neighbors n=", neighbor)
            sc.pp.neighbors(adata, n_neighbors=neighbor, n_pcs=npc)
            # sc.tl.umap(adata)

    print("clustering using leiden resolution=", resolution) 
    sc.tl.leiden(adata, resolution=resolution, key_added=add_key)

    adata.uns["process"] = {"smooth":s, "high variable":high_variables, "neighbor":neighbor, "npc":npc, "resolution":resolution}

    return adata
    
def saveadata(adata, name):
    global LogPath
    global Logs
    Logs.loc[name] = adata.uns["process"]
    saveLogs() 
    adata.write_h5ad(os.path.join(LogPath, name+".h5ad"))

    
