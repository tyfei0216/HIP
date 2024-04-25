import rpy2.robjects as ro 
import anndata2ri
anndata2ri.activate()
# %load_ext rpy2.ipython
ro.r('library(Seurat)')
ro.r('library(SingleCellExperiment)')

def doSCT(adata):
    # adata = sc.read_h5ad(path)
    try:
        adata.obs["x"] = adata.obsm["spatial"][:, 0]
        adata.obs["y"] = adata.obsm["spatial"][:, 1]
    except:
        print("no spatial field")
        pass
    ro.globalenv['res'] = adata
    adata = ro.r('''
    seurat_obj <- as.Seurat(res, counts=\"X\", data = NULL)
    seurat_obj <- SCTransform(object=seurat_obj, assay = "originalexp", return.only.var.genes = FALSE)
    as.SingleCellExperiment(seurat_obj)
    ''')
    adata.obsm["spatial"] = adata.obs[["x","y"]].values
    adata.raw = adata
    return adata 

from rpy2.robjects import pandas2ri
pandas2ri.activate()

def saveSeurat(adata, path):
    adata.obs["x"] = adata.obsm["spatial"][:, 0]
    adata.obs["y"] = adata.obsm["spatial"][:, 1]
    adata.uns = {}
    ro.globalenv['res'] = adata
    ro.r("seurat_obj <- as.Seurat(res, counts=\"X\", data = NULL)\nsaveRDS(seurat_obj, \"%s\")"%path)

def ReadDataframeFromRDS(path):
    df = ro.r("readRDS(\"%s\")"%path)
    return df 

def ReadAdataFromRDS(path):
    adata = ro.r("d1 <- readRDS(\"%s\")\nas.SingleCellExperiment(d1)"%path)
    return adata
