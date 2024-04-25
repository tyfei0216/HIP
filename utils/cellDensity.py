import pandas as pd 
import numpy as np 
import anndata


def df2adata(df, add=None):
    df["xbin"] = df["x"].astype(int)
    df["ybin"] = df["y"].astype(int)
    df["loci"] = df["xbin"]*100000+df["ybin"]
    g = df.groupby(["loci", "Cell_Type"])["cell"].count()
    # print(g.head())
    g = g.to_frame().reset_index()
    if add is not None:
        # loci = np.array(add.obsm["spatial"][:, 0]*100000+add.obsm["spatial"][:, 1])
        # loci = pd.DataFrame({"loci":loci})
        # loci["cell"] = 1
        # loci["Cell_Type"] = "na"
        # loci = loci.groupby(["loci", "Cell_Type"])["cell"].count()
        loci = np.array(add.obs["x"]*100000+add.obs["y"])
        add.obs["loci"] = loci
        loci = pd.DataFrame({"loci":loci})
        loci["cell"] = 1
        loci["Cell_Type"] = "na"
        # loci = loci.groupby(["loci", "Cell_Type"])["cell"].count()
        # print(loci.head())
        g = pd.concat([g, loci])
        # for i in loci:
        #     g[(i, "temp")] = 1
    
    g = g.pivot(index=['loci'],
                        columns=['Cell_Type'],
                        values=['cell']).fillna(0)
    g = g.droplevel(0, axis=1)
    g["sum"] = g.apply(sum, axis=1)
    # print(len(df), len(add))
    if add is not None:
        g /= (len(df)/len(add))
    adata = anndata.AnnData(g)
    adata.obs["loci"] = adata.obs.index.astype(int)
    adata.obs["x"] = adata.obs["loci"] // 100000
    adata.obs["y"] = adata.obs["loci"] % 100000
    adata.obsm["spatial"] = adata.obs[["x", "y"]].values
    adata.var_names
    return adata 