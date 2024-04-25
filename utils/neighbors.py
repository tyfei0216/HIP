from scipy.stats import gmean
import scanpy as sc 
import numpy as np
from sklearn.neighbors import NearestNeighbors 
from tqdm import tqdm 
import pandas as pd

def markGLu(x):
    if x.startswith("Glu"):
        return 1 
    else:
        return 0

def markGABA(x):
    if x.startswith("GABA"):
        return 1 
    else:
        return 0

def filtercells(excitory, inhibitory):

    def xx(x):
        if x.startswith("Glu"):
            if x in excitory:
                return x 
            else:
                return "NAN"
            
        if x.startswith("GABA"):
            if x in inhibitory:
                return x
            else:
                return "NAN"
            # return x[:x.rfind(" ")]
        return x
    return xx

def getType(x):
    if x.startswith("Glu"):
        return "Glu"
    elif x.startswith("GABA"):
        return "GABA"
    return "NONE"

def _getNeighbors(adata, celltype, n_neighbors):
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='ball_tree').fit(adata.obsm["spatial"])
    sub1 = adata[adata.obs["cellType"] == celltype] 
    _, indices = nbrs.kneighbors(sub1.obsm["spatial"])
    indices = indices[:, 1:]
    return sub1.obs.index, np.array(adata.obs.index)[indices]

from scipy.stats import binomtest

def test(adata:sc.AnnData, f:str, t:str, gene1:str, gene2:str, pairs, method="binomtest"):
    """
    test whether a ligand-receptor pair id significant between two cell types, 
    Args:
        adata (sc.AnnData): the anndata object for all cells
        f (str): the ligand cell type name
        t (str): the receptor cell type name
        gene1 (str): ligand gene name
        gene2 (str): receptor gene name
        pairs (_type_): cell neighbors pairs for testing
        method (str, optional): method for testing. Usually bionomtest have better result than 
        permutation tests. Defaults to "binomtest".

    Returns:
        the significance p value
    """
    sub1 = adata[adata.obs["cellType"] == f]
    p1 = sub1.X[:, adata.var_names.get_loc(gene1)].toarray().squeeze()
    p1 = (p1 > 0.001).sum() / len(p1) 
    sub2 = adata[adata.obs["cellType"] == t] 
    p2 = sub2.X[:, adata.var_names.get_loc(gene2)].toarray().squeeze()
    p2 = (p2 > 0.001).sum() / len(p2)
    fromgene = adata.X[:, adata.var_names.get_loc(gene1)].toarray().squeeze()
    fromgene = pd.Series(fromgene, index=adata.obs.index)
    togene = adata.X[:, adata.var_names.get_loc(gene2)].toarray().squeeze()
    togene = pd.Series(togene, index=adata.obs.index)
    cnt = 0 
    for i, j in pairs:
        if fromgene[i] > 0 and togene[j] > 0:
            cnt += 1 
    # print(cnt, len(pairs), p1*p2)
    if method == "binomtest":
        res = binomtest(cnt, len(pairs), p1*p2, "greater")
        return res.pvalue
    elif method == "permutation":
        TRIALS = 100000 
        success = 0 
        t1 = np.array(sub1.obs.index)
        t2 = np.array(sub2.obs.index) 
        for _ in tqdm(range(TRIALS)):
            t1_shuffle = t1.copy()
            t2_shuffle = t2.copy()
            np.random.shuffle(t1_shuffle)
            np.random.shuffle(t2_shuffle)
            mapt1 = {} 
            mapt2 = {}
            for i, j in zip(t1, t1_shuffle):
                mapt1[i] = j 
            for i, j in zip(t2, t2_shuffle):
                mapt2[i] = j 
            tcnt = 0 
            for i, j in pairs:
                if fromgene[mapt1[i]] > 0 and togene[mapt2[j]] > 0:
                    tcnt += 1
            if tcnt > cnt:
                success += 1 
        return success / TRIALS 

def testPair(adatas, cell1, cell2, genepairs, cell1_thres=30, cell2_thres=10, pair_thres = 20):
    res = {} 
    for i in genepairs:
        res[i] = []
    for adata in tqdm(adatas):
        # print(adata)
        adata = adatas[adata]
        cells = adata.obs["cellType"].value_counts() 
        if cells[cell1] > cell1_thres and cells[cell2] > cell2_thres:
            f, t = _getNeighbors(adata, cell1, 32)
            t2 = adata.obs.index[adata.obs["cellType"] == cell2]
            pairs = [] 
            for i in range(t.shape[0]):
                for j in range(t.shape[1]):
                    if t[i][j] in t2:
                        pairs.append((f[i], t[i][j]))
            # print(len(pairs))
            if len(pairs) < pair_thres:
                continue
            for i in genepairs:
                pval = test(adata, cell1, cell2, i[0], i[1], pairs, method="binomtest")
                res[i].append(pval)
    return res 

def toOne(t):
    ret = {} 
    for i in t:
        if len(t[i]) > 0:
            ret[i] = gmean(t[i])
        else:
            ret[i] = 1
    return ret

def getNeighbors(adatacells, n_neighbors = [32, 16, 32, 32]):
    results = {}
    for i in adatacells:
        adata = adatacells[i]
        celltypes = [] 
        for celltype in np.unique(adatacells[i].obs["cellType"]):
            if celltype.startswith("Glu"):
                celltypes.append(celltype)
        # print(celltypes)
        results[i] = {}
        for maintype, num in zip(["Glu", "GABA", "NONE", "Total"], n_neighbors):
            if maintype != "Total":
                subadata = adata[adata.obs["type"] == maintype]
            else:
                subadata = adata
            nbrs = NearestNeighbors(n_neighbors=num, algorithm='ball_tree').fit(subadata.obsm["spatial"])
            t = {}
            for celltype in celltypes:
                sub1 = adata[adata.obs["cellType"] == celltype]
                # if len(sub1) < 50:
                #     print(i, celltype, "<50", len(sub1))
                    # continue
                _, indices = nbrs.kneighbors(sub1.obsm["spatial"])
                indices = indices[:, 1:]
                tt = np.array(subadata.obs["cellType"])
                unique, counts = np.unique(tt[indices.flatten()], return_counts=True)
                # print(len(sub1), (num-1)*len(sub1), counts.sum())
                t[celltype] = pd.Series(counts, index=unique)
            df = pd.DataFrame(t)
            # sumdf = df.sum(axis=1)
            results[i][maintype] = df
    return results

def getNeighborEnrichmentScore(neighbors, cell_list):
    enrich_score = {}
    for maintype in ["Glu", "GABA", "NONE", "Total"]:
        # enrich_score[maintype] = {}
        base = None
        for i in neighbors:
            neighbors[i][maintype] = neighbors[i][maintype].fillna(0)
            
            for j in cell_list:
                # if j == 'Glu SUB int. 8':
                #     continue
                if j not in neighbors[i][maintype].columns:
                    neighbors[i][maintype][j] = 0
            
            neighbors[i][maintype] += 1
            neighbors[i][maintype] = neighbors[i][maintype][cell_list]
            if base is None:
                base = neighbors[i][maintype].copy()
            else: 
                base = base.add(neighbors[i][maintype], fill_value=0)
            # if (results[i][maintype].sum(axis=1) < 500).any():
            #     print(i, maintype, results[i][maintype].sum(axis=1))
            # results[i][maintype] = results[i][maintype].loc[results[i][maintype].sum(axis=1) > 500]
            neighbors[i][maintype] = neighbors[i][maintype].div(neighbors[i][maintype].sum(axis=0), axis=1)
            # enrich_score[maintype][i] = results[i][maintype].mean(axis=1)
        # enrich_score[maintype] = pd.DataFrame(enrich_score[maintype])
        base = base.div(base.sum(axis=0), axis=1)
        enrich_score[maintype] = base.mean(axis=1)
    return enrich_score

def normalizeScore(enrich_score, neighbors):
    for maintype in ["Glu", "GABA", "NONE", "Total"]:
        for i in neighbors:
            subdf = enrich_score[maintype].loc[neighbors[i][maintype].index]
            neighbors[i][maintype] = neighbors[i][maintype].div(subdf, axis=0)
            neighbors[i][maintype] += 1e-6
            neighbors[i][maintype] = np.log(neighbors[i][maintype])

    return neighbors

def fdr(p_vals):

    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 0.9999999

    return fdr