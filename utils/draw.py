import numpy as np
import scipy as sp
import scipy.ndimage
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def generateColor():
    ret = {'CA1': ('#f67088', ['#cc5d70', '#994554', '#7f3a46', '#ff748d', '#da7184', '#ae5069', '#f56077', '#f18099', '#e26265', '#ff8cae', '#eb7e70', '#fc5ea5']), 'CA3': ('#ad9c31', ['#cbb83a', '#998a2b', '#7f7324', '#ffe648', '#dad44f', '#ecc63f', '#b39d24', '#ffff58', '#eee22e', '#ece45b', '#e1c222', '#b8c449']), 'SUB': ('#33b07a', ['#3bcb8d', '#2c996a', '#257f58', '#4affb1', '#33ae7f',
                                                                                                                                                                                                                                                                                                                                                   '#3ee99e', '#37c577', '#5ef3a0', '#59fbcc', '#4ed49c', '#40dbbc', '#47fc8f']), 'DG': ('#38a8c5', ['#39aecc', '#2b8399', '#246d7f', '#48daff', '#42bde6', '#3a95b7', '#38d2dc', '#289ca5', '#4bffff', '#51a1c6', '#2ff1f4', '#4deade']), 'paraS': ('#cc79f4', ['#aa65cc', '#7f4c99', '#6a3f7f', '#d57eff', '#915ab3', '#c36ef2', '#ab75eb', '#c287dd', '#ef92ff', '#e671f4', '#a863a5', '#b678bb'])}
    ret["proS"] = ret["SUB"]
    ret["preS"] = ret["paraS"]
    ret["CA2"] = ret["CA3"]
    ret["CA4"] = ret["CA3"]
    ret["FC"] = ret["SUB"]
    return ret

def draw_sizebar(ax, length, r, dark=False):
    """
    Draw a horizontal bar with length of 0.1 in data coordinates,
    with a fixed label underneath.
    """
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    if dark:
        asb = AnchoredSizeBar(ax.transData,
                            length,
                            r,
                            loc='lower right', frameon=False, color="white")
    else:
        asb = AnchoredSizeBar(ax.transData,
                            length,
                            r,
                            loc='lower right', frameon=False)
    # pad=0.1, borderpad=0.5, sep=5,
    # frameon=False)
    ax.add_artist(asb)

def centralize(adata, palette, col):
    df = adata.obs.copy()
    df["color"] = df[col].map(palette)

    df["X"] = adata.obsm["spatial"][:, 1]-(np.max(adata.obsm["spatial"][:, 1]) + np.min(adata.obsm["spatial"][:, 1])) / 2
    df["Y"] = adata.obsm["spatial"][:, 0]-(np.max(adata.obsm["spatial"][:, 0]) + np.min(adata.obsm["spatial"][:, 0])) / 2
    # print(df.columns)
    return df, np.max(df["X"]), np.max(df["Y"])

def drawSlice(adata, palette, col="mapped", margin=20, dark=False, len=40, r="1mm"):
    df, x, y = centralize(adata, palette, col)
    # plt.axis('scaled')
    fig = plt.figure(figsize=(5,5))
    axis = fig.gca()
    axis.spines['top'].set_visible(False)
    axis.spines['bottom'].set_visible(False)
    axis.spines['left'].set_visible(False)
    axis.spines['right'].set_visible(False)
    lim = max(x, y) + margin
    axis.set_xlim(-lim, lim)
    axis.set_ylim(-lim, lim)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.scatter(df["Y"], df["X"], c=df["color"], s=2)
    draw_sizebar(axis, len, r, dark)

def drawSliceCmap(adata, cmap, col="mapped", dark=False, len=40, r="1mm"):
    # plt.axis('equal')
    x = np.max(adata.obs["x"])
    y = np.max(adata.obs["y"])
    fig = plt.figure()
    axis = fig.gca()

    axis.set_aspect('equal', adjustable='box')
    axis.spines['top'].set_visible(False)
    axis.spines['bottom'].set_visible(False)
    axis.spines['left'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.set_xlim(0, x+2)
    axis.set_ylim(0, y+2)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.scatter(adata.obs["x"], adata.obs["y"], c=adata.X[:, adata.var_names.get_loc(col)], s=2, cmap=cmap)
    draw_sizebar(axis, len, r, dark)

    return axis

def toLongTable(df: pd.DataFrame, row="row", col="col", val="val"):
    rows = []
    cols = []
    vals = []
    df = df.T
    for i, r in df.iterrows():
        for j in r.index:
            vals.append(r[j])
            cols.append(i)
            rows.append(j)
    longdf = pd.DataFrame({row: rows, col: cols, val: vals})
    return longdf


def drawFacet(df:pd.DataFrame, palette):
    # longdf = toLongTable(df)
    longdf = df
    sns.set(rc={'axes.facecolor':(0, 0, 0, 0)})
    # labels = list(df.index)
    # palette = []
    # for i in labels:
    #     if i.startswith("Glu"):
    #         palette.append((0.2528186785662627, 0.6322661966470429, 0.9586861264495917))
    #     else:
    #         palette.append((0.9692894417585417, 0.4522225702495641, 0.4261820543616833))

    g = sns.FacetGrid(longdf, row="row", hue="row", aspect=15,height=1, palette=palette)
    g.map(sns.lineplot, "EBZ", "val", linewidth=4)
    g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)
    labels = list(df.index)

    def label(x, label="A", color=None, ):
        ax = plt.gca()
        # global labels
        # label = labels[0]
        # labels = labels[1:]
        ax.text(-0.13, 0, label, fontweight="bold", color="black",
                ha="left", va="center", fontsize=15, transform=ax.transAxes)
    g.map(label, "row")

    def fill(x, val, EBZ, color, label):
        ax = plt.gca()
        # print(val)
        # print(EBZ)
        # print(color)
        # print(label)
        ax.fill_between(EBZ, val.values, alpha=0.5, color=color)

    g.map(fill, "row", "val", "EBZ")

    g.figure.subplots_adjust(hspace=-0.5)
    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.set(xticks=[], xlabel="")
    g.set_xlabels("")
    g.despine(bottom=True, left=True)

def prepareLongTable(df, dfEBZ):
    longdf = toLongTable(df)
    longdf["EBZ"] = -longdf["col"].map(dfEBZ["EBZ"])
    longdf = longdf.dropna()
    return longdf

def scatterCelltype(df: pd.DataFrame, celltype, seq, x="x", y="y"):
    num = len(seq)
    rows = (num+5) // 6
    plt.figure(figsize=(24, rows*4))
    for ii in range(len(seq)):
        i = seq[ii]
        dfsub = df[df["Slice"] == i]
        dfsub = dfsub[dfsub["Cell_Type"] == celltype]

        # print(df.head())
        plt.subplot(rows, 6, ii+1)
        plt.title(i)
        plt.scatter(dfsub[x], dfsub[y], s=1)

    plt.show()
    plt.close()


def toHex(c):
    ret = "#"
    for i in c:
        if i > 1.0:
            i = 1.0
        t = hex(int(i*255))
        if len(t) == 3:
            ret += "0"+t[-1]
        else:
            ret += t[-2:]
    return ret[:7]


def toDec(c):
    ret = np.zeros((3))
    ret[0] = int(c[1:3], 16) / 256
    ret[1] = int(c[3:5], 16) / 256
    ret[2] = int(c[5:], 16) / 256
    return ret 

def getPaletteColor(c):
    c = toDec(c)
    return c/np.max(c)


def ajust(v, a, b):
    v = v*a+b
    v[v > 1.0] = 1.0
    v[v < 0.0] = 0.0
    return v


def prepare(adata, l, quant=0.9):
    maxl = {}
    for i in l:
        ll = adata.var_names.get_loc(i)
        ll_max = np.quantile(adata.X[:, ll], quant)
        maxl[i] = ll_max
        adata.obs[i] = adata.X[:, ll]
    return maxl


def getc(genes, max_val, thres=2, norm=None):
    def xx(x):
        a = np.zeros((3))
        for i in genes:
            weight = min(max(0, (x[i] - max_val[i]/thres) /
                         (max_val[i] - max_val[i]/thres)), 1)
            if norm is not None:
                if x[norm] > 1:
                    weight /= x[norm]
            a += weight*genes[i]
        return toHex(a)
    return xx


def getcc(genes, max_val, thres=2):
    def xxx(x):
        a = np.zeros((3))
        for i in genes:
            a += min(max(0, (x[i] - max_val[i]/thres) /
                     (max_val[i] - max_val[i]/thres)), 1)*genes[i]
        return np.max(a)
    return xxx


def getcolor(adata:sc.AnnData, genes:dict, max=None, thres=0.1, thres2=2, quant=0.9):
    """calculate the color of each spot for an Anndata object given a list of genes to draw

    Args:
        adata (sc.AnnData): data anndata
        genes: dictionary object recording color for each gene
        max (optional): dictionary object for max value for each gene
        thres (float, optional): omit color less than . Defaults to 0.1.
        thres2 (int, optional): map (1/thres2, 1) -> (0, 1) for color code. Gene
         expression level of less then 1/thres2 will be mapped to  Defaults to 2.
        quant (float, optional): The quantile to calculate max value if max is not given. Defaults to 0.9.

    Returns:
        _type_: _description_
    """
    maxl = prepare(adata, genes.keys(), quant=quant)
    if max is not None:
        for i in maxl:
            maxl[i] = max
    adata.obs["cc"] = adata.obs.apply(getcc(genes, maxl, thres=thres2), axis=1)
    adata.obs["c"] = adata.obs.apply(
        getc(genes, maxl, thres=thres2, norm=None), axis=1)
    p = {}
    for i in adata.obs["c"]:
        p[i] = i
    sub = adata[adata.obs["cc"] > thres]
    return sub, p

import matplotlib as mpl
from matplotlib.colors import ListedColormap
gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))
def plot_color_gradients(category, cmap_list):
    # Create figure and adjust figure height to number of colormaps
    nrows = len(cmap_list)
    figh = 0.35 + 0.15 + (nrows + (nrows - 1) * 0.1) * 0.22
    fig, axs = plt.subplots(nrows=nrows + 1, figsize=(6.4, figh))
    fig.subplots_adjust(top=1 - 0.35 / figh, bottom=0.15 / figh,
                        left=0.2, right=0.99)
    axs[0].set_title(f'{category} colormaps', fontsize=14)

    for ax, name in zip(axs, cmap_list):
        if isinstance(name, str):
            ax.imshow(gradient, aspect='auto', cmap=mpl.colormaps[name])
        else:
            ax.imshow(gradient, aspect='auto', cmap=name)
        if isinstance(name, str):
            ax.text(-0.01, 0.5, name, va='center', ha='right', fontsize=10,
                    transform=ax.transAxes)
        else:
            ax.text(-0.01, 0.5, "test", va='center', ha='right', fontsize=10,
                    transform=ax.transAxes)

    # Turn off *all* ticks & spines, not just the ones with colormaps.
    for ax in axs:
        ax.set_axis_off()


def getColorMap():
    first = 30 
    next = 30 
    last = 200
    jet = mpl.colormaps["gnuplot2"]
    jets = jet(np.linspace(0, 1, 256))

    hot = mpl.colormaps["hot"]
    hots = hot(np.linspace(0, 1, 256))

    need = np.zeros([first+next+last, 4])
    need[:first] = jets[40:40+first]

    for i in range(next):
        need[first+i] = ((next-i-1)*jets[40+first] + (i+1)*hots[-last]) / next 
    need[-last:] = hots[-last:]

    # print(jets)
    newcmp = ListedColormap(need[:200])

    return newcmp

from scipy.ndimage import gaussian_filter1d
def drawCellProportion(x:np.array, y:np.array, axis, label=None, c="blue"):
    axis.scatter(x, y, c=c)
    # sq = np.array(y)
    sq = y**2
    sq = gaussian_filter1d(sq, sigma=2) 
    t = gaussian_filter1d(y, sigma=2) 
    sq = sq - t**2
    sq = np.sqrt(sq)
    axis.plot(x, t, linestyle="--", label=label, color=c)
    axis.fill_between(x, (t-sq), (t+sq), color=c, alpha=.1)