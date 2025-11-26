import scanpy as sc
import pandas as pd
from scipy.io import mmread

X = mmread(
    "GSE114412_counts.mtx"
).T.tocsr()  # transpose because AnnData expects cells × genes
obs_cluster = pd.read_csv("meta_cluster.csv", index_col=0)
var = pd.read_csv("var.csv")

adata = sc.AnnData(X=X, obs=obs_cluster, var=var)
adata.write_h5ad("GSE114412_cluster.h5ad")

obs_subcluster = pd.read_csv("meta_subcluster.csv", index_col=0)

adata = sc.AnnData(X=X, obs=obs_subcluster, var=var)
adata.write_h5ad("GSE114412_subcluster.h5ad")
