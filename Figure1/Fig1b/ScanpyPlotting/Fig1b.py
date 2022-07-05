##loading libraries
import numpy as np
import pandas as pd
import scanpy as sc

###loading data and metadata
exp = pd.read_csv("markers_for_scanpy.csv",index_col=0)
CellTypes = pd.read_csv("groups.csv",index_col=0)

##plotting upregulated genes in Cancer and NK cells
adata = sc.AnnData(X = exp, obs = CellTypes)
gn_up = ["HMGA1","TACSTD2","ANKRD11"]
sc.pl.matrixplot(adata, gn_up, groupby="CellTypes",var_group_rotation=0,save=True,cmap='Blues')

gn_up = ["KLRD1","LAIR1","TNFRSF9","CCR6"]
sc.pl.matrixplot(adata, gn_up, groupby="CellTypes",var_group_rotation=0,save=True,cmap='Blues')
