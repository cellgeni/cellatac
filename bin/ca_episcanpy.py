#!/usr/bin/python3

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sci
import scipy as ss
import episcanpy as epi
import argparse

min_score_value = 0.515
nrof_features   = 120000
##ncpus           = 1 (if needed change this to be this to be equal to cpus submitted for job)

my_parser = argparse.ArgumentParser()

my_parser.add_argument("-m", "--min-score-value", type=float)
my_parser.add_argument("-f", "--nrof-features", type=int)
my_parser.add_argument("-p", "--ncpus", type=int)

args = my_parser.parse_args()

if args.min_score_value is not None:
  min_score_value = args.min_score_value

if args.nrof_features is not None:
  nrof_features = args.nrof_features

if args.ncpus is not None:
  ncpus = args.ncpus


adata = ad.read_mtx("mmtx/filtered_window_bc_matrix.mmtx.gz")
adata = adata.T
region_names = pd.read_table("mmtx/regions.names", header = None)
cell_names = pd.read_table("mmtx/cell.names", header = None)
adata.X = ss.sparse.csr_matrix(adata.X)
adata.obs.index = cell_names[0]
adata.var.index = region_names[0]
adata.var.index.name = 'region_names'
adata.obs.index.name = 'cell_names'
epi.pp.pca(adata)
epi.pp.neighbors(adata, method="umap")
epi.tl.umap(adata)
epi.pp.normalize_total(adata)
epi.pp.cal_var(adata, save='episcanpy_variability.pdf')
epi.pl.variability_features(adata,log=None, min_score=min_score_value, nb_features=nrof_features, save='episcanpy_varfeat.pdf')
epi.pl.variability_features(adata,log='log10', min_score=min_score_value, nb_features=nrof_features, save='episcanpy_log_varfeat.pdf')
adata = epi.pp.select_var_feature(adata, nb_features=nrof_features, show=False, copy=True)
epi.tl.leiden(adata)
# epi.tl.getNClusters(adata, n_cluster=N_clusters, method='leiden')   # note constant.
# umap has to be saved as .png as .pdf produces empty pdf files
epi.pl.umap(adata, color='leiden', wspace=0.4, save='.png')
adata.write("adata.h5ad", compression='gzip')
pd.DataFrame(adata.obs['leiden']).to_csv('Leiden.tsv', sep="\t", header=False)
