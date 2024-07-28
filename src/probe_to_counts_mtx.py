import pandas as pd
import scipy
from scipy.sparse import csr_matrix
from scipy import sparse
from scipy.io import mmwrite
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--probe_count')
args = parser.parse_args()

file_in = pd.read_csv(args.probe_count,delimiter="\t",header=0,names=["cellID","geneID"])
dat=pd.crosstab(file_in["geneID"],file_in["cellID"])
csr_matrix = csr_matrix(dat.astype(pd.SparseDtype("float64",0)).sparse.to_coo())
mmwrite("probe_count_matrix.mtx",csr_matrix,field='integer')
dat.index.to_frame().to_csv("features.tsv",mode="w",sep="\t")
dat.T.index.to_frame().to_csv("barcodes.tsv",mode="w",sep="\t")