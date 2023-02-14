import json
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad

# load config
config = json.load(
    open("MIRSORT_ANNOTATION_DF.json")
)

expression_path = config["assemble_anndata.py"]["expression_file"]
obs_path = config["assemble_anndata.py"]["obs_file"]
version = config["assemble_anndata.py"]["version"]
ad_out_path = config["assemble_anndata.py"]["ad_path"]
fa_out_path = config["assemble_anndata.py"]["fasta_path"]

def RPMNormalize(ad,use_log1p):
    ad_log = ad.copy()
    ad_log.obs["total_count"] = ad_log.obs["total_count_after_preprocessing"] 
    np.divide(ad_log.X, ad_log.obs["total_count"][:, None] / 1e6, out=ad_log.X)
    if use_log1p:
        sc.pp.log1p(ad_log, base=2)
    ad_log.raw = ad
    return ad_log


# load data and assemble as anndata
obs_df = pd.read_csv(obs_path, index_col=0, sep="\t")
new_ad = ad.read_csv(expression_path, delimiter="\t").T


# RPMNormalize and log1p
new_ad = RPMNormalize(new_ad,use_log1p=True)


# add version number to uns
new_ad.uns["version"] = version


# save anndata
new_ad.write_h5ad(ad_out_path)


# save seqs as fasta
varnames = new_ad.var_names
ofile = open(fa_out_path, "w")
for i in range(len(varnames)):
    ofile.write(">hbdx_seq_" + str(i) + "\n" +varnames[i] + "\n")
ofile.close()

