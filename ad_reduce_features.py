
import json
from datetime import datetime
import pandas as pd
import anndata as ad
from utils import UseRawExpression

# load config
config = json.load(
    open("MIRSORT_ANNOTATION_DF.json")
)

version_init =  config["assemble_anndata.py"]["version"]
version =  config["ad_reduce_features.py"]["version"]
input_path_1 = config["assemble_anndata.py"]["ad_path"]
input_path_2 = config["create_seq_annotation_df.py"]["ANNOTATION_DF"]
outputpath = config["ad_reduce_features.py"]["ad_reduced_path"]


# define if var anno shall be added
add_var = True


# load anndata
MS_ad = ad.read_h5ad(input_path_1)
print('number of sequences before reduction: ' + str(MS_ad.shape[1]))


# load annotation file for UNreduced sequences (depleted for no_annotation seqs)
annot_df = pd.read_csv(input_path_2, index_col=0)
annot_df


# subset to overlap (should be the same as len(annot_df))
seq_overlap = list(set.intersection(set(annot_df.index),set(MS_ad.var_names)))
print('seq_overlap check: ' + str(len(seq_overlap) == len(annot_df)))
annot_df = annot_df.loc[seq_overlap,:]
MS_ad = MS_ad[:,seq_overlap]
if add_var:
    MS_ad.var = annot_df
MS_ad


# reduce also raw
MS_ad_raw = UseRawExpression(MS_ad)
MS_ad.raw = MS_ad_raw
print('number of sequences after reduction: ' + str(MS_ad.shape[1]))


# adjust uns
MS_ad.uns["last_modified"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
MS_ad.uns["release_notes"] = f"miRSort_combined__ngs-{version_init}.h5ad reduced to features with miRSort-specific subclass_name annotation"
if add_var:
    MS_ad.uns['input_paths']['var_annotation'] = input_path_2
MS_ad.uns["version"] = version



# write reduced ad to file
MS_ad.write_h5ad(outputpath)

