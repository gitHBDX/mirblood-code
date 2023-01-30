#%%
import json
from datetime import datetime
import anndata as ad

# load config
config = json.load(
    open("MIRSORT_ANNOTATION_DF.json")
)

version =  config["ad_reduce_features.py"]["version"]
inputpath_1 = config["ad_reduce_features.py"]["ad_reduced_path"]
outputpath_2 = config["ad_aggregate.py"]["ad_aggregated_path"] + config["ad_aggregate.py"]["ad_aggregated_filename"]


#%%
MS_ad = ad.read_h5ad(inputpath_1)
annotation_df = MS_ad.var[['subclass_name','small_RNA_class_annotation']].drop_duplicates('subclass_name').set_index('subclass_name')
annotation_df


#%%
# add to aggregated anndata
MS_ad = ad.read_h5ad(outputpath_2)
feature_overlap = list(set.intersection(set(MS_ad.var_names),set(annotation_df.index)))
annotation_df = annotation_df.loc[feature_overlap,:]
MS_ad = MS_ad[:,feature_overlap]
MS_ad.var = annotation_df

MS_ad.uns['input_paths']['var_annotation'] = {'aggrgated_version_date': datetime.now().strftime("%Y-%m-%d")}
MS_ad.uns["last_modified"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
MS_ad.uns["release_notes"] = 'aggregated version of ' + version                   
MS_ad.write_h5ad(outputpath_2)

#%%