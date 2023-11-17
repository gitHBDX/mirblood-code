import json
from datetime import datetime
import anndata as ad
from utils import UseRawExpression, FeatureThreshold


# load config
config = json.load(
    open("MIRSORT_ANNOTATION_DF.json")
)

inputpath = config["ad_aggregate.py"]["ad_aggregated_path"] + config["ad_aggregate.py"]["ad_aggregated_filename"]
version = config["ad_reduce_features.py"]["version"]
outputpath = config["ad_reduce_features_further.py"]["ad_aggregated_bloodcomponents_path"]
outputpath_extra = config["ad_reduce_features_further.py"]["ad_aggregated_wholeblood_path"]



# load anndata
MS_ad = ad.read_h5ad(inputpath)
print('number of features before reduction: ' + str(MS_ad.shape[1]))


# subset to unblocked PAXgene and WholeBlood
extra_list = ['EDTA_WholeBlood']
MS_extra_ad = MS_ad[MS_ad.obs.mirsort_group.isin(extra_list),:]
# remove unused categories
MS_extra_ad.obs['mirsort_group'].cat.remove_unused_categories(inplace=True)


# subset to cells to plasma only (+ WholeBlood)
cell_list = ['Neutrophils','CD8','NK_cells','CD4','Thrombocytes','Monocytes','B_cells','Basophils','Eosinophils','Erythrocytes','Plasma_Norgen']
MS_ad = MS_ad[MS_ad.obs.mirsort_group.isin(cell_list),:]
# remove unused categories
MS_ad.obs['mirsort_group'].cat.remove_unused_categories(inplace=True)


# subset by expression threshold (over_2)
expression_TH = 1.0986122886681096
MS_ad = FeatureThreshold(MS_ad, expression_TH)
print('number of features after reduction: ' + str(MS_ad.shape[1]))


# reduce also raw
MS_ad_raw = UseRawExpression(MS_ad)
MS_ad.raw = MS_ad_raw


# adjust uns
MS_ad.uns["last_modified"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
MS_ad.uns["release_notes"] = MS_ad.uns["release_notes"] + "; reduced by the following filters: IncludeObservation: mirsort_group = [" + ', '.join(cell_list) + "] & " + "FeatureThreshold: " + str(expression_TH)
MS_ad.uns["version"] = version

print("blood component anndata shape: ")
print(MS_ad.X.shape)

# write reduced ad to file
MS_ad.write_h5ad(outputpath)



# subset extra ad alike
MS_extra_ad = MS_extra_ad[:,MS_ad.var_names]
# reduce also raw
MS_extra_ad_raw = UseRawExpression(MS_extra_ad)
MS_extra_ad.raw = MS_extra_ad_raw


# adjust uns
MS_extra_ad.uns["last_modified"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
MS_extra_ad.uns["release_notes"] = MS_extra_ad.uns["release_notes"] + "; same features like blood component ad; mirsort_group = [" + ', '.join(cell_list) + "]"
MS_extra_ad.uns["version"] = version

print("whole blood anndata shape: ")
print(MS_extra_ad.X.shape)

# write reduced ad to file
MS_extra_ad.write_h5ad(outputpath_extra)

