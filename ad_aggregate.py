
import json
from datetime import datetime
import anndata as ad
from utils import AnnDataAggregator


# load config
config = json.load(
    open("MIRSORT_ANNOTATION_DF.json")
)

version =  config["ad_reduce_features.py"]["version"]
input_path = config["ad_aggregate.py"]["ad_seq_path"] + config["ad_aggregate.py"]["ad_seq_filename"]
output_path = config["ad_aggregate.py"]["ad_aggregated_path"] + config["ad_aggregate.py"]["ad_aggregated_filename"]



# load AnnData 
MS_ad = ad.read_h5ad(input_path)



# aggregate on subclass_name
aggregate_obs_status = False
aggregate_on = 'subclass_name'
raw_status = True
MS_ad = AnnDataAggregator(MS_ad, aggregate_obs=aggregate_obs_status,by=aggregate_on,work_on_raw=raw_status)
MS_ad.uns["aggregator_params"] = {'aggregate_obs': aggregate_obs_status, 'by': aggregate_on, 'work_on_raw': raw_status}
MS_ad.uns["last_modified"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
MS_ad.uns["release_notes"] = 'aggregated version of ' + version                   



# convert string columns to categorical
str_columns = MS_ad.obs.columns[MS_ad.obs.dtypes == 'string']
for col in str_columns:
    MS_ad.obs[col] = MS_ad.obs[col].astype("category")

print(MS_ad.X.shape)

# write to file
MS_ad.write_h5ad(output_path)


