
import json
from datetime import datetime
import anndata as ad
from anndata import AnnData
import pandas as pd
import numpy as np
import scanpy as sc
import omegaconf
from scipy.stats import mode
from joblib import Parallel, delayed, parallel_backend
import tqdm
from tqdm import tqdm


# load config
config = json.load(
    open("MIRSORT_ANNOTATION_DF.json")
)

version =  config["ad_reduce_features.py"]["version"]
input_path = config["ad_aggregate.py"]["ad_seq_path"] + config["ad_aggregate.py"]["ad_seq_filename"]
output_path = config["ad_aggregate.py"]["ad_aggregated_path"] + config["ad_aggregate.py"]["ad_aggregated_filename"]



def UseRawExpression(ad):
    raw_ad = ad.copy()
    raw_ad.X = raw_ad.raw[ad.obs_names, ad.var_names].X
    return raw_ad

def RPMNormalize(ad,use_log1p):
    ad_log = ad.copy()
    np.divide(ad_log.X, ad_log.obs["total_count"][:, None] / 1e6, out=ad_log.X)
    if use_log1p:
        sc.pp.log1p(ad_log, base=2)
    ad_log.raw = ad
    return ad_log

def AnnDataAggregator(ad,aggregate_obs,by,work_on_raw):
    def _join_mapped(x):
        return "__".join([str(i) for i in x])

    # the rows are always aggregated, so aggregated annotations have to oriented along first axis (rows)
    anno_df, other_df = (ad.obs.copy(), ad.var.copy()) if aggregate_obs else (ad.var.copy(), ad.obs.copy())

    if isinstance(by, (list, omegaconf.listconfig.ListConfig)):
        anno_df.dropna(subset=by, inplace=True)
        new_temp_column = "__".join(by)
        anno_df[new_temp_column] = anno_df.loc[:, by].agg(_join_mapped, axis=1)
        temp_col = new_temp_column
    else:
        anno_df.dropna(subset=[by], inplace=True)
        temp_col = by

    categories = anno_df[temp_col].unique()

    X_df = ad.to_df() if aggregate_obs else ad.to_df().T

    X_df = X_df.loc[anno_df.index, :]

    raw_df = pd.DataFrame(
        ad.raw[ad.obs_names, ad.var_names].X,
        index=ad.obs_names,
        columns=ad.var_names,
    )
    if not aggregate_obs:
        raw_df = raw_df.T
    raw_df = raw_df.loc[anno_df.index, :]

    # break data into chunks where each chunk is a different category to be aggregated
    X_chunks = []
    raw_chunks = []
    anno_chunks = []
    for category in categories:
        anno_chunk = anno_df[anno_df[temp_col] == category].copy()
        anno_chunks.append(anno_chunk)
        X_chunks.append(X_df.loc[anno_chunk.index, :].copy())
        raw_chunks.append(raw_df.loc[anno_chunk.index, :].copy())

        assert raw_chunks[-1].shape == X_chunks[-1].shape
        assert len(anno_chunks[-1]) == len(raw_chunks[-1]) == len(X_chunks[-1])

    del X_df
    del raw_df
    del anno_df

    def _aggregate_df(df, category, agg_function=None):
        def _default_agg_fn(x):
            if str(x.dtype) in ["object", "string", "str", "category", "bool"]:
                try:
                    m = mode(x, nan_policy="omit")
                    aggregated = m.mode[0]
                except TypeError as t:
                    aggregated = None

            else:
                aggregated = np.mean(x)

            return aggregated

        a = pd.Series(
            df.agg(
                _default_agg_fn if agg_function is None else agg_function,
                axis="rows",
            ),
            name=category,
        )
        return a

    with parallel_backend("loky", n_jobs=10):
        aggregated_list_X = Parallel()(
            delayed(_aggregate_df)(chunk, category, np.mean)
            for chunk, category in tqdm(
                zip(X_chunks, categories),
                desc=f"Aggregating AnnData expression over {'.obs' if aggregate_obs else '.var'} column {temp_col} ",
            )
        )

        aggregated_list_raw = Parallel()(
            delayed(_aggregate_df)(raw_chunk, category, np.sum)
            for raw_chunk, category in tqdm(
                zip(raw_chunks, categories),
                desc=f"Aggregating AnnData raw expression over {'.obs' if aggregate_obs else '.var'} column {temp_col} ",
            )
        )

    new_X_df = pd.DataFrame(aggregated_list_X)

    new_raw_df = pd.DataFrame(aggregated_list_raw)

    assert new_X_df.shape == new_raw_df.shape


    new_anno_df = pd.DataFrame(index=new_X_df.index)

    if work_on_raw:

        new_ad = AnnData(X=new_raw_df.values, obs=new_anno_df, var=other_df, uns=ad.uns)
        new_ad.raw = new_ad

    else:
        new_ad = AnnData(X=new_X_df.values, obs=new_anno_df, var=other_df, uns=ad.uns)
        new_ad.raw = AnnData(X=new_raw_df.values, obs=new_anno_df, var=other_df, uns=ad.uns)

    if not aggregate_obs:
        new_ad = new_ad.T

    if work_on_raw:
        new_ad = RPMNormalize(new_ad,use_log1p=True)
    return new_ad




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


