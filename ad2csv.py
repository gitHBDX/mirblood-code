#%%
import pandas as pd
import numpy as np
import anndata as ad
from anndata import AnnData
from sklearn.manifold import TSNE
from scipy.stats import pearsonr
import scanpy as sc
import json
import omegaconf
from scipy.stats import mode
from joblib import Parallel, delayed, parallel_backend
import tqdm
from tqdm import tqdm

#%%
# define if plots should be created (requires installation of plotly in environment)
make_plots = False

if make_plots:
    import plotly
    import plotly.express as px
    import plotly.io as pio
    pio.templates.default = "simple_white"
    import dash_bio as dashbio

# load config
config = json.load(
    open("MIRSORT_ANNOTATION_DF.json")
)

ad_path = config["ad_reduce_features_further.py"]["ad_aggregated_bloodcomponents_path"]
ad_extra_path = config["ad_reduce_features_further.py"]["ad_aggregated_wholeblood_path"]

tsne_path = config["ad2csv.py"]["tsne_path"]
obs_path = config["ad2csv.py"]["obs_path"]
GEO_smp_record_info = config["ad2csv.py"]["GEO_smp_record_info"]
suppl_table_1_path = config["ad2csv.py"]["suppl_table_1_path"]
suppl_table_2_path = config["ad2csv.py"]["suppl_table_2_path"]
suppl_table_3_path = config["ad2csv.py"]["suppl_table_3_path"]
sRNAtype_path = config["ad2csv.py"]["sRNAtype_path"]
sRNA_subtype_path = config["ad2csv.py"]["sRNA_subtype_path"]
expr_path = config["ad2csv.py"]["expr_path"]
expr_extra_path = config["ad2csv.py"]["expr_extra_path"]
box_path = config["ad2csv.py"]["box_path"]
sRNA_aggreg_path = config["ad2csv.py"]["sRNA_aggreg_path"]
expr_mean_path = config["ad2csv.py"]["expr_mean_path"]
marker_path = config["ad2csv.py"]["marker_path"]
cell_prop_path = config["ad2csv.py"]["cell_prop_path"]
RPM_mean_path = config["ad2csv.py"]["RPM_mean_path"]
scaled_prop_path = config["ad2csv.py"]["scaled_prop_path"]
corr_path = config["ad2csv.py"]["corr_path"]
tsne_plot_all_path = config["ad2csv.py"]["tsne_plot_all_path"]
tsne_plot_filtered_path = config["ad2csv.py"]["tsne_plot_filtered_path"]
counts_plot_path = config["ad2csv.py"]["counts_plot_path"]

#%%
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


#%%
cell_types = [
    "Plasma",
    "Erythrocytes",
    "Thrombocytes",
    "Monocytes",
    "Neutrophils",
    "Eosinophils",
    "Basophils",
    "NK cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "B cells",
]

color_list = ['rgb(169,169,169)','rgb(255,139,139)','rgb(239,231,222)','rgb(178,188,207)','rgb(108,194,74)','rgb(79,150,51)','rgb(176,233,160)','rgb(0,157,217)','rgb(0,91,153)','rgb(0,118,163)','rgb(0,112,192)']
celltype_lut = dict(zip(cell_types, color_list))

sRNAtype_lut = {
    'miRNA': '#6cc24a',
    'rRNA': '#005b99',
    'YRNA': '#009dd9',
    'tRNA': '#0076a1',
    'lncRNA': '#8a9ab6',
    'snoRNA': '#526483',
    'snRNA': '#ccd3df',
    'piRNA': '#edf4fd',   
}

cell_counts_prev = [
    'Erythrocytes',
    'Thrombocytes',
    'Monocytes',
    'Neutrophils',
    'Eosinophils',
    'Basophils',
    'NK_cells',
    'CD4',
    'CD8',
    'B_cells',
    'Monocytes__intermediate',
    'Monocytes__classical',
    'Monocytes__non-classical',
    'NK_cells__Cd56dim',
    'NK_cells__Cd56high',
    'NK_cells__Cd56dim_mature',
    'CD4__TEMRA',
    'CD4__naive',
    'CD4__TCM',
    'CD4__TEM',
    'CD4__Th17',
    'CD4__Th1',
    'CD4__Th2',
    'CD4__Tregs',
    'CD8__TEMRA',
    'CD8__naive',
    'CD8__TCM',
    'CD8__TEM',
]

cell_counts = [
    'Erythrocytes',
    'Thrombocytes',
    'Monocytes',
    'Neutrophils',
    'Eosinophils',
    'Basophils',
    'NK cells',
    'CD4+ T cells',
    'CD8+ T cells',
    'B cells',
    'Monocytes intermediate',
    'Monocytes classical',
    'Monocytes non-classical',
    'NK cells CD56dim',
    'NK cells CD56high',
    'NK cells CD56dim mature',
    'CD4+ TEMRA',
    'CD4+ naive',
    'CD4+ TCM',
    'CD4+ TEM',
    'CD4+ Th17',
    'CD4+ Th1',
    'CD4+ Th2',
    'CD4+ Tregs',
    'CD8+ TEMRA',
    'CD8+ naive',
    'CD8+ TCM',
    'CD8+ TEM',
]

col_rename_dict = dict(zip(cell_counts_prev, cell_counts))

extended_color_list = ['rgb(255,139,139)','rgb(239,231,222)','rgb(178,188,207)','rgb(108,194,74)','rgb(79,150,51)','rgb(176,233,160)','rgb(0,157,217)','rgb(0,91,153)','rgb(0,118,163)','rgb(0,112,192)','rgb(178,188,207)','rgb(178,188,207)','rgb(178,188,207)','rgb(0,157,217)','rgb(0,157,217)','rgb(0,157,217)','rgb(0,91,153)','rgb(0,91,153)','rgb(0,91,153)','rgb(0,91,153)','rgb(0,91,153)','rgb(0,91,153)','rgb(0,91,153)','rgb(0,91,153)','rgb(0,118,163)','rgb(0,118,163)','rgb(0,118,163)','rgb(0,118,163)']
cellcounts_lut = dict(zip(cell_counts, extended_color_list))

obs_to_keep = [
    'Patient',
    'miRBlood_ID',
    'sequencer',
    'total_read_counts_before_preprocessing',
    'total_count_after_preprocessing',
    'total_count_ratio__after_to_before',
    'total_count',
    'RIN',
    'QC_RNA_ratio',
    'Lab_RNA_amount',
    'Gender',
    'Age',
    'Age_binned_5',
    'Age_binned_10',
    'mirsort_group',
    'cellsorting__cell_counts',
    'cellsorting__elution_volume',
    'cellsorting__small_concentration',
    'cellsorting__Purity',
    'blood_count_microliter',
    'Erythrocytes',
    'Thrombocytes',
    'Monocytes',
    'Neutrophils',
    'Eosinophils',
    'Basophils',
    'NK cells',
    'CD4+ T cells',
    'CD8+ T cells',
    'B cells',
    'Leukocytes',
    'Lymphocytes',
    'Monocytes_clinical',
    'T_cells',
    'Monocytes intermediate',
    'Monocytes classical',
    'Monocytes non-classical',
    'NK cells CD56dim',
    'NK cells CD56high',
    'NK cells CD56dim mature',
    'CD4+ TEMRA',
    'CD4+ naive',
    'CD4+ TCM',
    'CD4+ TEM',
    'CD4+ Th17',
    'CD4+ Th1',
    'CD4+ Th2',
    'CD4+ Tregs',
    'CD8+ TEMRA',
    'CD8+ naive',
    'CD8+ TCM',
    'CD8+ TEM'
]


def streamline_obs(anndata):

    # rename sequencer
    anndata.obs.sequencer.cat.add_categories(['NextSeq_2000'], inplace=True)
    anndata.obs.loc[anndata.obs.sequencer.isin(['VH00530']),'sequencer'] = 'NextSeq_2000'
    anndata.obs.sequencer.cat.remove_unused_categories(inplace=True)

    # rename columns
    anndata.obs.rename(columns=col_rename_dict, inplace=True)

    return anndata


# define sRNA and signature list of interest
name_OI = 'rRNA-28S_bin-162'
sRNAs_OI = ['miR-2115-3p', 'miR-224-5p', 'miR-218-5p', 'miR-4676-3p', 'miR-6503-5p']


#%%
# load anndata of blood components
MS_ad = ad.read_h5ad(ad_path)
MS_ad

#%%
# rename blood components and subpopulations
MS_ad.obs.mirsort_group.cat.add_categories(['Plasma', 'NK cells', 'CD4+ T cells', 'CD8+ T cells', 'B cells'], inplace=True)
MS_ad.obs.loc[MS_ad.obs.mirsort_group.isin(['Plasma_Norgen']),'mirsort_group'] = 'Plasma'
MS_ad.obs.loc[MS_ad.obs.mirsort_group.isin(['NK_cells']),'mirsort_group'] = 'NK cells'
MS_ad.obs.loc[MS_ad.obs.mirsort_group.isin(['CD4']),'mirsort_group'] = 'CD4+ T cells'
MS_ad.obs.loc[MS_ad.obs.mirsort_group.isin(['CD8']),'mirsort_group'] = 'CD8+ T cells'
MS_ad.obs.loc[MS_ad.obs.mirsort_group.isin(['B_cells']),'mirsort_group'] = 'B cells'
MS_ad.obs.mirsort_group.cat.remove_unused_categories(inplace=True)
MS_ad.obs.mirsort_group.unique()

#%%
# streamline obs
MS_ad = streamline_obs(MS_ad)
MS_ad.obs

#%%
# load anndata of unblocked whole blood + PAXgene
MS_extra_ad = ad.read_h5ad(ad_extra_path)
MS_extra_ad

#%%
# rename  whole blood + PAXgene
MS_extra_ad.obs.mirsort_group.cat.add_categories(['Whole Blood'], inplace=True)
MS_extra_ad.obs.loc[MS_extra_ad.obs.mirsort_group.isin(['EDTA_WholeBlood']),'mirsort_group'] = 'Whole Blood'
MS_extra_ad.obs.mirsort_group.cat.remove_unused_categories(inplace=True)
MS_extra_ad.obs.mirsort_group.unique()

#%%
# streamline obs
MS_extra_ad = streamline_obs(MS_extra_ad)
MS_extra_ad.obs


#%%
# write data for supplementary table 1 to file
MS_ad.obs[['Patient'] + cell_counts].drop_duplicates().reset_index(drop=True).to_csv(suppl_table_1_path, index=False)

#%%
# write data for supplementary table 2 to file
obs_subs = MS_ad.obs.copy()
obs_subs[['miRBlood_ID', 'mirsort_group', 'cellsorting__cell_counts', 'cellsorting__elution_volume', 'cellsorting__small_concentration','cellsorting__Purity']].reset_index(drop=True).to_csv(suppl_table_2_path, index=False)

#%%
# write info for GEO sample record
GEO_smp_record = pd.concat([MS_ad.obs, MS_extra_ad.obs],axis=0)
GEO_smp_record[['mirsort_group','sequencer']].to_csv(GEO_smp_record_info)

#%%
################################################################################################################################################################
# tSNE DATA

#%%
# quality control with purity filter

# set 'cellsorting__Purity' of Plasma to 1
MS_ad.obs.loc[MS_ad.obs.mirsort_group.str.match('Plasma'),'cellsorting__Purity'] = 1
# apply filter
MS_ad = MS_ad[MS_ad.obs.cellsorting__Purity > 0.7,:]
MS_ad

#%%
# quality control with tSNE
tsne = TSNE(n_components=2, random_state=0)
projections = tsne.fit_transform(MS_ad.X)
tsne_df = pd.concat([pd.DataFrame(projections, index=MS_ad.obs_names, columns=['tSNE vector 1', 'tSNE vector 2']),MS_ad.obs.mirsort_group],axis=1)

if make_plots:
    fig_1 = px.scatter(
        tsne_df, x='tSNE vector 1', y='tSNE vector 2',
        color='mirsort_group',
        color_discrete_map=celltype_lut,
        hover_name=tsne_df.index
    )
    fig_1.update_layout(legend_title_text='cell type', font_family='Arial').show()
    fig_1.write_image(tsne_plot_all_path, width=600, height=600)
    fig_1.write_image(tsne_plot_all_path, width=600, height=600)


#%%
# apply tSNE filter
tSNE_outliers = [
    'Neutrophils__P51',
    'Plasma_Norgen__P02',
    'Basophils__P37',
    'B_cells__P38',
    ]

MS_ad = MS_ad[~MS_ad.obs_names.isin(tSNE_outliers),:] 
MS_ad


#%%
# quality control with tSNE
from sklearn.manifold import TSNE

tsne = TSNE(n_components=2, random_state=0)
projections = tsne.fit_transform(MS_ad.X)
tsne_df = pd.concat([pd.DataFrame(projections, index=MS_ad.obs_names, columns=['tSNE vector 1', 'tSNE vector 2']),MS_ad.obs.mirsort_group],axis=1)

if make_plots:
    fig_2 = px.scatter(
        tsne_df, x='tSNE vector 1', y='tSNE vector 2',
        color='mirsort_group',
        color_discrete_map=celltype_lut,
        hover_name=tsne_df.index
    )
    fig_2.update_layout(legend_title_text='cell type', font_family='Arial').show()
    fig_2.write_image(tsne_plot_filtered_path, width=600, height=600)
    fig_2.write_image(tsne_plot_filtered_path, width=600, height=600)

#%%
# write tSNE to file
tsne_df.to_csv(tsne_path)

#%%
################################################################################################################################################################
# OBS DATA
#%%
# write obs to file
obs_df = MS_ad.obs
obs_df[obs_to_keep].to_csv(obs_path)

#%%
# barplot of patient data (age and sex)
if make_plots:
    patients = px.histogram(obs_df.drop_duplicates(subset=["miRBlood_ID"]), x='Age', color='Gender', barmode='group')
    patients

#%%
# boxplot of total_count_after_preprocessing / RIN / total_count_ratio__after_to_before / QC_RNA_ratio / Lab_RNA_amount
if make_plots:
    fig_3 = px.box(MS_ad.obs, x='mirsort_group', y='total_count_after_preprocessing', color='mirsort_group', color_discrete_map=celltype_lut)
    fig_3.update_layout(xaxis_title="", yaxis_title="total read counts after preprocessing", showlegend=False, font_family='Arial')
    fig_3.show()
    fig_3.write_image(counts_plot_path, width=1200, height=600)
    fig_3.write_image(counts_plot_path, width=1200, height=600)

#%%
################################################################################################################################################################
# SRNA TYPE DATA
#%%
sRNAtype_df = MS_ad.var
sRNAtype_df

#%%
# write to file 
sRNAtype_df.to_csv(sRNAtype_path)

#%%
# plot pie chart of sRNA type attribution of signature
sRNA_type_label = sRNAtype_df.loc[sRNAs_OI,'small_RNA_class_annotation']
sRNA_type_label = sRNA_type_label.cat.remove_unused_categories()
if make_plots:
    signature_sRNAannot = px.pie(sRNA_type_label, values=sRNA_type_label.value_counts().values, names=sRNA_type_label.value_counts().index, color=sRNA_type_label.value_counts().index, hole = 0.55, width=500, height=500, color_discrete_map=sRNAtype_lut)
    signature_sRNAannot.update_layout(showlegend=False).update_traces(textposition="inside", textinfo="percent+label", textfont_size=15)
    signature_sRNAannot

#%%
################################################################################################################################################################
# SRNA SUBTYPE DATA
#%%
sRNA_subtypes_df = sRNAtype_df.reset_index()
sRNA_subtypes_df.loc[sRNA_subtypes_df.small_RNA_class_annotation == 'rRNA','sRNA_subclass'] = sRNA_subtypes_df.subclass_name.str.split('_').str[0]
sRNA_subtypes_df.loc[sRNA_subtypes_df.small_RNA_class_annotation == 'YRNA','sRNA_subclass'] = sRNA_subtypes_df.subclass_name.str.split('_').str[0]
sRNA_subtypes_df.loc[sRNA_subtypes_df.small_RNA_class_annotation == 'tRNA','sRNA_subclass'] = sRNA_subtypes_df.subclass_name.str.split('__').str[0]
sRNA_subtypes_df.loc[sRNA_subtypes_df.small_RNA_class_annotation == 'piRNA','sRNA_subclass'] = 'piRNA'
sRNA_subtypes_df.loc[sRNA_subtypes_df.small_RNA_class_annotation == 'lncRNA','sRNA_subclass'] = 'lncRNA'
sRNA_subtypes_df.loc[sRNA_subtypes_df.subclass_name.str.contains('SNORA'),'sRNA_subclass'] = 'H/ACA Box'
sRNA_subtypes_df.loc[sRNA_subtypes_df.subclass_name.str.contains('SNORD'),'sRNA_subclass'] = 'C/D Box'
sRNA_subtypes_df.loc[((sRNA_subtypes_df.small_RNA_class_annotation == 'snoRNA') & (~sRNA_subtypes_df.subclass_name.str.contains(r'SNOR[A|D]'))),'sRNA_subclass'] = 'misc'
sRNA_subtypes_df.loc[sRNA_subtypes_df.small_RNA_class_annotation == 'snRNA','sRNA_subclass'] = sRNA_subtypes_df.subclass_name.str.split('-').str[0]
sRNA_subtypes_df.loc[((sRNA_subtypes_df.small_RNA_class_annotation == 'miRNA') & (~sRNA_subtypes_df.subclass_name.str.contains(r'-[5|3]p'))),'sRNA_subclass'] = 'misc'
sRNA_subtypes_df.loc[sRNA_subtypes_df.subclass_name.str.contains('-5p'),'sRNA_subclass'] = '5p'
sRNA_subtypes_df.loc[sRNA_subtypes_df.subclass_name.str.contains('-3p'),'sRNA_subclass'] = '3p'
sRNA_subtypes_df = pd.DataFrame(sRNA_subtypes_df.groupby('small_RNA_class_annotation').sRNA_subclass.value_counts())
sRNA_subtypes_df


#%%
# write to file 
sRNA_subtypes_df.to_csv(sRNA_subtype_path)

#%%
################################################################################################################################################################
# EXPRESSION DATA (for download)
expr_df = MS_ad.to_df()
expr_df.to_csv(expr_path)

#%%
################################################################################################################################################################
# EXTRA EXPRESSION DATA (for download)
expr_extra_df = MS_extra_ad.to_df()
expr_extra_df.to_csv(expr_extra_path)



#%%
################################################################################################################################################################
# CELLTYPE MARKER DATA
#%%
sc.tl.rank_genes_groups(MS_ad, 'mirsort_group', use_raw=False, method='wilcoxon', pts=True)
sc.pl.rank_genes_groups_dotplot(MS_ad, n_genes=5, standard_scale='var')

#%%
# get ordered marker list (positive log2FC and 100% expression in blood component)

marker_dict = {}

for celltype in cell_types:
    markers = sc.get.rank_genes_groups_df(MS_ad, group=celltype).set_index('names')
    markers = markers.loc[((markers.pvals_adj < 0.05) & (markers.logfoldchanges > 0) & (markers.pct_nz_group == 1)),:].index
    marker_dict.update({celltype: list(markers)})

marker_dict.keys()

#%%
for celltype in cell_types:
    print(celltype)
    print(len(marker_dict[celltype]))

#%%
# save to file
with open(marker_path, "w") as outfile:
    json.dump(marker_dict, outfile)


#%%
################################################################################################################################################################
# BOXPLOT DATA
#%%
# calculate boxplot data for log2(RPM+1) plot
box_df = pd.concat([MS_ad.to_df(), MS_ad.obs.mirsort_group], axis=1)
box_df = box_df.groupby(by='mirsort_group').describe()
box_df = box_df.loc[cell_types,:]
box_df

#%%
# write to file
box_df.to_csv(box_path)

#%%
# plot boxplot for sRNA of interest
plot_df = pd.melt(box_df[name_OI][['min','25%','50%','75%','max']].T)
if make_plots:
    single_box = px.box(plot_df, x='mirsort_group', y='value', color='mirsort_group', color_discrete_map=celltype_lut).update_layout(xaxis_title="", yaxis_title="log2(RPM+1)", showlegend=False)
    single_box

#%%
################################################################################################################################################################
# SRNA CLASS AGGREGATED DATA
#%%
MS_ad_sRNAclass_aggreg = AnnDataAggregator(MS_ad, aggregate_obs=False, by='small_RNA_class_annotation', work_on_raw=True)
MS_ad_sRNAclass_aggreg

#%%
# calculate RPM from raw counts
MS_ad_sRNAclass_aggreg.X = (MS_ad_sRNAclass_aggreg.raw.X.T / MS_ad_sRNAclass_aggreg.raw.X.T.sum(axis=0) * 1e6).T

# calculate mean RPM
res_l_mean = []

for i in MS_ad_sRNAclass_aggreg.obs.mirsort_group.unique():
    res_l_mean.append(
        [
            np.mean(MS_ad_sRNAclass_aggreg[MS_ad_sRNAclass_aggreg.obs.mirsort_group.isin([i]), :].X, axis=0)
        ]
    )

res_l_mean = np.array(res_l_mean)[:, 0, :]

res_l_mean = pd.DataFrame(
    res_l_mean,
    columns=MS_ad_sRNAclass_aggreg.var_names,
    index=MS_ad_sRNAclass_aggreg.obs.mirsort_group.unique()
)

res_l_mean = res_l_mean.loc[cell_types,:]
res_l_mean

#%%
# write to file
res_l_mean.to_csv(sRNA_aggreg_path)

#%%
# plot proportional sRNA class (RPM) expression per blood component (heatmap)
if make_plots:
    px.imshow(res_l_mean / res_l_mean.sum(axis=1)[:,None], color_continuous_scale="blues")

#%%
# plot RPM proportion for a given blood component (pie)
celltype_OI = 'B cells'

pie_data = res_l_mean.copy()

# subset to celltype
pie_data = pie_data.loc[celltype_OI,:]

# calculate fraction
pie_data = pie_data / pie_data.sum()

# combine all sRNAs < 2% as others
others = pd.Series(data=pie_data[pie_data < 0.02].sum(), index=['others (<2%)'])
pie_data = pie_data[pie_data > 0.02]
pie_data = pd.concat([pie_data, others],axis=0)
if make_plots:
    px.pie(pie_data, values=pie_data.values, names=pie_data.index, color_discrete_sequence=extended_color_list)

#%%
################################################################################################################################################################
# PIECHART/HEATMAP DATA
#%%
# calculate log2(RPM+1) mean
res_i_mean = []

for i in MS_ad.obs.mirsort_group.unique():
    res_i_mean.append(
        [
            np.mean(MS_ad[MS_ad.obs.mirsort_group.isin([i]), :].X, axis=0)
        ]
    )

res_i_mean = np.array(res_i_mean)[:, 0, :]

res_i_mean = pd.DataFrame(
    res_i_mean,
    columns=MS_ad.var_names,
    index=MS_ad.obs.mirsort_group.unique()
)

res_i_mean = res_i_mean.loc[cell_types,:]
res_i_mean

#%%
# write to file
res_i_mean.to_csv(expr_mean_path)

#%%
# MARKER: plot log2(RPM+1) expression per celltype (top2_union only)
if make_plots:
    top2_union = []
    for celltype in marker_dict.keys():
        top2_union = top2_union + marker_dict[celltype][0:2]

    top2_union = list(set(top2_union))

    heat_df = res_i_mean.loc[:,top2_union].copy().dropna().T
    columns = list(heat_df.columns)
    rows = list(heat_df.index)
    top2_union_plot = dashbio.Clustergram(
        data=heat_df.values,
        column_labels=columns,
        row_labels=rows,
        center_values=False,
        color_threshold={"row": 250, "col": 700},
        cluster="row",
        height=800,
        width=700,
        color_map="blues",
    )
    top2_union_plot.update_layout({'paper_bgcolor': 'rgba(255, 255, 255, 255)'})#.update_traces(dict(showscale=False, coloraxis=None), selector={"type": "heatmap"})


#%%
# TOPEXPRESSED: plot log2(RPM+1) expression per celltype (top20_union only)
if make_plots:
    top20_union = []
    for celltype in res_i_mean.index:
        top20_union = top20_union + list(res_i_mean.loc[celltype,:].sort_values(ascending=False)[0:20].index)

    top20_union = list(set(top20_union))

    heat_df = res_i_mean.loc[:,top20_union].copy().dropna().T
    columns = list(heat_df.columns)
    rows = list(heat_df.index)
    top20_union_plot = dashbio.Clustergram(
        data=heat_df.values,
        column_labels=columns,
        row_labels=rows,
        center_values=False,
        color_threshold={"row": 250, "col": 700},
        cluster="row",
        height=800,
        width=700,
        color_map="blues",
    )
    top20_union_plot.update_layout({'paper_bgcolor': 'rgba(255, 255, 255, 255)'})#.update_traces(dict(showscale=False, coloraxis=None), selector={"type": "heatmap"})

#%%
# calculate RPM from raw counts
MS_ad.X = (MS_ad.raw.X.T / MS_ad.raw.X.T.sum(axis=0) * 1e6).T

# calculate RPM proportion
res_j_mean = []
res_j_std = []

for i in MS_ad.obs.mirsort_group.unique():
    res_j_mean.append(
        [
            np.mean(MS_ad[MS_ad.obs.mirsort_group.isin([i]), :].X, axis=0)
        ]
    )
    res_j_std.append(
        [
            np.nanstd(MS_ad[MS_ad.obs.mirsort_group.isin([i]), :].X, axis=0)
        ]
    )

res_j_mean = np.array(res_j_mean)[:, 0, :]
res_j_std = np.array(res_j_std)[:, 0, :]

res_j_mean = pd.DataFrame(
    res_j_mean,
    columns=MS_ad.var_names,
    index=MS_ad.obs.mirsort_group.unique()
)
res_j_std = pd.DataFrame(
    res_j_std,
    columns=MS_ad.var_names,
    index=MS_ad.obs.mirsort_group.unique()
)

res_j_prop = res_j_mean / res_j_mean.sum(axis=0)
res_j_prop.index = MS_ad.obs.mirsort_group.unique()


res_j_prop = res_j_prop.loc[cell_types,:]
res_j_prop

#%%
# calculate prop per celltype
cell_prop = res_j_mean / res_j_mean.sum(axis=1)[:,None]
cell_prop

#%%
# write to file
res_j_mean.to_csv(RPM_mean_path)
cell_prop.to_csv(cell_prop_path)

#%%
# plot RPM proportion per celltype
classes_OI = ['miRNA']
celltype_OI = 'B cells'

# OPTIONAL: subset to miRNAs
sRNAtype_df = pd.read_csv(sRNAtype_path, index_col=0)
sRNA_selected = sRNAtype_df.loc[sRNAtype_df.small_RNA_class_annotation.isin(classes_OI),:].index
pie_data = cell_prop[sRNA_selected]

# subset to celltype
pie_data = pie_data.loc[celltype_OI,:]

# calculate fraction
pie_data = pie_data / pie_data.sum()

# combine all sRNAs < 2% as others
others = pd.Series(data=pie_data[pie_data < 0.02].sum(), index=['others (<2%)'])
pie_data = pie_data[pie_data > 0.02]
pie_data = pd.concat([pie_data, others],axis=0)
if make_plots:
    cell_proportion = px.pie(pie_data, values=pie_data.values, names=pie_data.index, color_discrete_sequence=extended_color_list)
    cell_proportion


#%%
# calculate scaled RPM proportion
# set 'cellsorting__cell_counts' of Plasma to 1
MS_ad.obs.loc[MS_ad.obs.mirsort_group.str.match('Plasma'),'cellsorting__cell_counts'] = 1
MS_ad.obs.cellsorting__cell_counts

# set 'blood_count_microliter' of Plasma to 0.5
MS_ad.obs.loc[MS_ad.obs.mirsort_group.str.match('Plasma'),'blood_count_microliter'] = 0.5
MS_ad.obs.blood_count_microliter

# set cell count of 0 to nan (otherwise fails to calculate this sample)
MS_ad.obs.cellsorting__cell_counts[MS_ad.obs.cellsorting__cell_counts == 0] = np.nan


# calculate total small RNA eluted [pg]
MS_ad.obs['sRNA_mass_eluted'] = np.multiply(MS_ad.obs.cellsorting__small_concentration, MS_ad.obs.cellsorting__elution_volume)

# calculate small RNA mass per sorted cell [pg/cell]
MS_ad.obs['sRNA_mass_per_cells'] = np.divide(MS_ad.obs.sRNA_mass_eluted, MS_ad.obs.cellsorting__cell_counts)

# calculate sRNA mass of celltype per bloodvolume [pg/µl]
MS_ad.obs['sRNA_mass_per_bloodvol'] = np.multiply(MS_ad.obs.sRNA_mass_per_cells, MS_ad.obs.blood_count_microliter)


# create new adata with scaled RPM values
MS_ad_scaled = ad.AnnData(
    X=MS_ad.X * MS_ad.obs.sRNA_mass_per_bloodvol[:, np.newaxis],
    obs=MS_ad.obs,
    var=MS_ad.var,
)
MS_ad_scaled

# create dataframe with mean, std and prop of expression per cell type
res_k_mean = []
res_k_std = []

for i in MS_ad_scaled.obs.mirsort_group.unique():
    res_k_mean.append(
        [
            np.nanmean(MS_ad_scaled[MS_ad_scaled.obs.mirsort_group.isin([i]), :].X, axis=0)
        ]
    )
    res_k_std.append(
        [
            np.nanstd(MS_ad_scaled[MS_ad_scaled.obs.mirsort_group.isin([i]), :].X, axis=0)
        ]
    )

res_k_mean = np.array(res_k_mean)[:, 0, :]
res_k_std = np.array(res_k_std)[:, 0, :]

res_k_mean = pd.DataFrame(
    res_k_mean,
    columns=MS_ad_scaled.var_names,
    index=MS_ad_scaled.obs.mirsort_group.unique()
)
res_k_std = pd.DataFrame(
    res_k_std,
    columns=MS_ad_scaled.var_names,
    index=MS_ad_scaled.obs.mirsort_group.unique()
)

res_k_prop = res_k_mean / res_k_mean.sum(axis=0)
res_k_prop.index = MS_ad_scaled.obs.mirsort_group.unique()


res_k_prop = res_k_prop.loc[cell_types,:]
res_k_prop

#%%
# write to file
res_k_prop.to_csv(scaled_prop_path)
MS_ad.obs[['miRBlood_ID','mirsort_group','sRNA_mass_eluted','sRNA_mass_per_cells','sRNA_mass_per_bloodvol']].reset_index(drop=True).to_csv(suppl_table_3_path, index=False)


#%%
# plot pie of sRNA OI
if make_plots:
    single_deconvolution = px.pie(res_k_prop, values=name_OI, names=res_k_prop.index, color=res_k_prop.index, color_discrete_map=celltype_lut, hole = 0.55, width=500, height=500).update_layout(showlegend=False).update_traces(textposition="inside", textinfo="percent+label", textfont_size=15) # title="Cell type contribution to abundance in blood", title_x=0.5, 
    single_deconvolution

#%%
# plot heatmap of signature OI
if make_plots:
    heat_df = res_k_prop.copy().dropna(axis=1).T
    columns = list(heat_df)
    rows = list(sRNAs_OI)

    # clustermap of deconvolution
    signature_deconvolution = dashbio.Clustergram(
        data=heat_df.loc[rows].values,
        column_labels=columns,
        row_labels=rows,
        color_threshold={"row": 250, "col": 700},
        center_values=False,
        cluster="row",
        height=800,
        width=700,
        color_map="blues",
    )
    signature_deconvolution.update_layout({'paper_bgcolor': 'rgba(255, 255, 255, 255)'})
    signature_deconvolution.update_traces({'colorbar': dict(title="log2(RPM+1)", orientation="h", xpad = 10)}, selector={"type": "heatmap"})
    #.update_traces(dict(showscale=False, coloraxis=None), selector={"type": "heatmap"})



#%%
################################################################################################################################################################
# BLOODCOUNT CORRELATION DATA
#%%
# calculate Pearson correlation
def calc_correlation(X:pd.DataFrame, y:pd.Series, corr_fn:callable=pearsonr, get_p:bool=False):
    correlations=[ corr_fn(X.iloc[:,i], y) for i in range(X.shape[1]) ]
    r,p = zip(*correlations)
    
    correlation = pd.Series(p if get_p else r, index=X.columns)
    return correlation

def rank_by_correlation(X:pd.DataFrame, y:pd.Series, corr_fn:callable=pearsonr, rank_by_p:bool=False):
    correlations=[ corr_fn(X.iloc[:,i], y) for i in range(X.shape[1]) ]
    r,p = zip(*correlations)
    
    ranking = pd.Series(p if rank_by_p else r, index=X.columns).sort_values(ascending=True if rank_by_p else False)
    return ranking


#%%
corr_df = pd.DataFrame()

# calculate RPM from raw counts
MS_extra_ad.X = (MS_extra_ad.raw.X.T / MS_extra_ad.raw.X.T.sum(axis=0) * 1e6).T

# correlate against 'Whole Blood' expression
MS_extra_ad_subs = MS_extra_ad[MS_extra_ad.obs.mirsort_group.isin(['Whole Blood']),:]

for celltype in cell_counts:
    corr_df[celltype] = calc_correlation(MS_extra_ad_subs[MS_extra_ad_subs.obs[celltype].dropna().index,:].to_df(),MS_extra_ad_subs.obs[celltype].dropna())
corr_df

#%%
# write to file
corr_df.to_csv(corr_path)

#%%
# plot blood count correlation heatmap for signature sRNAs
if make_plots:
    signature_correlation = px.imshow(corr_df.loc[sRNAs_OI,:])#, color_continuous_scale='bluered', zmin=-0.3, zmax=0.3)
    signature_correlation

#%%
# plot ranked correlation for defined small RNA name
pearson_coef = corr_df.loc[name_OI,:].sort_values(ascending=False)
if make_plots:
    single_correlation = px.bar(pearson_coef, x=pearson_coef.index, y=pearson_coef.values, color=pearson_coef.index, color_discrete_map=cellcounts_lut).update_layout(yaxis_title="Pearson Coefficient", showlegend=False) # title='Blood count correlation', title_x=0.5, 
    single_correlation

#%%