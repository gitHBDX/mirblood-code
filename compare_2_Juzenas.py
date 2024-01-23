#%%
import pandas as pd
import numpy as np
import anndata as ad
from scipy.stats import pearsonr
import scanpy as sc
import json

import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

import plotly
import plotly.express as px
import plotly.io as pio
pio.templates.default = "simple_white"

from upsetplot import from_contents
from upsetplot import UpSet

from utils import AnnDataAggregator

# load config
config = json.load(
    open("MIRSORT_ANNOTATION_DF.json")
)


output_path = config["compare_2_Juzenas.py"]["output_folder"]
sc.settings.figdir = config["compare_2_Juzenas.py"]["output_folder"]

MS_ad_path = config["ad_reduce_features_further.py"]["ad_aggregated_bloodcomponents_path"]
Kiel_ad_path = config["compare_2_Juzenas.py"]["Kiel_ad_path"]


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
    #"EDTA_WholeBlood"
]

color_list = ['rgb(169,169,169)','rgb(255,139,139)','rgb(239,231,222)','rgb(178,188,207)','rgb(108,194,74)','rgb(79,150,51)','rgb(176,233,160)','rgb(0,157,217)','rgb(0,91,153)','rgb(0,118,163)','rgb(0,112,192)', 'rgb(0,0,0)']
celltype_lut = dict(zip(cell_types, color_list))

cell_types_b = [
    "Serum",
    "CD235a",
    "tbd_1",
    "CD14",
    "CD15",
    "tbd_2",
    "tbd_3",
    "CD56",
    "CD4",
    "CD8",
    "CD19",
    #"WB"
]
celltype_lut_b = dict(zip(cell_types_b, color_list))


#%%
# dictionary of group names
group_name_dict = {
  'Neutrophils': {'Kiel': 'CD15', 'miR-Blood': 'Neutrophils'},
  'Monocytes': {'Kiel': 'CD14', 'miR-Blood': 'Monocytes'},
  'NK cells': {'Kiel': 'CD56', 'miR-Blood': 'NK cells'},
  'CD4': {'Kiel': 'CD4', 'miR-Blood': 'CD4+ T cells'},
  'CD8': {'Kiel': 'CD8', 'miR-Blood': 'CD8+ T cells'},
  'B cells': {'Kiel': 'CD19', 'miR-Blood': 'B cells'},
  'Erythrocytes': {'Kiel': 'CD235a', 'miR-Blood': 'Erythrocytes'},
  'Plasma': {'Kiel': 'Serum', 'miR-Blood': 'Plasma'},
  #'EDTA_WholeBlood': {'Kiel': 'WB', 'miR-Blood': 'EDTA_WholeBlood'},
}


#%%
# load aggregated miR-Blood data
MS_sc_aggreg_ad = ad.read_h5ad(MS_ad_path)
# rename blood components and subpopulations
MS_sc_aggreg_ad.obs.mirsort_group.cat.add_categories(['Plasma', 'NK cells', 'CD4+ T cells', 'CD8+ T cells', 'B cells'], inplace=True)
MS_sc_aggreg_ad.obs.loc[MS_sc_aggreg_ad.obs.mirsort_group.isin(['Plasma_Norgen']),'mirsort_group'] = 'Plasma'
MS_sc_aggreg_ad.obs.loc[MS_sc_aggreg_ad.obs.mirsort_group.isin(['NK_cells']),'mirsort_group'] = 'NK cells'
MS_sc_aggreg_ad.obs.loc[MS_sc_aggreg_ad.obs.mirsort_group.isin(['CD4']),'mirsort_group'] = 'CD4+ T cells'
MS_sc_aggreg_ad.obs.loc[MS_sc_aggreg_ad.obs.mirsort_group.isin(['CD8']),'mirsort_group'] = 'CD8+ T cells'
MS_sc_aggreg_ad.obs.loc[MS_sc_aggreg_ad.obs.mirsort_group.isin(['B_cells']),'mirsort_group'] = 'B cells'
MS_sc_aggreg_ad.obs.mirsort_group.cat.remove_unused_categories(inplace=True)
MS_sc_aggreg_ad.obs.mirsort_group.unique()
MS_sc_aggreg_ad

#%%
# quality control with purity filter

# set 'cellsorting__Purity' of Plasma to 1
MS_sc_aggreg_ad.obs.loc[MS_sc_aggreg_ad.obs.mirsort_group.str.match('Plasma'),'cellsorting__Purity'] = 1
# apply filter
MS_sc_aggreg_ad = MS_sc_aggreg_ad[MS_sc_aggreg_ad.obs.cellsorting__Purity > 0.7,:]
MS_sc_aggreg_ad

#%%
# apply tSNE filter
tSNE_outliers = [
    'Neutrophils__P51',
    'Plasma_Norgen__P02',
    'Basophils__P37',
    'B_cells__P38',
    ]

MS_sc_aggreg_ad = MS_sc_aggreg_ad[~MS_sc_aggreg_ad.obs_names.isin(tSNE_outliers),:] 
MS_sc_aggreg_ad

#%%
# load Kiel blood catalogue data
Kiel_ad = ad.read_h5ad(Kiel_ad_path)
# remove WB and exosomes
Kiel_ad = Kiel_ad[~Kiel_ad.obs.group.isin(['WB', 'exosomes']),:]
# rename group 'serum' to 'Serum'
Kiel_ad.obs.group.cat.add_categories(['Serum'], inplace=True)
Kiel_ad.obs.loc[Kiel_ad.obs.group.isin(['serum']),'group'] = 'Serum'
Kiel_ad.obs.group.cat.remove_unused_categories(inplace=True)
# aggregate on subclass name
Kiel_sc_aggreg_ad = AnnDataAggregator(Kiel_ad, aggregate_obs=False, by='subclass_name', work_on_raw=True)
Kiel_sc_aggreg_ad

#%%
sc.tl.rank_genes_groups(MS_sc_aggreg_ad, 'mirsort_group', use_raw=False, method='wilcoxon', pts=True)
sc.pl.rank_genes_groups_dotplot(MS_sc_aggreg_ad, n_genes=5, standard_scale='var', size_title='Fraction of samples\nin group (%)', title='miR-Blood', save='__miRBlood__marker_sRNAs.pdf')

sc.tl.rank_genes_groups(Kiel_sc_aggreg_ad, 'group', use_raw=False, method='wilcoxon', pts=True)
sc.pl.rank_genes_groups_dotplot(Kiel_sc_aggreg_ad, n_genes=5, standard_scale='var', size_title='Fraction of samples\nin group (%)', title='IKMB catalogue', save='__Kiel__marker_sRNAs.pdf')


#%%
fig_3 = px.box(MS_sc_aggreg_ad.obs, x='mirsort_group', y='total_count_after_preprocessing', color='mirsort_group', color_discrete_map=celltype_lut, category_orders={'mirsort_group': ['Plasma', 'Erythrocytes', 'Thrombocytes', 'Monocytes', 'Neutrophils', 'Eosinophils', 'Basophils', 'NK cells', 'CD4+ T cells', 'CD8+ T cells', 'B cells']})
fig_3.update_layout(xaxis_title="", yaxis_title="total read counts after preprocessing", title='miR-Blood', showlegend=False, font_family='Arial', font_size=16)
fig_3.update_xaxes(tickangle=0)
fig_3.show()
fig_3.write_image(output_path + 'boxplot__totalcounts__miRBlood.pdf', width=1200, height=600)

fig_3b = px.box(Kiel_sc_aggreg_ad.obs, x='group', y='total_count', color='group', color_discrete_map=celltype_lut_b, category_orders={'group': ['Serum', 'CD235a', 'CD14', 'CD15', 'CD56', 'CD4', 'CD8', 'CD19']})
fig_3b.update_layout(xaxis_title="", yaxis_title="total read counts after preprocessing", title='IKMB catalogue (previous benchmark)', showlegend=False, font_family='Arial', font_size=16)
fig_3b.update_xaxes(tickangle=0)
fig_3b.show()
fig_3b.write_image(output_path + 'boxplot__totalcounts__Kiel.pdf', width=1200, height=600)


#%%
# check number of non-zero seqs per group
def get_non_zero_seqs(AnnData, groupby):
  counts_df = pd.DataFrame(data=AnnData.raw.X, index=AnnData.obs[groupby], columns=AnnData.var_names)
  counts_df = counts_df > 0
  counts_df = counts_df.groupby(counts_df.index).sum()
  return counts_df

def get_non_zero_seqs_fraction(AnnData, groupby):
  counts_df = pd.DataFrame(data=AnnData.raw.X, index=AnnData.obs[groupby], columns=AnnData.var_names)
  counts_df = counts_df > 0
  # calculate fraction
  counts_df = counts_df.groupby(counts_df.index).mean()
  return counts_df

MS_sc_aggreg_nz_counts_df = get_non_zero_seqs(MS_sc_aggreg_ad, 'mirsort_group')
Kiel_sc_aggreg_nz_counts_df = get_non_zero_seqs(Kiel_sc_aggreg_ad, 'group')

MS_sc_aggreg_nz_fract_df = get_non_zero_seqs_fraction(MS_sc_aggreg_ad, 'mirsort_group')
Kiel_sc_aggreg_nz_fract_df = get_non_zero_seqs_fraction(Kiel_sc_aggreg_ad, 'group')

#%%
if False:
    # per celltype (row) get list of features with non-zero expression in any sample
    MS_sc_aggreg_nz_counts_dict = MS_sc_aggreg_nz_counts_df.apply(lambda x: list(x[x>0].index), axis=1).to_dict()
    Kiel_sc_aggreg_nz_counts_dict = Kiel_sc_aggreg_nz_counts_df.apply(lambda x: list(x[x>0].index), axis=1).to_dict()
    # combine dicts, add suffix to keys
    combined_dict = {'miR-Blood -' + k: v for k, v in MS_sc_aggreg_nz_counts_dict.items()}
    combined_dict.update({'IKMB catalogue -' + k: v for k, v in Kiel_sc_aggreg_nz_counts_dict.items()})
    combined_dict.keys()

#%%
# per celltype (row) get list of features with non-zero fraction
MS_sc_aggreg_nz_fract_dict = MS_sc_aggreg_nz_fract_df.apply(lambda x: list(x[x==1].index), axis=1).to_dict()
Kiel_sc_aggreg_nz_fract_dict = Kiel_sc_aggreg_nz_fract_df.apply(lambda x: list(x[x==1].index), axis=1).to_dict()
# combine dicts, add suffix to keys
combined_fract_dict = {'miR-Blood - ' + k: v for k, v in MS_sc_aggreg_nz_fract_dict.items()}
combined_fract_dict.update({'IKMB catalogue - ' + k: v for k, v in Kiel_sc_aggreg_nz_fract_dict.items()})
combined_fract_dict.keys()

#%%
# generate upset plot for all 
if False:
    content_data = from_contents(combined_fract_dict)

    upset = UpSet(content_data, subset_size='auto', show_counts='%d')
    upset.plot()

    plt.savefig(output_path + 'upset_plot__nonzerofraction__miRBlood_vs_Kiel.png', dpi=300)

#%%
# plot per celltype
for celltype in group_name_dict.keys():
    
  tmp_dict = {k: v for k, v in combined_fract_dict.items() if k in [f"miR-Blood - {group_name_dict[celltype]['miR-Blood']}", f"IKMB catalogue - {group_name_dict[celltype]['Kiel']}"]}

  content_data = from_contents(tmp_dict)

  upset = UpSet(content_data, subset_size='auto', show_counts='%d')
  plt.rcParams['font.size'] = 15
  upset.plot()
  # increase width of plot slightly to show labels fully
  if celltype == 'Plasma':
    plt.gcf().set_size_inches(plt.gcf().get_size_inches()[0]+0.2, plt.gcf().get_size_inches()[1])
  plt.savefig(output_path + f'upset_plot__nonzerofraction__miRBlood_vs_Kiel__{celltype}.png', dpi=300)
  plt.clf()


#%%
# subset to subclass_name overlap
sc_overlapp = list(set.intersection(set(MS_sc_aggreg_ad.var_names),set(Kiel_sc_aggreg_ad.var_names)))
print(len(sc_overlapp))
MS_sc_aggreg_ad = MS_sc_aggreg_ad[:,sc_overlapp]
Kiel_sc_aggreg_ad = Kiel_sc_aggreg_ad[:,sc_overlapp]



#%%
# dictionary of group names
group_name_dict = {
  'Neutrophils': {'Kiel': 'CD15', 'miR-Blood': 'Neutrophils'},
  'Monocytes': {'Kiel': 'CD14', 'miR-Blood': 'Monocytes'},
  'NK cells': {'Kiel': 'CD56', 'miR-Blood': 'NK cells'},
  'CD4': {'Kiel': 'CD4', 'miR-Blood': 'CD4'},
  'CD8': {'Kiel': 'CD8', 'miR-Blood': 'CD8'},
  'B cells': {'Kiel': 'CD19', 'miR-Blood': 'B cells'},
  'Erythrocytes': {'Kiel': 'CD235a', 'miR-Blood': 'Erythrocytes'},
  'Plasma': {'Kiel': 'Serum', 'miR-Blood': 'Plasma'},
  #'EDTA_WholeBlood': {'Kiel': 'WB', 'miR-Blood': 'EDTA_WholeBlood'},
}

#%%
# compare shared expression of miR-Blood to Kiel (aggregated on subclass_name)
for celltype in group_name_dict.keys():
    print(celltype)
    x_data = MS_sc_aggreg_ad[MS_sc_aggreg_ad.obs.mirsort_group.str.match(group_name_dict[celltype]['miR-Blood']),:].X.mean(axis=0)
    y_data = Kiel_sc_aggreg_ad[Kiel_sc_aggreg_ad.obs.group.str.match(group_name_dict[celltype]['Kiel']),:].X.mean(axis=0)

    nas = np.logical_or(np.isnan(x_data), np.isnan(y_data))
    x_non_zero_shared = len(x_data[~np.isnan(x_data)])
    y_non_zero_shared = len(y_data[~np.isnan(y_data)])
    x_data = x_data[~nas]
    y_data = y_data[~nas]

    pearson_coef, pearson_pval = pearsonr(x_data, y_data)

    plt.rcParams['font.size'] = 15

    sns.scatterplot(x=x_data, y=y_data, s=10).set(title= f"{celltype} - Pearson coefficient = {str(round(pearson_coef, 2))}", xlabel=f"miR-Blood [log2(RPM+1)]", ylabel=f"IKMB catalogue [log2(RPM+1)]")
    #sns.scatterplot(x=x_data, y=y_data, s=10).set(title= f"{celltype} - Pearson coefficient = {str(round(pearson_coef, 2))}", xlabel=f"miR-Blood {group_name_dict[celltype]['miR-Blood']} expression [log2(RPM+1)]", ylabel=f"IKMB catalogue {group_name_dict[celltype]['Kiel']} expression [log2(RPM+1)]")
    # set y and x axis limits
    plt.xlim(0, 19)
    plt.ylim(0, 19)
    # same width and height
    plt.gca().set_aspect('equal', adjustable='box')
    # add dashed diagonal line
    plt.plot([0, 19], [0, 19], 'k--')
    # color miR-16-5p red
    plt.scatter(x_data[MS_sc_aggreg_ad.var_names[~nas].isin(['miR-16-5p'])], y_data[MS_sc_aggreg_ad.var_names[~nas].isin(['miR-16-5p'])], color='red')
    # make sure that x-axis title is displayed fully (currently cut off)
    plt.tight_layout()
    plt.savefig(output_path + 'scatterplot__sharedscs__miRBlood_vs_Kiel__' + celltype +'.png', facecolor='white', dpi=300)
    plt.clf()

#%%
# check highest expressed sRNAs for defined celltype
celltype = 'Neutrophils'
tmp_df = pd.DataFrame(data=MS_sc_aggreg_ad[MS_sc_aggreg_ad.obs.mirsort_group.str.match(group_name_dict[celltype]['miR-Blood']),:].X.mean(axis=0), index=MS_sc_aggreg_ad.var_names)
tmp_df.sort_values(0,ascending=False)

#%%
celltype = 'Plasma'
tmp_df = pd.DataFrame(data=Kiel_sc_aggreg_ad[Kiel_sc_aggreg_ad.obs.group.str.match(group_name_dict[celltype]['Kiel']),:].X.mean(axis=0), index=Kiel_sc_aggreg_ad.var_names)
tmp_df.sort_values(0,ascending=False)












#%%