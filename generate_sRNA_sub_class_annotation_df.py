import json
import numpy as np
import pandas as pd

# load config
config = json.load(
    open("MIRSORT_ANNOTATION_DF.json")
)

inputpath = config["merge_sRNAclass_annotations.py"]["sRNA_class_annot_path"]
inputpath_4 = config["generate_sRNA_sub_class_annotation_df.py"]["snoDB_info_path"]
inputpath_5 = config["rRNA_position_classification.py"]["rRNA_name_annot_path"]
inputpath_6 = config["rRNA_position_classification.py"]["YRNA_name_annot_path"]
outputpath = config["generate_sRNA_sub_class_annotation_df.py"]["sRNA_sub_class_annot_path"]

# load merged annotations
annotation_merge = pd.read_csv(inputpath, index_col=0)


# generate annotation matrix for small RNA classes 
# based on the annotations from unitas (snoDB included) and sports
# sequential mode: miRNA > tRNA > rRNA > YRNA > snoRNA > lncRNA > snRNA > piRNA
annotation_merge["small_RNA_class_annotation"] = np.where(
    (
        annotation_merge.Annotation_unitas.str.contains("mir-", case=False, na=False)
        | annotation_merge.Annotation_unitas.str.contains("let-", na=False)
    ), 'miRNA', 
    np.where(
    (
        annotation_merge.Annotation_unitas.str.contains("tRF", na=False)
        | annotation_merge.Annotation_unitas.str.contains("tR-half", na=False)
    ),'tRNA',
    np.where(
        annotation_merge.Annotation_sports.str.contains("S-rRNA", na=False), 
    'rRNA',
    np.where(
    annotation_merge.Annotation_sports.str.contains("YRNA", na=False), 
    'YRNA',
    np.where(
    (
        annotation_merge.Annotation_unitas.str.contains("snoRNA", na=False)
        | annotation_merge.Annotation_unitas.str.contains("snoDB", na=False)
    ),'snoRNA',
    np.where(
    (
        annotation_merge.Annotation_unitas.str.contains("lncRNA", na=False)
    ),
    'lncRNA',
    np.where(
    (
        annotation_merge.Annotation_unitas.str.contains("snRNA", na=False)
    ),'snRNA',
    np.where(
    (
        annotation_merge.Annotation_unitas.str.contains("piR", na=False)
    ),
    'piRNA',
    'no_annotation',
    ))))))))
# NOTE: ignore vault_RNA, lincRNA, scaRNA, sRNA, ribozyme, protein_coding and DNA labels



# add miRNA_names
## concat miR names of unitas and MuK
miR_unitas = (
    annotation_merge[annotation_merge.small_RNA_class_annotation == 'miRNA']
    .Annotation_unitas.str.split(";", expand=True)
    .apply(lambda x: x.str.split("(", 1, expand=True)[0])
)
miR_unitas = miR_unitas.apply(lambda x: x.str.cat(sep=";"), axis=1)
## remove redundant precursor names
def redundantPrecursorRemover(string_list):
    out = []
    for s in string_list:
        if not any([s.lower() in r.lower() for r in string_list if s != r]):
            out.append(s)
    return out
miR_unitas = miR_unitas.apply(lambda x: ";".join(redundantPrecursorRemover(list(set(x.split(";"))))))
## add miR name information
annotation_merge.loc[annotation_merge.small_RNA_class_annotation == 'miRNA', "subclass_name"] = miR_unitas
miR_unitas


# add tRNA_names
tRNA_unitas = (
    annotation_merge[annotation_merge.small_RNA_class_annotation == 'tRNA'].Annotation_unitas.str.replace(' from ','__', regex=False).str.replace(r' \([0-9]*-[0-9]*\)','', regex=True).str.split(";", expand=True).apply(lambda x: x.str.replace(r'.*\|.*', '', regex=True)).apply(lambda x: x.str.split("-ENSG|>ENST", 1, expand=True)[0]).apply(lambda x: x.str.replace(r'-[0-9]*-[0-9]$', '', regex=True)).apply(lambda x: ";".join(redundantPrecursorRemover(list(set(x.dropna())))),axis=1)
)
## add tRNA/MT name information
annotation_merge.loc[annotation_merge.small_RNA_class_annotation == 'tRNA', "subclass_name"] = tRNA_unitas
tRNA_unitas


# add snoRNA_names
sno_unitas = annotation_merge[annotation_merge.small_RNA_class_annotation == 'snoRNA'].Annotation_unitas.str.split(";", expand=True).apply(lambda x: x.str.split("snoRNA\|", 1, expand=True)[1]).apply(lambda x: x.str.replace(r'.*\|.*', '', regex=True)).apply(lambda x: x.str.split("-ENSG|>ENST", 1, expand=True)[0]).apply(lambda x: x.str.cat(sep=";"), axis=1).str.replace(r'_\[[0-9]*\]','', regex=True)
## use only SNODB names from snoDB (rest from standard unitas annot)
snoDB_unitas = annotation_merge[annotation_merge.small_RNA_class_annotation == 'snoRNA'].Annotation_unitas.str.split(";", expand=True).apply(lambda x: x.str.split("snoDB\|", 1, expand=True)[1]).apply(lambda x: x.str.replace(r'.*\|.*', '', regex=True)).apply(lambda x: x.str.split(" ", 1, expand=True)[0]).apply(lambda x: x.str.cat(sep=";"), axis=1).str.replace(r'_\[[0-9]*\]','', regex=True)
sno_concat = snoDB_unitas + ';' + sno_unitas
sno_concat = sno_concat.apply(lambda x: ";".join(redundantPrecursorRemover(list(set(x.split(";"))))))
annotation_merge.loc[annotation_merge.small_RNA_class_annotation == 'snoRNA', "subclass_name"] = sno_concat
sno_concat


# function to directly add name of unitas for certain small RNA class to annotation df
def add_unitas_name(sRNA_bool, unitas_flag):
    unitas_name = annotation_merge[annotation_merge.small_RNA_class_annotation == sRNA_bool].Annotation_unitas.str.split(";", expand=True).apply(lambda x: x.str.split(unitas_flag, 1, expand=True)[1]).apply(lambda x: x.str.replace(r'.*\|.*', '', regex=True)).apply(lambda x: x.str.split("(", 1, expand=True)[0]).apply(lambda x: x.str.replace('>','', regex=False).str.replace('-$','',regex=True).str.rstrip()).apply(lambda x: x.str.cat(sep=";"), axis=1)
    unitas_name = unitas_name.apply(lambda x: ";".join(redundantPrecursorRemover(list(set(x.split(";"))))))
    annotation_merge.loc[annotation_merge.small_RNA_class_annotation == sRNA_bool, "subclass_name"] = unitas_name


# add lncRNA_names
add_unitas_name('lncRNA', 'lncRNA\|')

# add snRNA_names
add_unitas_name('snRNA', 'snRNA\|')

# add piRNA_names
add_unitas_name('piRNA', 'piR-cluster:')



# remove -ENST/G affix (esp. impacts lncRNAs and snRNAs)
annotation_merge.subclass_name = annotation_merge.subclass_name.str.replace(r'-ENS(G|T)[0-9]*\.[0-9]*','', regex=True)


# remove ENST annotations
annotation_merge.subclass_name = annotation_merge.subclass_name.str.replace(r'ENST[0-9]*\.[0-9]*;','',regex=True).str.replace(r'ENST[0-9]*\.[0-9]*','',regex=True)


# set sequence as index
annotation_merge = annotation_merge.set_index("Sequence")


# read and add rRNA name based on binned start position within parental rRNA
rRNA_name_df = pd.read_csv(inputpath_5, index_col=0, names=['subclass_name'], skiprows=1)
rRNA_overlap = list(set.intersection(set(annotation_merge[annotation_merge.small_RNA_class_annotation=='rRNA'].index),set(rRNA_name_df.index)))
rRNA_name_df = rRNA_name_df.loc[rRNA_overlap,:]
rRNA_name_df.reindex(annotation_merge.loc[rRNA_name_df.index,:].index)
annotation_merge.loc[rRNA_name_df.index,'subclass_name'] = rRNA_name_df.subclass_name
annotation_merge[annotation_merge.small_RNA_class_annotation=='rRNA']

# read and add YRNA name based on binned start position within parental YRNA
YRNA_name_df = pd.read_csv(inputpath_6, index_col=0, names=['subclass_name'], skiprows=1)
YRNA_overlap = list(set.intersection(set(annotation_merge[annotation_merge.small_RNA_class_annotation=='YRNA'].index),set(YRNA_name_df.index)))
YRNA_name_df = YRNA_name_df.loc[YRNA_overlap,:]
YRNA_name_df.reindex(annotation_merge.loc[YRNA_name_df.index,:].index)
annotation_merge.loc[YRNA_name_df.index,'subclass_name'] = YRNA_name_df.subclass_name
annotation_merge[annotation_merge.small_RNA_class_annotation=='YRNA']



# remove seqs without subclass_name
annotation_merge = annotation_merge[~annotation_merge.subclass_name.str.contains(r'^\s*$', na=True)]
annotation_merge


# remove duplicated subclass_names
annotation_merge.subclass_name = annotation_merge.subclass_name.apply(lambda x: ";".join(redundantPrecursorRemover(list(set(x.split(";"))))))


# remove multimaps
annotation_merge = annotation_merge[~annotation_merge.subclass_name.str.contains(';')] 


print(annotation_merge.small_RNA_class_annotation.value_counts())


print("number of unique subclass names: " + str(len(annotation_merge.subclass_name.unique())))


# check number of unique subclass_names per small_RNA_class
print("number of unique subclass names per small RNA class: ") 
print(annotation_merge.groupby('small_RNA_class_annotation').subclass_name.unique().apply(lambda x: len(x)))


# check number of subclass names that occure more than 10 times
print("number of subclass names that occure more than 10 times: ") 
print(pd.Series(annotation_merge.subclass_name.value_counts().values > 10).value_counts())


# write to file
annotation_merge.to_csv(outputpath)


