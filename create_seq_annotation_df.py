import json
import pandas as pd

# load config
config = json.load(
    open("MIRSORT_ANNOTATION_DF.json")
)

inputpath_1 = config["generate_sRNA_sub_class_annotation_df.py"]["sRNA_sub_class_annot_path"]
outputpath = config["create_seq_annotation_df.py"]["ANNOTATION_DF"]


# read small RNA (sub)class annotation dataframe
annotation_df = pd.read_csv(inputpath_1, index_col=0)

# add columns with basic sequence info: 'length', 'g_fraction', 'a_fraction', 't_fraction', 'c_fraction', 'gc_fraction'
annotation_df['length'] = annotation_df.index.str.len()
annotation_df['g_fraction'] = annotation_df.index.str.count('G')/annotation_df['length']
annotation_df['a_fraction'] = annotation_df.index.str.count('A')/annotation_df['length']
annotation_df['t_fraction'] = annotation_df.index.str.count('T')/annotation_df['length']
annotation_df['c_fraction'] = annotation_df.index.str.count('C')/annotation_df['length']
annotation_df['gc_fraction'] = annotation_df.index.str.count('G|C')/annotation_df['length']

# store number of columns and their names after reordering
col_num_pre = len(annotation_df.columns)
col_names_pre = set(annotation_df.columns) 

# reorder columns
annotation_df = annotation_df[['length',
 'g_fraction',
 'a_fraction',
 't_fraction',
 'c_fraction',
 'gc_fraction',
 'Annotation_unitas',
 'Annotation_sports',
 'small_RNA_class_annotation',
 'subclass_name']]

# check if all columns are present after reordering
col_num_after = len(annotation_df.columns)
col_names_after = set(annotation_df.columns) 
if col_num_pre != col_num_after:
    print(str(col_num_pre-col_num_after) + ' column(s) were lost during reordering:')
    print(col_names_pre - col_names_after)


# write dataframe to be used to annotate var of AnnData objects
annotation_df.to_csv(outputpath)
