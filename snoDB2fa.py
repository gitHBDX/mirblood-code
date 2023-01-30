#%%
import pandas as pd
import json

#%%
# load config
config = json.load(
    open("MIRSORT_ANNOTATION_DF.json")
)
input_path = config["generate_sRNA_sub_class_annotation_df.py"]["snoDB_info_path"]
outputpath = config["merge_sRNAclass_annotations.py"]["snoDB_annot_path"]

#%%
# load snoDB tsv file
snoDB_df = pd.read_csv(input_path, sep='\t', index_col=0)
snoDB_df

#%%
# generate fasta-file
ofile = open(outputpath, "w")

for i in range(len(snoDB_df)):
    ofile.write(">snoDB|" + snoDB_df.iloc[i,:].Symbol + "\n" + snoDB_df.iloc[i,:].Sequence + "\n")

ofile.close()

#%%