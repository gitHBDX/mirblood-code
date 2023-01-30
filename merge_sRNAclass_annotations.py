import json
import pandas as pd

# load config
config = json.load(
    open("MIRSORT_ANNOTATION_DF.json")
)

inputpath1 = config["merge_sRNAclass_annotations.py"]["unitas_annot_path"]
inputpath2 = config["merge_sRNAclass_annotations.py"]["sports_annot_path"]
outputpath = config["merge_sRNAclass_annotations.py"]["sRNA_class_annot_path"]


# load unitas annotation matrix
unitas_annot = pd.read_csv(
    inputpath1,
    sep="\t",
    names=[
        "Reads",
        "Annotation",
        "d",
        "e",
        "f",
        "g",
        "h",
        "i",
        "j",
        "k",
        "l",
        "m",
        "n",
        "o",
        "p",
        "q"
    ],
    skiprows=1,
)

# combine annotation in one column
unitas_annot = unitas_annot.reset_index()
unitas_annot = unitas_annot.rename(columns={"index": "Sequence"})
annot_concat = unitas_annot.iloc[:, 2:17].apply(lambda x: x.str.cat(sep=";"), axis=1)
unitas_annot = unitas_annot.iloc[:, 0:1]
unitas_annot["Annotation"] = annot_concat


# load sports annotation matrix
sports_annot = pd.read_csv(inputpath2, sep="\t")
# bring to same format like unitas
sports_annot = sports_annot[["Sequence", "Annotation"]]

# merge unitas and sports annotation
seq_overlap = set.intersection(set(unitas_annot.Sequence), set(sports_annot.Sequence))
annotation_merge = unitas_annot.merge(
    sports_annot,
    how="left",
    left_on="Sequence",
    right_on="Sequence",
    suffixes=("_unitas", "_sports"),
)


# write annotation_merge dataframe to be used to generate annotation matrix for small RNA classes and subclasses
annotation_merge.to_csv(outputpath)

