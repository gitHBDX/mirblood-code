import json
import os
import pandas as pd
from io import StringIO

# load config
config = json.load(
    open("MIRSORT_ANNOTATION_DF.json")
)

inputpath = config["rRNA_position_classification.py"]["seq_path"]
inputpath_2 = config["rRNA_position_classification.py"]["ref_path"]
outputpath = config["rRNA_position_classification.py"]["rRNA_name_annot_path"]
outputpath_2 = config["rRNA_position_classification.py"]["YRNA_name_annot_path"]


# function to generate rRNA/YRNA name based on binned start position
def sRNA_name(sRNA_type_list, col_name):
    annot_df = pd.DataFrame(columns=[col_name])

    for rRNA_type in sRNA_type_list:

        seq_path = inputpath + 'features_detected_sequences__publication_match_rRNA_' + rRNA_type + '_match_genome.fa' 
        ref_path =  inputpath_2 + 'human_rRNA_' + rRNA_type + '.fa'


        # use seqkit locate to identify position within parental rRNA (like in SPORTS allow for 2 mismatches and only positive strand)
        stream = os.popen('seqkit locate -m 2 -P -f ' + seq_path + ' ' + ref_path)
        output = stream.read()
        locate_df = pd.read_csv(StringIO(output), sep='\t', header=0, index_col=3)


        # divide maximal existing end position by 25 (= typical small RNA size) to set bins
        bins = int(locate_df.end.max()/25)
        print(rRNA_type + ' bins: ' + str(bins))

        # generate annotation based on binned start position
        locate_df[col_name] = pd.cut(locate_df['start'], bins, precision=0, labels=['rRNA-' + rRNA_type + '_bin-' + str(i+1) for i in range(bins)])
        print(pd.cut(locate_df['start'], bins, precision=0).sort_values().unique())
        locate_df = locate_df[[col_name]]
        annot_df = annot_df.append(locate_df)
    
    # remove duplicates (these are very repetitive sequences)
    annot_df = annot_df.loc[~annot_df.index.duplicated(),:]

    return annot_df

rRNA_annot_df = sRNA_name(['5S','5.8S','18S','28S','45S','12S','16S'], 'rRNA_name')
YRNA_annot_df = sRNA_name(['RNY1','RNY3','RNY4','RNY5'], 'YRNA_name')


# write to file
rRNA_annot_df.to_csv(outputpath)
YRNA_annot_df.to_csv(outputpath_2)


