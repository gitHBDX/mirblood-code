# Pre-processing by MuK
## TODO: add code

# Install conda environment
conda env create --file annot_env.yml

# Activate conda environment
conda activate annot_env

# Run sRNA annotation pipelines (unitas or sports)
## Use fasta file of sequences after pre-processing for sequence annotation
## Allow 1 missmatch in annotation pipelines (unitas and sports)
## Generate fasta file from snoDB tsv file (downloaded from https://bioinfo-scottgroup.med.usherbrooke.ca/snoDB/)
python snoDB2fa.py
## Run unitas (https://www.smallrnagroup.uni-mainz.de/software.html)
perl unitas_1.7.7.pl -i features_detected_sequences__publication.fa -species homo_sapiens -species_miR_only -tail 2 -intmod 1 -mismatch 1 -insdel 0 -refseq snoDB.fa -dump_prefix unitas_annotation/UNITAS
## Run sports (https://github.com/junchaoshi/sports1.1)
### NOTE: get sports pre-compiled 'Homo_sapiens' annotation database from https://ncrnainfo-my.sharepoint.com/personal/sports_ncrna_info/_layouts/15/guestaccess.aspx?docid=0773ed3d5f6b74f35bbd643e1af221c31&authkey=AcRxf8walnGUIEhgI--8CDc
perl sports.pl -i features_detected_sequences__publication.fa -p 4 -k -M 1 -g Homo_sapiens/genome/hg38/genome -m Homo_sapiens/miRBase/21/miRBase_21-hsa -r Homo_sapiens/rRNAdb/human_rRNA -t Homo_sapiens/GtRNAdb/hg19/hg19-tRNAs -w Homo_sapiens/piRBase/piR_human -e Homo_sapiens/Ensembl/release-89/Homo_sapiens.GRCh38.ncrna -f Homo_sapiens/Rfam/12.3/Rfam-12.3-human -o sports_annotation/

# Drop all sequences that do not have any annotation in unitas or sports
cd unitas_annotation/UNITAS_dd-mm-yyyy_features_detected_sequences__publication.fa_#1
awk 'NF>=3' unitas.full_annotation_matrix.txt | awk '$3 !~ "low_complexity" {print $0}' > unitas.full_annotation_matrix_justannoseqs.txt
cd sports_annotation/1_features_detected_sequences__publication/features_detected_sequences__publication_result
awk '$6 !~ "NO_Annotation" {print $0}' features_detected_sequences__publication_output.txt > features_detected_sequences__publication_output_justannoseqs.txt

# Merge annotations
python merge_sRNAclass_annotations.py

# Get subclassification for rRNAs and YRNAs
python rRNA_position_classification.py

# Generate sRNA subclass annotation 
python generate_sRNA_sub_class_annotation_df.py

# Create sequence-based sRNA annotation dataframe
python create_seq_annotation_df.py

# Reduce AnnData to features with 'subclass_name' annotation
python ad_reduce_features.py

# Aggregate AnnData based on 'subclass_name'
python ad_aggregate.py

# Create subclass-name-based annotation dataframe and add to aggregated AnnData
python create_aggregated_annotation_df.py 

# Reduce features by expression threshold and subset to blood components and whole blood samples
python ad_reduce_features_further.py

# Create csv files for dashboard
python ad2csv.py






