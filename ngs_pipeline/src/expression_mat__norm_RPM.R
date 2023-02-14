rm(list = ls())

######################## PACKAGES
library(argparse)


######################## HELPERS

source("src/helpers/helpers_NGS.R")



######################## ARGUMENTS

parser = ArgumentParser()
parser$add_argument("--input-file-mat", help="Expression matrix", required=T)
parser$add_argument("--output-file-mat", help="Expression matrix (norm RPM and log2)", required=T)
parser$add_argument("--input-file-total-count", help="File with total counts", required=F, default = "NA")

args = parser$parse_args(
  # list(
  #   "--input-file-mat=../data/preprocessed/expression_mat/raw_count__mapped_and_collapsed_mir.txt",
  #   "--output-file-mat=../data/preprocessed/expression_mat/norm_RPM__mapped_and_collapsed_mir.txt"
  #   "--input-file-total-count=./data/preprocessed/stats/read_length_distribution__UMI/0_all_samples.txt"
  # )
)

input_file = args$input_file_mat
output_file = args$output_file_mat
total_count_file = args$input_file_total_count
######################## MAIN

cat("Normalization: Calculating log2(RPM)-values ...", "\n")

input_path = dirname(input_file)
#output_file = paste0(input_path, "/", "norm_RPM_log2_", gsub("raw_", "", basename(input_file)))


raw_mat_file = input_file
raw_mat = read.delim(raw_mat_file, stringsAsFactors = F, check.names = F)
if(colnames(raw_mat)[1] == "") raw_mat = read.delim(raw_mat_file, stringsAsFactors = F, check.names = F, row.names = 1)

norm_mat = NA
norm_mat_non_log = NA
if(total_count_file == "NA"){
  norm_mat = norm_RPM(raw_mat)
  norm_mat_non_log = norm_RPM(raw_mat, is_log2_transformation = FALSE)
} else {
  total_count_df = read.delim(total_count_file, stringsAsFactors = F, check.names = F)
  total_count_vals = total_count_df$counts_all_smps
  names(total_count_vals) = total_count_df$smpID
  
  norm_mat = norm_RPM(raw_mat, is_log2_transformation = TRUE, total_counts_input = total_count_vals)
  
  norm_mat_non_log = norm_RPM(raw_mat, is_log2_transformation = FALSE, total_counts_input = total_count_vals)
}


cat("Output:", output_file, "\n")
write.table(norm_mat, output_file, row.names = T, sep = "\t", quote = F)

output_file_non_log = gsub("_log2", "", output_file)
write.table(norm_mat_non_log, output_file_non_log, row.names = T, sep = "\t", quote = F)
