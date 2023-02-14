rm(list = ls())

######################## PACKAGES
library(argparse)

######################## HELPERS
source("src/helpers/helpers_0_general.R")
source("src/helpers/helpers_NGS.R")

######################## FUNCTIONS

evaluate_threshold_feature_present_in_min_X_samples = function(){
  res = min(feature_present_in_min_X_samples, length(files))
  if(res <= 0) res = 1
  return(res)
}

######################## ARGUMENTS

parser = ArgumentParser()
parser$add_argument("--input-path", help="Path of the input fasta folder", required=T)
parser$add_argument("--output-path", help="Path of the output fasta folder", required=T)
parser$add_argument("--feature-present-in-min-X-samples", help="In how many samples a feautre should be present", required=T, type="integer", default = 1)
parser$add_argument("--feature-with-at-least-X-count", help="Sequences with at least X counts are considered", required=T, type="integer", default = 1)

args = parser$parse_args(
  # list("--input-path=../data/preprocessed/trimmed_and_collapsed/",
  #      "--output-path=../data/preprocessed/expression_mat/",
  #      "--feature-present-in-min-X-samples=3",
  #      "--feature-with-at-least-X-count=5"
  #      )
)

input_path = paste0(args$input_path, "/")
output_path = paste0(args$output_path, "/")
feature_present_in_min_X_samples = args$feature_present_in_min_X_samples
feature_with_at_least_X_count = args$feature_with_at_least_X_count

######################## MAIN


cat("Creating raw count matrix (features as sequences) ...", "\n")


files = list.files(input_path, full.names = T)
all_features = c()

file_extension = getFileExtension(files[1], with_dot = T)


cat("Extracting only feature IDs", "\n")
for(i in 1:length(files)){
  tmp_file = files[i]
  cat(i, "of", length(files), " -- extracting feature IDs -- ", basename(tmp_file), "\n")
  
  tmp_fastaq = create_fastAQ_table(tmp_file)
  if(nrow(tmp_fastaq) == 0) next
 #tmp_fastaq = tmp_fastaq[1:10,]
  all_features = c(all_features, tmp_fastaq$sequences)
}

table_features = table(all_features)
feature_present_in_min_X_samples = evaluate_threshold_feature_present_in_min_X_samples()
table_features_filtered = table_features[table_features >= feature_present_in_min_X_samples]
all_features = names(table_features_filtered)

all_features = unique(all_features)
all_features = all_features[order(all_features)]



expr_mat = data.frame(matrix(0, nrow = length(all_features), ncol = length(files)), stringsAsFactors = F)
rownames(expr_mat) = all_features
colnames(expr_mat) = unlist(lapply(files, getSampleName))


cat("Creating expression matrix", "\n")
for(i in 1:length(files)){
  tmp_file = files[i]
  cat(i, "of", length(files), " -- creating expression matrix -- ", basename(tmp_file), "\n")
  
  tmp_basename = getSampleName(tmp_file)
  tmp_fastaq = create_fastAQ_table(tmp_file)
  if(nrow(tmp_fastaq) == 0) next
  
  tmp_fastaq = tmp_fastaq[tmp_fastaq$sequences %in% all_features,]
  tmp_counts = unlist(lapply(tmp_fastaq$headers, get_count_from_header))
  
  if(feature_with_at_least_X_count > 1){
    tmp_pos = which(tmp_counts >= feature_with_at_least_X_count)
    tmp_fastaq = tmp_fastaq[tmp_pos,]
    tmp_counts = tmp_counts[tmp_pos]
  }
  
  expr_mat[tmp_fastaq$sequences, tmp_basename] = tmp_counts
}

if(feature_with_at_least_X_count > 1){
  tmp_row_sums = rowSums(expr_mat)
  tmp_row_sums_with_at_least_X_count = tmp_row_sums[tmp_row_sums >= feature_with_at_least_X_count]
  tmp_features_with_at_least_X_count = names(tmp_row_sums_with_at_least_X_count)
  expr_mat = expr_mat[tmp_features_with_at_least_X_count,]
  all_features = tmp_features_with_at_least_X_count
}

output = paste0(output_path, "/features_all_sequences.txt")
write(all_features, output)

all_features_header = paste0(">seq__", all_features)
all_features_fasta = c(rbind(all_features_header, all_features))
output = paste0(output_path, "/features_all_sequences.fa")
write(all_features_fasta, output)

output = paste0(output_path, "/raw_count__sequences.txt")
write.table(expr_mat, output, row.names = T, sep = "\t", quote = F)
