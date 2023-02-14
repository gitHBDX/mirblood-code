rm(list = ls())

######################## packages
library(argparse)

######################## HELPERS

source("src/helpers/helpers_0_general.R")
source("src/helpers/helpers_NGS.R")

######################## ARGUMENTS

parser = ArgumentParser()
parser$add_argument("--output-path", help="Path of the output folder", required=T)
parser$add_argument("--input-path-total-read-count-before-preprocessing", help="Path of the input folder; total read count, before preprocessing", required=T)
parser$add_argument("--input-file-total-count-after-preprocessing", help="Path of the input file; total  count, after preprocessing", required=T)
parser$add_argument("--input-path-flowcell-stats", help="Path of the input folder; flowcell", required=T)
parser$add_argument("--input-path-index-stats", help="Path of the input folder; index", required=T)
parser$add_argument("--input-path-sequencer-stats", help="Path of the input folder; sequencer", required=T)



args = parser$parse_args(
  # list("--output-path=../../data/preprocessed/0_info/",
  #      "--input-path-total-read-count-before-preprocessing=../data/preprocessed/stats/fastq_basic_infos__total_read_count_before_preprocessing/",
  #      "--input-file-total-count-after-preprocessing=../data/preprocessed/stats/read_length_distribution__UMI/0_all_samples.txt",
  #      "--input-path-flowcell-stats=../data/preprocessed/stats/fastq_basic_infos__flowcell_count/",
  #      "--input-path-index-stats=../data/preprocessed/stats/fastq_basic_infos__index_count/",
  #      "--input-path-sequencer-stats=../data/preprocessed/stats/fastq_basic_infos__sequencer_count/"
  #      )
)
output_path = paste0(args$output_path, "/")

input_path_total_read_count_before_preprocessing = paste0(args$input_path_total_read_count_before_preprocessing, "/")
input_file_total_count_after_preprocessing = paste0(args$input_file_total_count_after_preprocessing)
input_path_flowcell_stats = paste0(args$input_path_flowcell_stats, "/")
input_path_index_stats = paste0(args$input_path_index_stats, "/")
input_path_sequencer_stats = paste0(args$input_path_sequencer_stats, "/")


######################## METHODS

get_info = function(input_file){
  input = read.delim(input_file, header = F, stringsAsFactors = F)
  res_vals = unique(input$V2)
  res_vals = res_vals[order(res_vals)]
  res_string = paste(res_vals, collapse = ",")
  return(res_string)
}

get_count = function(input_file){
  input = read.delim(input_file, header = F, stringsAsFactors = F)
  
  res_string = paste(unique(input$V1), collapse = ",")
  return(res_string)
}
######################## 

cat("Creating sample info ...", "\n")


files = list.files(input_path_flowcell_stats, full.names = F, pattern = ".txt")
smpIDs = unlist(lapply(files,getSampleName))
hbdxIDs = unlist(lapply(smpIDs, get_smpID_from_NGS_file))


total_read_counts_before_preprocessing = c()
flowcells = c()
indices = c()
sequencer = c()

cat("Extracting flowcell counts", "\n")
for(i in 1:length(files)){
  tmp_file = files[i]
  cat(i, "of", length(files), " -- extracting infos -- ", tmp_file, "\n")
  
  tmp_file_full = paste0(input_path_total_read_count_before_preprocessing, "/", tmp_file)
  tmp_val = get_count(tmp_file_full)
  total_read_counts_before_preprocessing = c(total_read_counts_before_preprocessing, tmp_val)
  
  tmp_file_full = paste0(input_path_flowcell_stats, "/", tmp_file)
  tmp_val = get_info(tmp_file_full)
  flowcells = c(flowcells, tmp_val)

  tmp_file_full = paste0(input_path_index_stats, "/", tmp_file)
  tmp_val = get_info(tmp_file_full)
  indices = c(indices, tmp_val)
  
  tmp_file_full = paste0(input_path_sequencer_stats, "/", tmp_file)
  tmp_val = get_info(tmp_file_full)
  tmp_val = gsub("@", "", tmp_val)
  sequencer = c(sequencer, tmp_val)
}


total_count_info = read.delim(input_file_total_count_after_preprocessing, stringsAsFactors = F, check.names = F, sep = "\t")
total_count_after_preprocessing = total_count_info$counts_all_smps[match(smpIDs, total_count_info$smpID)]
total_count_ratio__after_to_before = as.numeric(total_count_after_preprocessing) / as.numeric(total_read_counts_before_preprocessing)

res = data.frame(smpID = smpIDs, HBDx_ID = hbdxIDs, flowcells, sequencer, total_read_counts_before_preprocessing, total_count_after_preprocessing, total_count_ratio__after_to_before,stringsAsFactors = F, check.names = F)




output = paste0(output_path, "/sample_info.txt")
write.table(res, output, row.names = F, sep = "\t", quote = F)
