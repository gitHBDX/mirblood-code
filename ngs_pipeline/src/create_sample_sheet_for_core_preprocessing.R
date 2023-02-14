rm(list = ls())

source("src/helpers/helpers_0_general.R")
source("src/helpers/helpers_NGS.R")
##################

if(T){
  input_folder = "<PATH>"
  output_folder = "<PATH>"
  dir.create(output_folder, recursive = T)
  
  
  files = list.files(input_folder, pattern = ".fastq.gz", full.names = T, recursive = T)
  files = files[!grepl("Undetermined", files)]
  
  input_files = c(files)
  output_sample_names = gsub(".fastq.gz", "", basename(input_files))

  res_df = data.frame(input_files, output_sample_names, stringsAsFactors = F, check.names = F)
  ofile = paste0(output_folder, "/sample_sheet_for_core_preprocessing.txt")
  write.table(res_df, ofile, row.names = F, quote = F, sep = "\t")
}

