rm(list = ls())

######################## PACKAGES
library(argparse)
library(RColorBrewer)
#library(WriteXLS)

######################## HELPERS
source("src/helpers/helpers_0_general.R")
source("src/helpers/helpers_NGS.R")

######################## ARGUMENTs

parser = ArgumentParser()
parser$add_argument("--input-path", help="Path of the input fasta/fastq folder", required=T)
parser$add_argument("--output-path", help="Path of the output folder", required=T)
parser$add_argument("--count-type", help="READ or UMI", required=T, default = "READ")
parser$add_argument("--total-count", help="TRUE for total-count (READ or UMI), FALSE for number of unique sequences (READ or UMI is ignored)", required=T, default = T)


args = parser$parse_args(
  # list(
  #  "--input-path=../data/preprocessed/trimmed_and_collapsed/",
  #  "--output-path=../data/preprocessed/stats/read_length_distribution__UMI/",
  #  "--count-type=NA",
  #  "--total-count=False"
  # )
)


# READ or UMI
count_type = args$count_type

# TRUE for total-count (READ or UMI), FALSE for number of unique sequences
is_total_count = args$total_count

# input folder: processed files
input_path = paste0(args$input_path, "/")

# output folder: read statistics
output_path = paste0(args$output_path, "/")



######################## METHODS

# get the read_length distribution from a single file
get_read_length_distribution = function(filename){
  all_lines = readLines(filename)
  sequences = NA
  file_type = get_NGS_file_type(all_lines)
  
  if(file_type == "fastq"){
    headers = all_lines[seq(1, length(all_lines), by = 4)]
    sequences = all_lines[seq(2, length(all_lines), by = 4)]
  } else if (file_type == "fasta") {
    headers = all_lines[seq(1, length(all_lines), by = 2)]
    sequences = all_lines[seq(2, length(all_lines), by = 2)]
  } else {
    if(!(getFileExtension(filename) %in% c("fa", "fasta"))){
      stop("filetype is wrong")
    } else {
      headers = c()
      sequences = c()
    }
  }
  
  counts = unlist(lapply(headers, get_count_from_header))
  if(is_total_count == F) counts = rep(1, length(sequences))
  
  read_length = unlist(lapply(sequences, FUN = function(x) nchar(x)))
  read_length_dist = rep(0, max(50, max(read_length)))
  
  for(i in 1:length(read_length_dist)){
    tmp_len = i
    tmp_pos = which(read_length == tmp_len)
    if(length(tmp_pos)>0){
      read_length_dist[tmp_len] = read_length_dist[tmp_len] + sum(counts[tmp_pos])
    }
  }
  
  return(read_length_dist)
}


# plot a read length distribution
plot_read_length_distribution = function(filename, smp, len_dist, total_count = NA, count_type = "READ"){
  
  x = 1:length(len_dist)
  y = len_dist
  names(len_dist) = c(1:length(len_dist))
  max_additional_counts = if(max(len_dist) < 10^6) 10^5 else 10^6
  
  ylab = paste(count_type, "count")
  subtitle = paste0("Total ", count_type ," count: ", total_count)
  if(is_total_count==F){
    ylab = "Number of unique sequences"
    subtitle = paste0("Total number of unique sequences: ", total_count)
  }
  
  pdf(filename, width = 9)
  bp = barplot(len_dist, ylim = c(0,max(len_dist)+max_additional_counts), xlim = c(0,length(len_dist)), width = rep(0.8, length(len_dist)), cex.names = 0.4, main = smp, xlab = "Read length", ylab = ylab)
  if(!is.na(total_count)) mtext(subtitle, line = 0)
  dev.off()
}

######################## MAIN

msg = NA

if(is_total_count %in% c(F, FALSE, "False", "FALSE", "false", "F")) is_total_count = F


if(is_total_count){
  msg = paste(count_type, "counts")
} else {
  msg = "number of unique sequences"
}
cat("Calculating read statistics ...", msg, "\n")

# get all preprocessed files
files = list.files(input_path, full.names=T)
smps = c()



# create the vectors for counts and read_lengths distribution, as separated files, as overall and as table with smps as separated rows
counts_all_smps = c()
read_length_dist_overall = c()
read_length_dist_all_smps_as_separated_rows = c()

for(i in 1:length(files)){
  tmp_filename_input = files[i]
  tmp_smp = getSampleName(tmp_filename_input)
  smps = c(smps, tmp_smp)
  tmp_filename_output = paste0(output_path, tmp_smp, ".pdf")
  cat(i, "of", length(files), "-- read stats -- ", tmp_smp, "\n")
  
  ##### calculate read length distribution
  tmp_read_len = get_read_length_distribution(filename = tmp_filename_input)
  tmp_total_count = get_total_count(filename = tmp_filename_input, is_total_count = is_total_count)
  counts_all_smps = c(counts_all_smps, tmp_total_count)
  if(i == 1){
    read_length_dist_overall = tmp_read_len
    read_length_dist_all_smps_as_separated_rows = tmp_read_len
  } else {
    read_length_dist_overall = read_length_dist_overall + tmp_read_len
    read_length_dist_all_smps_as_separated_rows = rbind(read_length_dist_all_smps_as_separated_rows, tmp_read_len)
  }
  # plot read length distribution
  plot_read_length_distribution(filename = tmp_filename_output, smp = tmp_smp, len_dist = tmp_read_len, total_count = tmp_total_count, count_type = count_type)
}



# plot overall read-length-distribution
ofile = paste0(output_path, "0_all_samples_overall.pdf")
plot_read_length_distribution(filename = ofile, smp = "Overall (all samples)", len_dist = read_length_dist_overall, total_count = sum(counts_all_smps), count_type = count_type)

# create read-length-distribution as table
read_length_dist_all_smps_as_separated_rows = data.frame(read_length_dist_all_smps_as_separated_rows, stringsAsFactors = F)
colnames(read_length_dist_all_smps_as_separated_rows) = paste0("len_", c(1:ncol(read_length_dist_all_smps_as_separated_rows)))
read_length_dist_all_smps_as_separated_rows = cbind(counts_all_smps, read_length_dist_all_smps_as_separated_rows)

# write out the tables
read_length_dist_all_smps_as_separated_rows = cbind(smpID = smps, read_length_dist_all_smps_as_separated_rows)
ofile = paste0(output_path, "0_all_samples.txt")
write.table(read_length_dist_all_smps_as_separated_rows, ofile, row.names = F, sep = "\t", quote = F)
ofile = paste0(output_path, "0_all_samples.xlsx")
res_xlsx = read_length_dist_all_smps_as_separated_rows
#WriteXLS(res_xlsx, ofile)




