rm(list = ls())

######################## packages
library(argparse)

######################## HELPERS

source("src/helpers/helpers_0_general.R")


######################## ARGUMENTS

parser = ArgumentParser()
parser$add_argument("--input", help="input file", required=T)
parser$add_argument("--input-smpID", help = "input smpID", required=T)
parser$add_argument("--output-path-fastq-stats", help="output path", required=T)
parser$add_argument("--output-path-index-stats", help="output path", required=T)
parser$add_argument("--output-path-total-read-count-stats", help="output path", required=T)
parser$add_argument("--output-path-sequencer-stats", help="output path", required=T)



args = parser$parse_args(
  # list(
  #   "--input=/data/hbdx/raw_data/ngs_data/test/fastq/30178-003.R1.fastq.gz",
  #   "--input-smpID=30178-003",
  #   "--output-path-fastq-stats=results/stats/raw_fastq_flowcell_count/",
  #   "--output-path-index-stats=results/stats/raw_fastq_index_count/"
  # )
)

input_file = args$input
output_file_basename = paste0(args$input_smpID, ".txt")
output_path_fastq_stats = args$output_path_fastq_stats
output_path_index_stats = args$output_path_index_stats
output_path_total_read_count_stats = args$output_path_total_read_count_stats
output_path_sequencer_stats = args$output_path_sequencer_stats

######################## 

# determine which pager-tool should be used: zcat (compressed) or cat (non-compressed)
pager_tool = "cat"
if(getFileExtension(input_file) %in% c("gz")) pager_tool = "zcat"

# determine the operating system to use zcat properly
# linux: zcat input_file
# mac-os: zcat < input_file
piping_sign_for_the_operating_system = ""
if(get_os() == "osx") {
  if(pager_tool == "zcat") piping_sign_for_the_operating_system = "<"
}

# basic read in command
cmd_read_in = paste(pager_tool, piping_sign_for_the_operating_system, input_file)

# determine in which lines are the sequences (each line, every first line for fa/fasta or every 4th line for fastq)
print_only_headers = ""
tmp_cmd = paste0(cmd_read_in, " | head -n 4")
first_4_lines = system(tmp_cmd, intern = T)
if(substr(first_4_lines[1], 1, 1) == "@"){
  if(substr(first_4_lines[3], 1, 1) == "+"){
    print_only_headers = " | sed -n '1~4p'" # to print each 4th line (fastq-format) starting from first line which is a header
  } else {
    print_only_headers = " | sed -n '1~2p'" # to print each 2nd line (fasta-format) starting from first line which is a header
  }
  
  if(get_os() == "osx") print_only_headers = gsub("sed", "gsed", print_only_headers)
  
} else {
  print_only_headers = "" # to print each line (assuming all lines are headers)
}

# start with command for read-in
cmd_begin = paste0(cmd_read_in, print_only_headers)

# flowcell count
cmd_flowcell = " | awk -F':' '{print $3}'" # split the header with colon; third string is the flowcell

# index count
cmd_index = " | awk -F':' '{print $10}'" # split the header with colon; 10th string is the index

# sequencer count
cmd_sequencer = " | awk -F':' '{print $1}'" # split the header with colon; first string is the sequencer

# cmd total read count without preprocessing
cmd_total_read_count = " | wc -l"


cmd_list = c(cmd_flowcell, cmd_index, cmd_total_read_count, cmd_sequencer)
output_path_list = c(output_path_fastq_stats, output_path_index_stats, output_path_total_read_count_stats, output_path_sequencer_stats)

# cmd general part
cmd_general_part = paste(" | sort | uniq -c | sort -k 1nr -k 2", # sort the sequences via command "sort" to collapse and count them via command "uniq -c" (sorting is required before collapsing); then sort 1st by counts decreasingly and 2nd by sequences in lexicographical ascending way
                         " | awk '{print $1\"\\t\"$2}'" # save the read/umi count ($1) in the header following in the next line by the sequence ($2))
                         )


for(i in 1:length(cmd_list)){
  # derive the command for getting infos: flowcell
  output_file = paste0(output_path_list[i], "/", output_file_basename)
  
  if(cmd_list[i] %in% c(cmd_total_read_count)){
    cmd = paste(cmd_begin, cmd_list[i])
  } else {
    cmd = paste(cmd_begin,
                cmd_list[i],
                cmd_general_part
                )
  }
  cmd = paste(cmd, " > ", output_file) # write the data to the output file
  
  cat(cmd, "\n")
  res = system(cmd, wait = T, intern = T) # execute the command
}
