rm(list = ls())

######################## packages
library(argparse)

######################## HELPERS

source("src/helpers/helpers_0_general.R")


######################## ARGUMENTS

parser = ArgumentParser()
parser$add_argument("--input", help="input file", required=T)
parser$add_argument("--output", help="output file", required=T)
parser$add_argument("--adapter-3p", help="3p adapter - The full adapter is required in the read", required=T)
parser$add_argument("--insert-length-min", help="Length of the biological insert", required=T, type="integer")
parser$add_argument("--use-umi-counts", help="Use TRUE for umi counts, FALSE for read counts.", required=T)
parser$add_argument("--umi-length", help="The expected UMI-length that should be used.", required=F, type ="integer")
parser$add_argument("--read-2-adapter", help="Read 2 adapter", required=F)
parser$add_argument("--read-2-adapter-min-overlap", help="Minimum overlap to read 2 adapter", required=F, type="integer")


args = parser$parse_args(
  # list(
  #   "--input=../../results_test/raw_ngs_preprocessing/pipeline/input/test_file_2.txt",
  #   "--output=../../results_test/raw_ngs_preprocessing/pipeline/output/results_UMI/trimmed_and_collapsed/test_file__debug_2.txt",
  #   "--adapter-3p=AAAAAAAAAA",
  #   "--insert-length-min=18",
  #   "--use-umi-counts=TRUE",
  #   "--umi-length=10",
  #   "--read-2-adapter=ACAGTCGAT",
  #   "--read-2-adapter-min-overlap=4"
  # )
)

input_file = args$input
output_file = args$output
adapter_3p = args$adapter_3p
insert_length_min = args$insert_length_min
umi_mode = args$use_umi_counts

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

# determine in which lines are the sequences (each line, every 2nd line for fa/fasta or every 4th line for fastq)
print_only_sequences = ""
tmp_cmd = paste0(cmd_read_in, " | head -n 4")
first_4_lines = system(tmp_cmd, intern = T)
if(substr(first_4_lines[1], 1, 1) == "@"){
  if(substr(first_4_lines[3], 1, 1) == "+"){
    print_only_sequences = " | sed -n '2~4p'" # to print each 4th line (fastq-format) starting from second line which is a sequence
  } else {
    print_only_sequences = " | sed -n '2~2p'" # to print each 2nd line (fasta-format) starting from second line which is a sequence
  }
  
  if(get_os() == "osx") print_only_sequences = gsub("sed", "gsed", print_only_sequences)
  
} else {
  print_only_sequences = "" # to print each line (assuming all lines are sequences)
}

# start with command for read-in
cmd = paste0(cmd_read_in, print_only_sequences)
count_tag_for_the_header = "SEQ_read_count" # this tag for the header; it is dependent on the UMI-mode

# derive the command for UMI-count or read-count
if(umi_mode){
  ### UMI-count
  umi_length = args$umi_length
  read_2_adapter = args$read_2_adapter
  read_2_adapter_min_overlap = args$read_2_adapter_min_overlap
  read_2_adapter_partial = substr(read_2_adapter, 1, min(read_2_adapter_min_overlap, nchar(read_2_adapter)))
  count_tag_for_the_header = "UMI_count"
  
  cmd = paste0(cmd,
               " | grep -E '(", adapter_3p, ")[ACGTN]{", umi_length,"}(", read_2_adapter_partial, ")'", # get all reads with insert-min-length and the 3p-adapter, UMI with fixed length and min-overlap of the Read-2-adapter
               " | awk -F'", adapter_3p, "' '{u=substr($2,1,", umi_length,"); print $1,u}'", # split the read between the 3p-adapter and return the subsequence/insert on the left-side of the 3p-adapter and the UMI on the right-side: [$1][adapter_3p][substr($2,1,umi_length)]
               " | sort | uniq -c | awk '{print $2}'" # collapse to unique insert-UMI-combinations (1st "sort", then 2nd "uniq -c") so that each insert-UMI-combinations occur exactly once and return only the insert ($2): [$1 counts][$2 insert][$3 umi]
  )
} else {
  ### read-count
  cmd = paste0(cmd,
               " | grep -E '(", adapter_3p, ")'", # get all reads with insert-min-length and the 3p-adapter
               " | awk -F'", adapter_3p, "' '{print $1}'" # split the read between the 3p-adapter and return the subsequence/insert on the left-side of the 3p-adapter: [$1][adapter_3p][$2]
  )
}

# complete the final command
cmd = paste0(cmd,
             " | grep -v 'N'", # discard reads with N
             " | awk '{if(length($0)>=", insert_length_min, "){print $0}}'", # get only sequences that are at least long as the insert_length_min
             " | sort | uniq -c | sort -k 1nr -k 2", # sort the sequences via command "sort" to collapse and count them via command "uniq -c" (sorting is required before collapsing); then sort 1st by counts decreasingly and 2nd by sequences in lexicographical ascending way
             " | awk '{print \">",count_tag_for_the_header, "__x\"$1\"\\n\"$2}'", # save the read/umi count ($1) in the header following in the next line by the sequence ($2))
             " > ", output_file # write the data to the output-file
)


cat(cmd, "\n")
res = system(cmd, wait = T, intern = T) # execute the command
cat(res, "\n")
