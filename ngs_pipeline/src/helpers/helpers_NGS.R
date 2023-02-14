

# return the file type: fastq or fasta
get_NGS_file_type = function(lines){
  file_type = "UNKNOWN"
  if(substr(lines[1],1,1) %in% c("@", ">")){
    if(substr(lines[3],1,1) == "+"){
      file_type = "fastq"
    } else {
      file_type = "fasta"
    }
  }
  return(file_type)
}

# get the single count from the header
get_count_from_header = function(header){
  return(as.numeric(unlist(strsplit(header, "_x"))[2]))
}

# get the total read/umi/number-of-unique-sequences count from a file
get_total_count = function(filename, is_total_count = T){
  seq_ids = NA
  
  all_lines = readLines(filename)
  file_type = get_NGS_file_type(all_lines)
  
  
  if(file_type == "fastq"){
    headers = all_lines[seq(1, length(all_lines), by = 4)]
    seq_ids = headers
  } else if (file_type == "fasta") {
    headers = all_lines[seq(1, length(all_lines), by = 2)]
    seq_ids = headers
  } else {
    if(!(getFileExtension(filename) %in% c("fa", "fasta"))){
      stop("filetype is wrong")
    } else {
     return(0)
    }
  }
  
  counts = unlist(lapply(headers, get_count_from_header))
  
  res = sum(counts)
  
  if(is_total_count == F){
    # number of unique sequences
    res = length(headers)
  }
  
  return(res)
}



create_fastAQ_table = function(input_file){
  
  cat(input_file, "\n")
  all_lines = readLines(input_file)
  file_type = get_NGS_file_type(all_lines)
  
  fastaq_lines = NA
  if(file_type == "fastq"){
    fastaq_lines = 4
  } else if(file_type == "fasta"){
    fastaq_lines = 2
  } else {
    if(!(getFileExtension(input_file) %in% c("fa", "fasta"))){
      msg = paste("no fasta or fastq file:", input_file)
      stop(msg)
    } else {
      result = data.frame(matrix(NA, nrow = 0, ncol = 4))
      return(result)
    }
  }
  
  
  headers = all_lines[seq(1, length(all_lines), fastaq_lines)]
  sequences = all_lines[seq(2, length(all_lines), fastaq_lines)]
  sequence_length = unlist(lapply(sequences, nchar))
  counts = unlist(lapply(headers, FUN = function(x) unlist(strsplit(x, "_x"))[2]))
  
  result = data.frame(headers,
                      sequences,
                      sequence_length,
                      counts,
                      stringsAsFactors = F)
  return(result)
}



write_fastaq_file = function(filename, fastaq_df){
  
  fastaq_lines = NA
  if(grepl(".fastq", filename)){
    fastaq_lines = 4
  } else if(grepl(".fa|.fasta", filename)){
    fastaq_lines = 2
  } else {
    msg = paste("no fasta or fastq file:", filename)
    stop(msg)
  }
  
  all_lines = c()
  
  tmp_df = fastaq_df[,c(1:fastaq_lines)]
  
  if(fastaq_lines == 4){
    all_lines = c(rbind(tmp_df$headers, tmp_df$sequences, tmp_df$comments, tmp_df$qualities))
  } else {
    all_lines = c(rbind(tmp_df$headers, tmp_df$sequences))
  }
  
  write(all_lines, filename, append = F)
}



norm_RPM = function(df, is_log2_transformation = TRUE, total_counts_input = NA){
  
  colsms = NA
  if(length(total_counts_input)==1){
    if(is.na(total_counts_input)){
      colsms = colSums(df)
    }
  } else {
    if(!is.null(names(total_counts_input)) & all(colnames(df) %in% names(total_counts_input))){
      colsms = as.numeric(total_counts_input[colnames(df)])
    } else {
      msg = "Providing wrong sample names for colsums in normaliziation via RPM."
      stop(msg)
    }
  }
  
  for(i in 1:ncol(df)){
    
    df[,i] = df[,i]/colsms[i]
    df[,i] = df[,i] * 10^6
    #df[,i] = df[,i] + (1 - min(df[,i])) ### to check
    
    if(is_log2_transformation){
      df[,i] = df[,i] + 1
      df[,i] = log2(df[,i])
    }
  }
  
  return(df)
}

norm_RPM_without_log2_and_shift = function(df, total_counts_input = NA){
  
  colsms = NA
  if(length(total_counts_input)==1){
    if(is.na(total_counts_input)){
      colsms = colSums(df)
    }
  } else {
    if(!is.null(names(total_counts_input)) & all(colnames(df) %in% names(total_counts_input))){
      colsms = as.numeric(total_counts_input[colnames(df)])
    } else {
      msg = "Providing wrong sample names for colsums in normaliziation via RPM."
      stop(msg)
    }
  }
  
  for(i in 1:ncol(df)){
    
    df[,i] = df[,i]/colsms[i]
    df[,i] = df[,i] * 10^6
    #df[,i] = df[,i] + (1 - min(df[,i])) ### to check
    #df[,i] = df[,i] + 1
    #df[,i] = log2(df[,i])
  }
  
  return(df)
}
                         

get_flowcell_ID_from_folder_name = function(input){
  tmp = input
  res = rev(unlist(strsplit(tmp, "_")))[1]
  return(res)
}


## geth hbdx id of different filenames created by the lab, EMBL or other
method_list_for_deriving_final_filename = c("hbdxID", "hbdxID_extended", "full_name")
get_smpID_from_NGS_file = function(input, correct_with_underscores=T, method_for_deriving_final_filename = method_list_for_deriving_final_filename[1]){
  tmp_input = input
  tmp_input = basename(tmp_input)
  
  if (!(method_for_deriving_final_filename %in% method_list_for_deriving_final_filename)){
    msg = paste("Following method is not defined for deriving filenames/smpIDs:", method_for_deriving_final_filename)
    stop(msg)
  }
  
  if(method_for_deriving_final_filename == method_list_for_deriving_final_filename[1]){
    if(grepl("lane1", tmp_input)){
      # EMBL_Data
      # example: HC2CGBGXH_HBDx_463_Pool_1_20s004658-1-1_Rasti_lane14630019_sequence.txt
      tmp_input = unlist(strsplit(tmp_input, "lane1"))[2]
      tmp_input = gsub("_sequence.txt.gz", "", tmp_input)
    } else if(grepl("R1_001", tmp_input)){
      #
      # example: 454_A_0027_C_S68_R1_001.fastq.gz, 473_0025_S12_L001_R1_001.fastq.gz
      tmp_input = unlist(strsplit(tmp_input,"_S[0-9]{1,}_" ))[1]
      
      
    } else if(grepl("__", tmp_input)){
      #
      # example: 468_0064__AAAGJ55M5__S1_L001_R1_001.fastq.gz
      # when "hbdxID" is chosen, you get the HBDx_ID, in case of "full_name" there will be no modification in this part
      tmp_input = unlist(strsplit(tmp_input,"__" ))[1]
    }
    
    if(correct_with_underscores){
      if(!grepl("_", tmp_input)){
        
        tmp_input_split = unlist(strsplit(tmp_input, ""))
        
        
        ID = ""
        for(i in 1:length(tmp_input_split)){
          if(i == 4){
            ID = paste0(ID, "_")
          }
          
          if(grepl("[[:alpha:]]", tmp_input_split[i])){
            if(i > 7){
              ID = paste0(ID, "_")
            }
          }
          
          ID = paste0(ID, tmp_input_split[i])
          
          if(grepl("[[:alpha:]]", tmp_input_split[i])){
            if(i == 4){
              ID = paste0(ID, "_")
            }
          }
          
        }
        tmp_input = ID
      }
    }
  }
  
  if(method_for_deriving_final_filename == method_list_for_deriving_final_filename[3]){
    
    # example: 454_A_0027_C_S68_R1_001.fastq.gz, 473_0025_S12_L001_R1_001.fastq.gz
    tmp_input = unlist(strsplit(tmp_input,"_S[0-9]{1,}_" ))[1]
  }
  
   tmp_input = gsub(".fastq.gz|.fastq|.fasta.gz|.fasta|.fa.gz|.fa", "", tmp_input)
    
  return(tmp_input)
}
