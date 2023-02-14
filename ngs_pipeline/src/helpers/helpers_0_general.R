

# get operating system
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
  return(os)
}


# get file extension
getFileExtension <- function(file, with_dot = F){ 
  ex <- rev(unlist(strsplit(basename(file), split="\\.")))[1]
  if(with_dot) ex = paste0(".", ex)
  return(ex)
}


# get only sample name without file extension
getSampleName = function(file){
  file_extension = getFileExtension(file = file, with_dot = T)
  smp_name = gsub(file_extension, "", basename(file))
  smp_name = gsub(".fastq|.fasta", "", smp_name)
  return(smp_name)
}



