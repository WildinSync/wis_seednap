#!/usr/bin/env Rscript
# generate cutadapt.R 
# usage : Rscript generate_cutadapt_ligation.R lab_tags_files.csv

# ------------------------- # 
# FUNCTIONS

# Reverse complement
rc <- function(z){
  rc1 <- function(zz){
    s <- strsplit(zz, split = "")[[1]]
    s <- rev(s)
    dchars <- strsplit("ACGTMRWSYKVHDBNI", split = "")[[1]]
    comps <- strsplit("TGCAKYWSRMBDHVNI", split = "")[[1]]
    s <- s[s %in% dchars] # remove spaces etc
    s <- dchars[match(s, comps)]
    s <- paste0(s, collapse = "")
    return(s)
  }
  z <- toupper(z)
  tmpnames <- names(z)
  res <- unname(sapply(z, rc1))
  if(!is.null(attr(z, "quality"))){
    strev <- function(x) sapply(lapply(lapply(unname(x), charToRaw), rev), rawToChar)
    attr(res, "quality") <- unname(sapply(attr(z, "quality"), strev))
  }
  names(res) <- tmpnames
  return(res)
}

# convert df to fasta file
df_to_fasta <- function(file, output_file_path){
  fa <- character(2 * nrow(file))
  fa[c(TRUE, FALSE)] = sprintf("> %s", file[,1])
  fa[c(FALSE, TRUE)] = as.character(file[,2])
  writeLines(fa, output_file_path)
}

# ENF OF FUNCTIONS 
# ------------------------- # %   

# Library
library(tidyverse)
suppressWarnings(suppressMessages(library("tidyverse", quietly = TRUE)))

args<-commandArgs(T)

metadata_lab <- args[1] # Must be a csv file - the corr tags
# directory_output <- args[2] # Must a directory name

dir.create("outputs/00_demultiplex_ligation/cutadapt_tags/", recursive = TRUE)
directory_output <- "outputs/00_demultiplex_ligation/cutadapt_tags/"

#file_metadata <- "raw/metadata/Corr_tags_Teleo.csv"
#directory_output <- "raw/metadata/cutadapt_tags/Teleo/"

file_metadata_df <- read.csv(metadata_lab, sep=",", h=T, stringsAsFactors=F)
file_metadata_df <- file_metadata_df[,1:3] # sometimes extra column - messes with things
colnames(file_metadata_df) <- c("eventID", "tag_demultiplex", "library")

# Split by library 
list_file_metadata_df <- split(file_metadata_df, file_metadata_df$library)

# Generate a file per library 
files_to_write <- lapply(list_file_metadata_df, function(x){
  file_output <- x %>% 
    dplyr::mutate(cutadapt = paste0(tag_demultiplex, ";min_overlap=8...", tolower(rc(tag_demultiplex)), ";min_overlap=8")) %>% 
    dplyr::select(eventID, cutadapt, library)
})

# export these files 
lapply(files_to_write, function(x){
  name_output <- paste0(directory_output, "/", unique(x$library), ".fasta")
  df_to_fasta(x, name_output)
})
