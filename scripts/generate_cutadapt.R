#!/usr/bin/env Rscript
# generate cutadapt.R 
# usage : Rscript generate_cutadapt corr.tags.files output_directory

source("scripts/functions.R")

# Library
library(tidyverse)
suppressWarnings(suppressMessages(library("tidyverse", quietly = TRUE)))

args<-commandArgs(T)

file_metadata <- args[1] # Must be a csv file - the corr tags
directory_output <- args[2] # Must a directory name

#file_metadata <- "raw/metadata/Corr_tags_Teleo.csv"
#directory_output <- "raw/metadata/cutadapt_tags/Teleo/"

file_metadata_df <- read.csv(file_metadata, sep=",", h=T, stringsAsFactors=F)
file_metadata_df <- file_metadata_df[,1:6] # sometimes extra column - messes with things
colnames(file_metadata_df) <- c("plaque", "run", "sample_name", "tag", "primerF", "primerR")

# Split by library 
list_file_metadata_df <- split(file_metadata_df, file_metadata_df$run)

# Generate a file per library 
files_to_write <- lapply(list_file_metadata_df, function(x){
  file_output <- x %>% 
    dplyr::mutate(cutadapt = paste0(tag, ";min_overlap=8...", tolower(rc(tag)), ";min_overlap=8")) %>% 
    dplyr::select(sample_name, cutadapt, run)
})

# export these files 
lapply(files_to_write, function(x){
  name_output <- paste0(directory_output, "/", unique(x$run), ".fasta")
  df_to_fasta(x, name_output)
})