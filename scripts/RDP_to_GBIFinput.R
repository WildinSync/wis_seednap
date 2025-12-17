#From dada2 output to ready for GBIF input 
suppressPackageStartupMessages({
  library(tidyverse)
})

# Get the arguments
args <- commandArgs(T)

# Check file exist
if (length(args) == 0) {message("Please enter an valid file name.")} else
{  
  file_path <- args[1] #Must be character 
}

# Debug
# file_path <- "outputs/mamm07_dada2RDP.csv"

source("scripts/functions.R")

# Start processing
cat(paste0("Processing to convert dada2 output to GBIF input compatible for file ", file_path, "\n"))

# Open file
file <- read.csv(file_path) 
col_keep <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "sequence")

# Clean X column if present
if ("X" %in% colnames(file)) {
  file <- file %>% select(-X)
}

# Format properly - long format, remove 0, add rank and taxon, order columns properly
file_out <- file %>%
  # reformat
  pivot_longer(!all_of(col_keep), names_to = "filter_code", values_to = "nb_reads") %>%
  filter(nb_reads > 0) %>%
  # Add infos
  add_rank_dada() %>%
  add_taxon_dada() %>%
  # Final ordering
  dplyr::select(kingdom, phylum, class,order, family, genus, species, taxon, rank, sequence, nb_reads, filter_code)%>%
  rename(eventID = filter_code) 

# Export file
prefix <- sub("\\.csv$", "", file_path)
write.csv(file_out, paste0(prefix, "_gbif_input.csv"), row.names=FALSE)

# Message
cat(paste0("Wrote processed file to  ", paste0(prefix, "_gbif_input.csv \n")))