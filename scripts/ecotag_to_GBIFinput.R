#From dada2 output with ecotag to ready for GBIF input 
suppressPackageStartupMessages({
  library(tidyverse)
})

source("scripts/functions.R")

# Get the arguments
args <- commandArgs(T)

# Check file exist
if (length(args) == 0) {message("Please enter an valid file name.")} else
{  
  file_path <- args[1] #Must be character 
}

# Debug
# file_path <- "outputs/teleo_ecotag.csv"

### FUNCTIONS 
# Function to rename specific columns
rename_columns <- function(df) {
  df %>%
    rename_with(~ str_replace_all(., c("family_name" = "family", 
                                       "species_name" = "species", 
                                       "genus_name" = "genus",
                                       "order_name" = "order")))
}
#### END OF FUNCTIONS 

# Start processing
cat(paste0("Processing to convert dada2 output from ecotag to GBIF input compatible for file ", file_path, "\n"))

# Open file
file <- read.csv(file_path) 
file$kingdom <- "NA"
file$phylum <- "NA"
file$class <- "NA"
col_keep <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "sequence", "rank")

# Clean X column if present
if ("X" %in% colnames(file)) {
  file <- file %>% select(-X)
}

# Format properly - long format, remove 0, add rank and taxon, order columns properly
file_out <- file %>%
  # Clean columns
  dplyr::select(-id, -definition, -count, -family, -genus, -order, -scientific_name, -species) %>% 
  dplyr::select(-matches("^best_identity|^best_match|^match_count|^species_list|^taxid")) %>%
  rename_columns() %>%
  # reformat
  pivot_longer(!all_of(col_keep), names_to = "filter_code", values_to = "nb_reads") %>%
  filter(nb_reads > 0) %>%
  # Add infos
  add_rank_dada() %>%
  add_taxon_dada() %>%
  # Final ordering
  dplyr::select(kingdom, phylum, class,order, family, genus, species, taxon, rank, sequence, nb_reads, filter_code)

# Export file
prefix <- sub("\\.csv$", "", file_path)
write.csv(file_out, paste0(prefix, "_gbif_input.csv"), row.names=FALSE)

# Message
cat(paste0("Wrote processed file to  ", paste0(prefix, "_gbif_input.csv \n")))
