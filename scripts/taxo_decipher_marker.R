# Taxo using dada2

# Libs
library(dada2); packageVersion("dada2")
library(Biostrings)
library(dplyr)
library(DECIPHER)

## FUNCTIONS

# Decipher output unpractical file format, change it back to dataframe
decipher_list_to_df <- function(file, confidence = FALSE){
  
  require(tidyverse)
  file_df <- bind_rows(lapply(file, function(x){dplyr::as_data_frame(x)}), .id = "seq_uniq") %>% 
    select(-confidence) %>% 
    pivot_wider(names_from = rank, values_from = taxon, values_fn = list) %>% 
    as.data.frame() %>% 
    # remove column lists 
    transmute_all(~sapply(., toString, sep=";"))
  
  if(confidence == TRUE){
    
    file_df_confidence <- bind_rows(lapply(file, function(x){dplyr::as_data_frame(x)}), .id = "seq_uniq") %>% 
      select(-taxon) %>% 
      pivot_wider(names_from = rank, values_from = confidence, names_prefix = "confidence_", values_fn = list) %>% 
      as.data.frame() %>% 
      # remove column lists 
      transmute_all(~sapply(., toString, sep=";"))# %>% 
      #mutate(seq_uniq = as.numeric(seq_uniq))%>%
      #arrange(seq_uniq)
    
    file_df <- merge(file_df, file_df_confidence) 
  }
  return(file_df)
}
## END OF FUNCTIONS

# Get the arguments
args <- commandArgs(T)

if (length(args) < 2) {
  stop("Please provide marker and trained classifier path arguments.")
}

marker <- tolower(args[1]) # Must be character
path_decipher_trained <- args[2] # Must be character
threshold <- if (length(args) >= 3) as.integer(args[3]) else 60
processors <- if (length(args) >= 4) as.integer(args[4]) else 8
output_dir <- if (length(args) >= 5) args[5] else "outputs"

# Debug
# marker <- "teleo"
# path_decipher_trained <- "utils/teleo_trained.rds"

# Open seqtab - abundance table output of dada2
seqtab <- readRDS(file.path(output_dir, "02_dada2", marker, "seqtab_clean.rds"))

dna <- DNAStringSet(getSequences(as.matrix(seqtab))) # Create a DNAStringSet from the ASVs
dna@ranges@NAMES <- colnames(seqtab)
trainingset <- readRDS(path_decipher_trained)
ids <- IdTaxa(dna, trainingset, strand="both", processors=processors, verbose=TRUE, threshold=threshold) # use configured processors/threshold
ids_df <- decipher_list_to_df(ids, confidence = TRUE)
colnames(ids_df)[1] <- "sequence"
 
write.csv(ids_df, file.path(output_dir, "02_dada2", marker, "taxo_assigned_decipher.csv"), row.names= FALSE)

# Now link it to abundance table
# Open data 
samples <- read.csv(file.path(output_dir, "02_dada2", marker, "seqtab_clean_t.csv"))
colnames(samples)[1] <- "sequence"

# Now link 
decipher_taxo_complete <- dplyr::left_join(ids_df, samples) 

# Now write the output
write.csv(decipher_taxo_complete, file.path(output_dir, paste0(marker, "_decipher.csv")), row.names=TRUE)
