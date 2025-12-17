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

# Check file exist
if (length(args)==0) {message("Please enter an valid file name.")} else
{  
  marker <- tolower(args[1]) #Must be character 
  path_decipher_trained <- (args[2]) #Must be character 
}

# Debug
# marker <- "teleo"
# path_decipher_trained <- "utils/teleo_trained.rds"

# Open seqtab - abundance table output of dada2
seqtab <- readRDS(paste0("outputs/02_dada2/", marker, "/seqtab_clean.rds"))

dna <- DNAStringSet(getSequences(as.matrix(seqtab))) # Create a DNAStringSet from the ASVs
dna@ranges@NAMES <- colnames(seqtab)
trainingset <- readRDS(paste0("utils/", marker, "_trained.rds")) 
ids <- IdTaxa(dna, trainingset, strand="both", processors=8, verbose=TRUE, threshold = 60) # use all processors
ids_df <- decipher_list_to_df(ids, confidence = TRUE)
colnames(ids_df)[1] <- "sequence"
 
write.csv(ids_df, paste0("outputs/02_dada2/", marker, "/taxo_assigned_decipher.csv"), row.names= FALSE)

# Now link it to abundance table
# Open data 
samples <- read.csv(paste0("outputs/02_dada2/", marker, "/seqtab_clean_t.csv"))
colnames(samples)[1] <- "sequence"

# Now link 
decipher_taxo_complete <- dplyr::left_join(ids_df, samples) 

# Now write the output
write.csv(decipher_taxo_complete, paste0("outputs/", marker, "_decipher.csv"), row.names=TRUE)