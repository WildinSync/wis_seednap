# Link taxo with files 
library(tidyverse)

# Get the arguments
args <- commandArgs(T)

# Check file exist
if (length(args)==0) {message("Please enter an valid file name.")} else
{  
  marker <- tolower(args[1]) #Must be character 
}
# Debug
# marker <- "teleo"

# Open data 
samples <- read.csv(paste0("outputs/02_dada2/", marker, "/seqtab_clean_t.csv"))
colnames(samples)[1] <- "sequence"

# Taxa with ecotag
ecotag_taxo <- read.csv(paste0("outputs/03_taxo/", marker, "/query_ecotag.tsv"), sep = "\t")
ecotag_taxo$sequence <- toupper(ecotag_taxo$sequence)
ecotag_taxo_complete <- dplyr::left_join(ecotag_taxo, samples) 

# Export 
write.csv(ecotag_taxo_complete, paste0("outputs/", marker, "_ecotag.csv"), row.names= FALSE)