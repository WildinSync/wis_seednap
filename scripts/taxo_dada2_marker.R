# Taxo using dada2

# Libs
suppressMessages(suppressWarnings({
  library(dada2)
  library(Biostrings)
  library(dplyr)
}))

# Get the arguments
args <- commandArgs(T)

# Check file exist
if (length(args)==0) {message("Please enter an valid file name.")} else
{  
  marker <- tolower(args[1]) #Must be character 
  path_dada_all <- (args[2]) #Must be character 
  path_dada_species <- (args[3]) #Must be character 
}

# Debug
# marker <- "teleo"
# path_dada_all <- "utils/teleo_clean_crabs_taxo.fasta"
# path_dada_species <- "utils/teleo_clean_crabs_species.fasta"

cat(paste0("marker is ", marker, " \n the path dada all is ", path_dada_all, " \n and the path dada species is ", path_dada_species))

# Open seqtab - abundance table output of dada2
seqtab <- readRDS(paste0("outputs/02_dada2/", marker, "/seqtab_clean.rds"))

# Assign taxo dada2
taxa <- assignTaxonomy(seqtab, path_dada_all, multithread = TRUE, tryRC = TRUE)
# Sometimes issue that Species column exist at this stage - remove it
if ("Species" %in% colnames(taxa)) {
  taxa <- taxa[, !(colnames(taxa) == "Species")]
}
taxa_plus <- addSpecies(taxa, path_dada_species, verbose = TRUE, allowMultiple = TRUE)
colnames(taxa_plus) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
taxa_plus <- as.data.frame(taxa_plus)
taxa_plus$sequence <- rownames(taxa_plus)

write.csv(taxa_plus, paste0("outputs/02_dada2/", marker, "/", "taxonomy_dada2RDP.csv"), row.names=TRUE)

# Now link it to abundance table
# Open data 
samples <- read.csv(paste0("outputs/02_dada2/", marker, "/seqtab_clean_t.csv"))
colnames(samples)[1] <- "sequence"

# Now link 
dada2_RDP_taxo_complete <- dplyr::left_join(taxa_plus, samples) 

# Now write the output
write.csv(dada2_RDP_taxo_complete, 
          paste0("outputs/", marker, "_dada2RDP.csv"), 
          row.names = FALSE)