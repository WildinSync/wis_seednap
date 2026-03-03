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
if (length(args) < 3) {
  stop("Please provide marker, RDP db path, and species db path arguments.")
}

marker <- tolower(args[1])
path_dada_all <- args[2]
path_dada_species <- args[3]
output_dir <- if (length(args) >= 4) args[4] else "outputs"
multithread <- if (length(args) >= 5) as.logical(args[5]) else TRUE

# Debug
# marker <- "teleo"
# path_dada_all <- "utils/teleo_clean_crabs_taxo.fasta"
# path_dada_species <- "utils/teleo_clean_crabs_species.fasta"

cat(paste0("marker is ", marker, " \n the path dada all is ", path_dada_all, " \n and the path dada species is ", path_dada_species))

# Open seqtab - abundance table output of dada2
seqtab <- readRDS(file.path(output_dir, "02_dada2", marker, "seqtab_clean.rds"))

# Assign taxo dada2
taxa <- assignTaxonomy(seqtab, path_dada_all, multithread = multithread, tryRC = TRUE)
# Sometimes issue that Species column exist at this stage - remove it
if ("Species" %in% colnames(taxa)) {
  taxa <- taxa[, !(colnames(taxa) == "Species")]
}
taxa_plus <- addSpecies(taxa, path_dada_species, verbose = TRUE, allowMultiple = TRUE)
colnames(taxa_plus) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
taxa_plus <- as.data.frame(taxa_plus)
taxa_plus$sequence <- rownames(taxa_plus)

write.csv(taxa_plus, file.path(output_dir, "02_dada2", marker, "taxonomy_dada2RDP.csv"), row.names=TRUE)

# Now link it to abundance table
# Open data 
samples <- read.csv(file.path(output_dir, "02_dada2", marker, "seqtab_clean_t.csv"))
colnames(samples)[1] <- "sequence"

# Now link 
dada2_RDP_taxo_complete <- dplyr::left_join(taxa_plus, samples) 

# Now write the output
write.csv(dada2_RDP_taxo_complete,
          file.path(output_dir, paste0(marker, "_dada2RDP.csv")),
          row.names = FALSE)