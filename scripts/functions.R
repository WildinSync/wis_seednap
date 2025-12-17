# Functions 

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

# Add the rank level & clean the species column - if the rank is genus due to multiple hits, put NA in species
add_rank_dada <- function(file){
    file_out <- file %>%
        mutate(rank = case_when(
           !grepl("/", species) & !is.na(species) ~ "species", 
           !is.na(genus) & (grepl("/", species) | is.na(species)) ~ "genus",
           !is.na(family) & is.na(genus) ~ "family",
           TRUE ~ "higher"
        )) %>%
        mutate(species = case_when(
            rank != "species" & grepl("/", species) ~ NA,
            TRUE ~ species
        ))
    return(file_out)
}

# Add the taxon information - this is the lowest taxonomic information available
add_taxon_dada <- function(file){
    file %>%
        mutate(taxon = case_when(
            rank == "species" ~ species, 
            rank == "genus" ~ genus, 
            rank == "family" ~ family, 
            rank == "higher" & !is.na(order) ~ order, 
            rank == "higher" & !is.na(class) ~ class, 
            rank == "higher" & !is.na(phylum) ~ phylum, 
            rank == "higher" & !is.na(kingdom) ~ kingdom, 
            TRUE ~ NA
        ))
}
## END OF FUNCTIONS 

# ENF OF FUNCTIONS 
# ------------------------- # 