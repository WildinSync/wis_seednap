# DECIPHER IdTaxa taxonomic assignment.
#
# Reads sequences from a FASTA file (query.fasta) so the script is
# cluster-method agnostic: it works after either DADA2 ASV or SWARM OTU
# clustering. Writes a per-sequence taxonomy CSV; the merge with the
# abundance table is performed on the Python side via
# seednap.utils.taxonomy.link_taxonomy_with_abundance so the output schema
# matches every other taxonomy method.

suppressMessages(suppressWarnings({
  library(Biostrings)
  library(DECIPHER)
}))

# Convert IdTaxa's nested list output into a flat DataFrame with one row
# per query sequence and columns kingdom..species (+ optional confidence).
decipher_list_to_df <- function(ids, taxonomic_ranks, confidence = FALSE) {
  rows <- lapply(seq_along(ids), function(i) {
    rec <- ids[[i]]
    out <- as.list(setNames(rep(NA_character_, length(taxonomic_ranks)), taxonomic_ranks))
    if (confidence) {
      conf_cols <- paste0("confidence_", taxonomic_ranks)
      out_conf <- as.list(setNames(rep(NA_real_, length(taxonomic_ranks)), conf_cols))
      out <- c(out, out_conf)
    }
    if (length(rec$taxon) == 0 || length(rec$rank) == 0) return(as.data.frame(out, stringsAsFactors = FALSE))
    for (j in seq_along(rec$rank)) {
      r <- tolower(rec$rank[j])
      if (r %in% taxonomic_ranks) {
        out[[r]] <- as.character(rec$taxon[j])
        if (confidence) {
          out[[paste0("confidence_", r)]] <- as.numeric(rec$confidence[j])
        }
      }
    }
    as.data.frame(out, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

args <- commandArgs(T)

if (length(args) < 4) {
  stop(paste(
    "Required args (in order):",
    "  marker, trained_classifier, query_fasta, taxonomy_output_csv,",
    "  [threshold (default 60)], [processors (default 8)]"
  ))
}

marker <- tolower(args[1])
path_decipher_trained <- args[2]
query_fasta <- args[3]
output_csv <- args[4]
threshold <- if (length(args) >= 5) as.integer(args[5]) else 60
processors <- if (length(args) >= 6) as.integer(args[6]) else 8

cat(paste0(
  "marker=", marker,
  " | trained_classifier=", path_decipher_trained,
  " | query=", query_fasta,
  " | threshold=", threshold,
  "\n"
))

if (!file.exists(query_fasta)) {
  stop(paste("Query FASTA not found:", query_fasta))
}
if (!file.exists(path_decipher_trained)) {
  stop(paste("Trained DECIPHER classifier not found:", path_decipher_trained))
}

dna <- readDNAStringSet(query_fasta)
if (length(dna) == 0) {
  stop(paste("Query FASTA has zero sequences:", query_fasta))
}
cat(paste("Loaded", length(dna), "sequences from query.fasta\n"))

trainingset <- readRDS(path_decipher_trained)
ids <- IdTaxa(
  dna, trainingset,
  strand = "both",
  processors = processors,
  verbose = TRUE,
  threshold = threshold
)

taxonomic_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
ids_df <- decipher_list_to_df(ids, taxonomic_ranks, confidence = TRUE)

# Attach the actual sequence (one row per sequence in input order)
ids_df$sequence <- as.character(dna)

# Reorder: sequence first, then ranks, then confidence columns
conf_cols <- paste0("confidence_", taxonomic_ranks)
out_cols <- c("sequence", taxonomic_ranks, conf_cols)
ids_df <- ids_df[, out_cols]

# Ensure output directory exists
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

write.csv(ids_df, output_csv, row.names = FALSE)
cat(paste("Wrote per-sequence taxonomy to", output_csv, "\n"))

# Note: the merge with the abundance table is performed on the Python side via
# seednap.utils.taxonomy.link_taxonomy_with_abundance so it goes through the
# same correctness pass as every other taxonomy method (LEFT-merge from
# abundance, cascade null, contaminant flag, BLAST-compatible schema).
