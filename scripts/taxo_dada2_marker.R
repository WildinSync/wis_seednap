# DADA2 RDP taxonomic assignment.
#
# Reads sequences from a FASTA file (query.fasta), produced by either the
# DADA2 ASV path or the SWARM OTU path -- the script is cluster-method
# agnostic by design. Applies a Wang 2007 RDP bootstrap threshold,
# cascade-nulls finer ranks below the threshold, and writes a per-sequence
# taxonomy CSV. The merge with the abundance table is performed on the
# Python side via seednap.utils.taxonomy.link_taxonomy_with_abundance so
# the output schema matches every other taxonomy method.

suppressMessages(suppressWarnings({
  library(dada2)
  library(Biostrings)
}))

args <- commandArgs(T)

if (length(args) < 5) {
  stop(paste(
    "Required args (in order):",
    "  marker, rdp_db, species_db, query_fasta, taxonomy_output_csv,",
    "  [multithread (default TRUE)], [bootstrap_threshold (default 80)]"
  ))
}

marker <- tolower(args[1])
path_dada_all <- args[2]
path_dada_species <- args[3]
query_fasta <- args[4]
output_csv <- args[5]
multithread <- if (length(args) >= 6) as.logical(args[6]) else TRUE
bootstrap_threshold <- if (length(args) >= 7) as.integer(args[7]) else 80

cat(paste0(
  "marker=", marker,
  " | rdp_db=", path_dada_all,
  " | species_db=", path_dada_species,
  " | query=", query_fasta,
  " | bootstrap_threshold=", bootstrap_threshold,
  "\n"
))

if (!file.exists(query_fasta)) {
  stop(paste("Query FASTA not found:", query_fasta))
}

# Read sequences from query.fasta. assignTaxonomy() accepts a character
# vector of sequences directly, so we don't need a sample-x-sequence matrix.
seqs <- as.character(readDNAStringSet(query_fasta))
if (length(seqs) == 0) {
  stop(paste("Query FASTA has zero sequences:", query_fasta))
}
cat(paste("Loaded", length(seqs), "sequences from query.fasta\n"))

# DADA2 RDP with bootstrap. outputBootstraps=TRUE returns a list with `tax`
# and `boot` matrices (rows=sequences, cols=Kingdom..Genus).
result <- assignTaxonomy(
  seqs, path_dada_all,
  multithread = multithread,
  tryRC = TRUE,
  outputBootstraps = TRUE
)
taxa <- result$tax
boot <- result$boot

# Apply Wang 2007 bootstrap threshold: where bootstrap < threshold, rank -> NA
threshold_mask <- boot < bootstrap_threshold
taxa[threshold_mask] <- NA

# Sometimes the RDP DB contains a Species column at this stage; remove it so
# addSpecies() controls species-level assignment.
if ("Species" %in% colnames(taxa)) {
  taxa <- taxa[, !(colnames(taxa) == "Species"), drop = FALSE]
  boot <- boot[, !(colnames(boot) == "Species"), drop = FALSE]
}

# Species-level assignment (exact match; not bootstrap-based by design)
taxa_plus <- addSpecies(taxa, path_dada_species, verbose = TRUE, allowMultiple = TRUE)
colnames(taxa_plus) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
taxa_plus <- as.data.frame(taxa_plus, stringsAsFactors = FALSE)

# Cascade null: if a coarse rank is NA, every finer rank is also NA. This
# eliminates orphan-rank rows like "kingdom=Animalia, phylum=NA, class=Mammalia".
ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
for (i in seq_len(nrow(taxa_plus))) {
  for (j in seq_along(ranks)) {
    if (is.na(taxa_plus[i, ranks[j]])) {
      if (j < length(ranks)) {
        taxa_plus[i, ranks[(j + 1):length(ranks)]] <- NA
      }
      break
    }
  }
}

# bootstrap_min: lowest bootstrap across the kept ranks (kingdom..genus only;
# species is exact-match). Acts as a per-OTU confidence summary, analogous to
# BLAST's pident. NA where no rank passed.
bootstrap_min <- apply(boot, 1, function(row) {
  passed <- row[row >= bootstrap_threshold]
  if (length(passed) == 0) NA_real_ else min(passed)
})
taxa_plus$bootstrap_min <- bootstrap_min

# rownames are the sequences themselves (DADA2 convention from assignTaxonomy)
taxa_plus$sequence <- rownames(taxa_plus)

# Reorder for downstream Python merge: sequence first, then ranks, then bootstrap
out_cols <- c("sequence", ranks, "bootstrap_min")
taxa_plus <- taxa_plus[, out_cols]

# Ensure output directory exists
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

write.csv(taxa_plus, output_csv, row.names = FALSE)
cat(paste("Wrote per-sequence taxonomy to", output_csv, "\n"))

# The merge with the abundance table is performed on the Python side via
# seednap.utils.taxonomy.link_taxonomy_with_abundance (LEFT-merge, cascade
# null already applied, contaminant flag, BLAST-compatible schema).
