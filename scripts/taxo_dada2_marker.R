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

# Drop any Species column the RDP DB emitted at this stage so addSpecies() owns species.
if ("Species" %in% colnames(taxa)) {
  keep <- colnames(taxa) != "Species"
  taxa <- taxa[, keep, drop = FALSE]
  boot <- boot[, keep, drop = FALSE]
}

# Species-level assignment FIRST, on the UNMASKED genus. addSpecies() is exact-match
# (100% identity) and is NOT bootstrap-based by design. Running it before the bootstrap
# mask is deliberate: masking the genus to NA first makes addSpecies's genus-consistency
# check fail, so a confident exact species hit (e.g. a 100%-identity Cervus ASV whose RDP
# *genus* bootstrap is only ~67) would be silently discarded. The old pipeline assigned
# these (it used minBoot=50 and added species on the kept genus); we preserve that for
# exact matches while keeping a stricter threshold for everything else.
taxa_plus <- addSpecies(taxa, path_dada_species, verbose = TRUE, allowMultiple = TRUE)
colnames(taxa_plus) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
taxa_plus <- as.data.frame(taxa_plus, stringsAsFactors = FALSE)

# Apply the Wang 2007 bootstrap threshold to kingdom..genus, EXCEPT on rows where
# addSpecies made an exact match: that 100%-identity match confirms the whole lineage and
# overrides the RDP bootstrap (so the species and its implied genus/family/... are kept).
boot_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
species_matched <- !is.na(taxa_plus$species)
for (j in seq_along(boot_ranks)) {
  drop_rank <- (boot[, j] < bootstrap_threshold) & !species_matched
  taxa_plus[drop_rank, boot_ranks[j]] <- NA
}

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

# bootstrap_min: lowest RDP bootstrap across the ranks actually kept (kingdom..genus;
# species is exact-match). A per-OTU confidence summary analogous to BLAST's pident. NA
# when no kingdom..genus rank was kept. Computed over the KEPT ranks (not just those that
# passed the threshold) so an exact-match-overridden rank reports its real, lower boot.
bootstrap_min <- vapply(seq_len(nrow(taxa_plus)), function(i) {
  kept <- !is.na(taxa_plus[i, boot_ranks])
  if (!any(kept)) NA_real_ else min(boot[i, kept])
}, numeric(1))
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
