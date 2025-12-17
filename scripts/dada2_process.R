# dada2 code
options(bitmapType = "cairo")

suppressMessages(suppressWarnings({
  library(dada2); packageVersion("dada2")
  library(Biostrings)
  library(DECIPHER)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
}))

print("[INFO] Start DADA2 processing.")

# Get the arguments
args <- commandArgs(T)

# Check file exist
if (length(args)==0) {message("Please enter an valid file name.")} else
{  
  marker <- tolower(args[1]) #Must be character 
}
# Debug
# marker <- "teleo"

# ---------------------------------- # 
# FUNCTIONS 

# Convert df to fasta file
df_to_fasta <- function(file, output_file_path){
  fa <- character(2 * nrow(file))
  fa[c(TRUE, FALSE)] = sprintf(">%s", file[,1])
  fa[c(FALSE, TRUE)] = as.character(file[,2])
  writeLines(fa, output_file_path)
}

extract_pattern_samplename <- function(input_string) {
  sub(".*/([^/]+)\\.[Rr][12].*", "\\1", input_string)
}

# END of FUNCTIONS
# ------------------------------ #

# Create directory
dir.create(paste0("outputs/02_dada2/", marker, "/QC/"))

# File parsing
pathFR <- paste0("outputs/01_trim/", marker, "/") # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
filtpathFR <- file.path(pathFR, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
fastqFs <- sort(list.files(pathFR, pattern="R1.fastq"))
fastqRs <- sort(list.files(pathFR, pattern="R2.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# Explore quality 
fnFs <- sort(list.files(pathFR, pattern="R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(pathFR, pattern="R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)

# Check if file is not empty
valid_indices <- which(file.exists(fnFs) & file.exists(fnRs) & file.info(fnFs)$size > 0 & file.info(fnRs)$size > 0)

# Generate QC images
invisible(mclapply(valid_indices, function(i) {
  name_sample <- extract_pattern_samplename(fnFs[i])
  
  # Try-catch to avoid errors stopping the loop
  p_f <- try(plotQualityProfile(fnFs[i]), silent = TRUE)
  p_r <- try(plotQualityProfile(fnRs[i]), silent = TRUE)
  
  # Check if plots were generated successfully
  if (!inherits(p_f, "try-error") & !inherits(p_r, "try-error")) {
    p_i <- p_f + p_r
    png(paste0("outputs/02_dada2/", marker, "/QC/", name_sample, "_dada2QC.png"), width = 700, height = 700)
    print(p_f + p_r)
    dev.off()
  } else {
    message("Skipping ", name_sample, " due to an issue with plotQualityProfile.")
  }
}, mc.cores = detectCores() - 30))

# File parsing
filtpathFR <- paste0("outputs/01_trim/", marker, "/filtered")

# Filter
filterAndTrim(fwd=file.path(pathFR, fastqFs), filt=file.path(filtpathFR, fastqFs),
              rev=file.path(pathFR, fastqRs), filt.rev=file.path(filtpathFR, fastqRs),
              maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=FALSE, verbose=TRUE, multithread=TRUE)

# Generate QC images - after filtering
invisible(mclapply(seq_along(paste0(filtpathFR, "/", fastqFs)), function(i) {
  name_sample <- extract_pattern_samplename(paste0(filtpathFR, "/", fastqFs)[i])
  
  # Try-catch to avoid errors stopping the loop
  p_f <- try(plotQualityProfile(paste0(filtpathFR, "/", fastqFs)[i]), silent = TRUE)
  p_r <- try(plotQualityProfile(paste0(filtpathFR, "/", fastqFs)[i]), silent = TRUE)
  
  # Check if plots were generated successfully
  if (!inherits(p_f, "try-error") & !inherits(p_r, "try-error")) {
    p_i <- p_f + p_r
    png(paste0("outputs/02_dada2/", marker, "/QC/", name_sample, "_cleaned_dada2QC.png"), width = 700, height = 700)
    print(p_f + p_r)
    invisible(dev.off())
  } else {
    message("Skipping ", name_sample, " due to an issue with plotQualityProfile.")
  }
}, mc.cores = detectCores() - 30))

# FIltered reads
filtFs <- list.files(filtpathFR, pattern="R1.fastq", full.names = TRUE)
filtRs <- list.files(filtpathFR, pattern="R2.fastq", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "\\."), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "\\."), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = 20)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, paste0("outputs/02_dada2/", marker, "/seqtab.rds"))

# Merge multiple runs (if necessary)
st1 <- readRDS(paste0("outputs/02_dada2/", marker, "/seqtab.rds"))
# Remove chimeras
seqtab <- removeBimeraDenovo(st1, method="consensus", multithread=TRUE)
write.csv(seqtab, paste0("outputs/02_dada2/", marker, "/seqtab_clean.csv"), row.names = TRUE)
saveRDS(seqtab,  paste0("outputs/02_dada2/", marker, "/seqtab_clean.rds")) 
#saveRDS(t(seqtab),  paste0("outputs/02_dada2/", marker, "/seqtab_clean_t.rds"))

write.csv(t(seqtab), paste0("outputs/02_dada2/", marker, "/seqtab_clean_t.csv"), row.names = TRUE)

# Output a fasta for ecotag and blast
df_seq <- data.frame(sequence = colnames(seqtab))
df_seq <- df_seq %>%
  mutate(ASV_n = paste0("ASV", row_number())) %>%
  select(ASV_n, sequence)
write.csv(df_seq, paste0("outputs/02_dada2/", marker, "/corresp_seq.csv"), row.names = FALSE)

# fasta convertion
df_to_fasta(file = df_seq, 
  output_file_path = paste0("outputs/02_dada2/", marker, "/query.fasta"))