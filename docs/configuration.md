# Configuration Reference

SeeDNAP uses a single YAML file per marker to configure the entire pipeline. Example configs are in `config/markers/`.

## Generating a Config

```bash
seednap init --marker teleo --output config/markers/teleo.yaml
```

## Validating a Config

```bash
seednap validate config/markers/teleo.yaml
```

This checks YAML syntax, field types, required values, and reports any errors.

**Note:** All config models use strict validation (`extra="forbid"`). Typos in field names will be rejected at load time with a clear error message.

## Full Configuration Structure

### `version`

```yaml
version: "0.1.0"
```

Config format version. Must match the pipeline version.

### `marker`

```yaml
marker:
  name: "teleo"                              # Marker name (lowercase)
  description: "Teleost fish 12S rRNA"       # Optional description
  primers:
    forward: "ACACCGCCCGTCACTCT"             # Forward primer (5' to 3')
    reverse: "CTTCCGGTACACTTACCATG"          # Reverse primer (5' to 3')
    name: "Teleo"                            # Primer set name
    target: "12S rRNA"                       # Target gene region
    amplicon_length: [40, 100]               # Expected amplicon range [min, max]
```

Primer sequences are validated for valid IUPAC DNA bases (A, C, G, T, R, Y, M, K, S, W, H, B, V, D, N). Minimum length: 10 bp.

### `paths`

```yaml
paths:
  raw_data: "/path/to/raw/fastq/files"       # Input FASTQ directory
  output: "outputs"                          # Base output directory
  logs: "logs"                               # Log files directory
  references: "/path/to/reference/databases" # Reference databases
```

Relative paths are resolved to absolute paths. `~` is expanded. Output and log directories are created automatically.

### `demultiplex`

```yaml
demultiplex:
  enabled: false                             # Enable/disable demultiplexing
  protocol: "none"                           # "ligation", "standard", or "none"
  metadata: "/path/to/metadata.csv"          # Required if enabled
  skip: false                                # True if raw inputs are already
                                             #   demultiplexed (one FASTQ per
                                             #   sample); the orchestrator
                                             #   records the step as skipped
                                             #   instead of running it.
  max_sample_failure_rate: 0.5               # Abort the demultiplex step if
                                             #   more than this fraction of
                                             #   samples fail; otherwise log
                                             #   the failures and continue.
```

The `ligation` protocol processes one bad sample at a time with
`try`/`except` and only fails the whole library when the per-sample
failure rate crosses `max_sample_failure_rate` (default 50%). The
`standard` protocol is reserved for future work and currently raises a
`NotImplementedError` with a pointer to the ligation path.

### `trimming`

```yaml
trimming:
  tool: "cutadapt"                           # Only cutadapt supported
  min_length: 20                             # Min read length after trimming (bp)
  max_error_rate: 0.1                        # Max error rate for primer matching
  cores: 12                                  # CPU cores for cutadapt
  discard_untrimmed: true                    # Discard reads without primers
  overlap: 3                                 # Min overlap for primer detection (bp)
```

### `swarm`

```yaml
swarm:
  merge:
    fastq_maxdiffs: 10                       # Max differences in overlap region
    fastq_minovlen: 10                       # Min overlap length for merging
    allow_stagger: false                     # Allow staggered read merging
  clustering:
    d: 1                                     # SWARM distance threshold
    fastidious: true                         # Refine singletons
    boundary: 3                              # Min mass for large OTUs (fastidious)
    threads: 4                               # CPU threads
  chimera:
    method: "denovo"                         # "denovo" or "none"
  min_sequence_length: 20                    # Min sequence length after merging
```

### `dada2`

```yaml
dada2:
  filter:
    max_ee: 2.0                              # Maximum expected errors
    trunc_q: 11                              # Truncate at first base with Q <= trunc_q
    max_n: 0                                 # Max N bases allowed
    rm_phix: true                            # Remove PhiX contamination
    min_len: null                            # Optional min read length
    max_len: null                            # Optional max read length
  merge:
    min_overlap: 20                          # Min overlap for merging (bp)
    max_mismatch: 0                          # Max mismatches in overlap
  chimera:
    method: "consensus"                      # "consensus", "pooled", or "none"
  pool: false                                # Pool samples for denoising
  multithread: true                          # Use multithreading
```

### `taxonomy`

```yaml
taxonomy:
  method: "blast"                            # "blast", "dada2", "decipher", "ecotag"

  # Marker-level contaminant species. Matched against the assigned `species`
  # column in the post-processor; matching rows get
  # `is_contaminant_candidate=True`. Rows are NEVER deleted -- downstream
  # decides what to do with the flag. Use the underscore-separated CRABS
  # format. (Whitmore et al. 2023, Nat. Ecol. Evol.)
  contaminants:
    - "Homo_sapiens"
    - "Bos_taurus"
    - "Sus_scrofa"

  databases:
    blast:
      fasta: "/path/to/blast_db.fasta"       # Reference FASTA (required)
      perc_identity: 80.0                    # Min percent identity (default: 80.0)
      qcov_hsp_perc: 80.0                    # Min query coverage per HSP (default: 80.0)
      evalue: 1.0e-25                        # Max e-value (default: 1.0e-25)
      max_target_seqs: 5                     # Max hits retained (default: 5)
      task: "megablast"                      # blastn task; "megablast" (default)
                                             #   for short, high-identity vertebrate
                                             #   amplicons against curated DBs;
                                             #   "blastn" for divergent references.
      # Per-rank cascade thresholds. If percent identity falls below the
      # threshold for a rank, that rank AND every finer rank are nulled
      # (cascade null), so the output never shows orphan ranks like
      # "kingdom set, phylum None, class Mammalia". Defaults follow
      # Pappalardo 2025 (Methods Ecol. Evol. 16:2380-2394) with rRNA-marker
      # tweaks (family raised vs eDNAFlow).
      threshold_species: 99.0                # (default: 99.0)
      threshold_genus: 96.0                  # (default: 96.0)
      threshold_family: 90.0                 # (default: 90.0)
      threshold_order: 80.0                  # (default: 80.0)
      threshold_class: 70.0                  # (default: 70.0)
      # MEGAN-LR style top-bitscore band: hits within this percent of the
      # best bitscore are considered together for LCA resolution.
      # 0.0 = exact ties only.
      top_bitscore_pct: 10.0                 # (default: 10.0)

    dada2:
      all: "/path/to/dada2_all.fasta"        # RDP-format database (required)
      species: "/path/to/dada2_species.fasta" # Species-level database (optional)
      # Naive Bayesian (Wang 2007) bootstrap confidence threshold. Below
      # this value, the rank is nulled AND every finer rank cascades to
      # null. 80% is the published recommendation for short rRNA reads.
      bootstrap_threshold: 80                # (default: 80)

    ecotag:
      tree: "/path/to/taxonomy/"             # NCBI taxonomy tree
      fasta: "/path/to/reference.fasta"      # Reference sequences

    decipher:
      trained: "/path/to/trained.rds"        # Trained classifier
      threshold: 60                          # Confidence threshold (0-100)
      processors: 8                          # CPU cores
```

You only need to provide the database section for the method you selected.

**Defaults vs CLI shortcuts.** The pipeline (YAML) defaults above are the
field-standard cascade thresholds used by the production runs. The
standalone `seednap blast` and `seednap assign-taxonomy` commands keep
their historical CLI defaults (`threshold_species=98.0`,
`threshold_family=86.5`); pass `--threshold-*` explicitly or run via
`run-pipeline` with a YAML config to use the cascade defaults.

### `export`

```yaml
export:
  formats:
    - "csv"
  gbif:
    enabled: true                            # Generate GBIF-format output
    add_rank: true                           # Add taxonomic rank column
    add_taxon: true                          # Add lowest taxon column
```

### `metrics`

```yaml
metrics:
  generate_plots: true                       # Generate QC plots
  plot_format: "png"                         # "png", "pdf", or "svg"
  metrics:
    - "read_counts"
    - "quality_scores"
    - "length_distribution"
```

### `logging`

```yaml
logging:
  level: "INFO"                              # "DEBUG", "INFO", "WARNING", "ERROR"
  format: "detailed"                         # "simple", "detailed", "json"
  file: true                                 # Write to log file
  console: true                              # Write to console
```

### `pipeline`

```yaml
pipeline:
  steps:
    - "trim"
    - "swarm"                                # or "dada2"
    - "taxonomy"
  skip: []                                   # Steps to skip (e.g., ["trim"])
```

Valid step names: `demultiplex`, `trim`, `dada2`, `swarm`, `taxonomy`, `export`.

## Example Configs

See `config/markers/` for complete working examples:

- `teleo.yaml` -- Teleo 12S fish marker (Namibia dataset, ligation demux)
- `mifish.yaml` -- MiFish-U 12S fish marker (Argentina dataset)
- `mam07.yaml` -- MamP007 16S mammal marker (Switzerland dataset)
- `teleo_rhone.yaml` -- Teleo 12S (Switzerland Rhone dataset)
