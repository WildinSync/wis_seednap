# Configuration Reference

SeeDNAP uses a single YAML file per marker to configure the entire pipeline. Example configs are in `config/markers/`.

## Generating a Config

```bash
seednap init --marker teleo --output config/markers/my_marker.yaml   # minimal (required fields only)
seednap init --marker teleo --output config/markers/my_marker.yaml --full   # fully-annotated template
```

`init` writes a minimal config (just the required fields) by default; pass `--full` for the
fully-annotated reference template. A standalone minimal example also lives at
`config/markers/minimal.example.yaml`.

## Validating a Config

```bash
seednap validate config/markers/teleo.yaml
```

This checks YAML syntax, field types, required values, reports which `taxonomy.databases.<method>`
block is actually used, and flags any referenced database or `raw_data` path that is missing on disk.

> [!NOTE]
> All config models use strict validation (`extra="forbid"`), so a typo in any field name is
> rejected at load time with a clear error, including inside the `taxonomy.databases.<method>`
> blocks. You learn about a misspelled key from `seednap validate`, not hours into a run.

## Config at a glance

```
PipelineConfig
  version       free-form string (default "0.1.0")
  marker        (REQUIRED)  name + primers.{forward, reverse}
  paths         set raw_data; output/logs default
  demultiplex   off by default
  trimming      Cutadapt 2-pass (defaults)
  dada2         [ASV path only]   filter / merge / chimera / per_library
  swarm         [OTU path only]   merge / clustering / chimera
  taxonomy      (REQUIRED)  method + ONLY the matching databases.<method> block
  export        GBIF / DarwinCore (defaults)
  metrics       QC metrics (defaults)
  report        read-tracking + HTML report (on by default)
  cleaning      control decontamination (off by default)
  logging       (defaults)
  pipeline      steps: pick the dada2 OR the swarm path
```

The **dada2** and **swarm** sections are the two mutually-exclusive clustering paths; you fill
in only the one named in `pipeline.steps`. Likewise under `taxonomy` you fill only
`databases.<method>` for your chosen `method`; the other method blocks are ignored.

**Required keys** are exactly: `marker.name`, `marker.primers.forward`/`reverse`,
`taxonomy.method`, and the required path(s) in the selected database block (`blast.fasta`;
`dada2.all`; `ecotag.tree` + `fasta`; `decipher.trained`). A minimal config setting only these
loads and runs (see `config/markers/minimal.example.yaml`).

## How config merging works

Your YAML is merged over the model defaults, so you only specify what differs: any field with a
default may be omitted. Nested sections merge recursively, but **lists and scalars are replaced
wholesale, not appended**, so to change one entry of a list (e.g. `pipeline.steps` or
`taxonomy.contaminants`) you restate the whole list.

## Full Configuration Structure

### `version`

```yaml
version: "0.1.0"
```

Config format version (informational; free-form string, default `"0.1.0"`, not enforced at load time).

### `marker`

```yaml
marker:
  name: "teleo"                              # Marker name (lowercase)
  description: "Teleost fish 12S rRNA"       # Optional description
  primers:
    forward: "ACACCGCCCGTCACTCT"             # Forward primer (5' to 3')
    reverse: "CTTCCGGTACACTTACCATG"          # Reverse primer (5' to 3')
```

Primer sequences are validated for valid IUPAC DNA bases (A, C, G, T, R, Y, M, K, S, W, H, B, V, D, N). Minimum length: 10 bp.

### `paths`

```yaml
paths:
  raw_data: "/path/to/raw/fastq/files"       # Input FASTQ directory
  output: "outputs"                          # Base output directory
  logs: "logs"                               # Log files directory
```

Reference databases are not set here; each taxonomy method points at its own database under
`taxonomy.databases.<method>` (see below).

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
  min_length: 20                             # Min read length after trimming (bp)
  max_error_rate: 0.1                        # Max error rate for primer matching
  cores: 1                                   # CPU cores for cutadapt (default 1; shipped marker configs raise to 12)
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
  per_library: false                         # Learn a separate error model per sequencing
                                             #   library (grouped from the manifest seq_run_id),
                                             #   then merge + collapse. Default false = one
                                             #   pooled model. Use for multi-run datasets.
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
      # An in-band hit must also be within this many percent-identity points of
      # the best in-band hit to count toward the cascade LCA. 0 disables.
      lca_pident_delta: 1.0                  # (default: 1.0)
      # LCA algorithm: "cascade" (default), "collapsed_taxonomy", or
      # "fishbase_tiered". "cascade" is the default header-derived resolver and
      # uses top_bitscore_pct + lca_pident_delta plus the per-rank threshold_*
      # values above. "collapsed_taxonomy" is the eDNAFlow/OceanOmics
      # %identity-window collapse-to-LCA: header-based (reads the CRABS lineage
      # from the reference FASTA headers), needs no NCBI taxids/taxdump, runs
      # offline, and does NOT apply the per-rank threshold_* cascade. It is tuned
      # by lca_pid / lca_diff below. "fishbase_tiered" is not implemented and
      # raises if selected.
      lca_algorithm: "cascade"               # (default: "cascade")
      # Used only when lca_algorithm: "collapsed_taxonomy".
      lca_pid: 90.0                          # Hard %identity floor (default: 90.0)
      lca_diff: 1.0                          # Identity-window width within which
                                             #   disagreeing hits collapse to their
                                             #   LCA (default: 1.0)

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

**CRABS "NA" sentinel.** The 2025 CRABS reference DBs write the literal string
`"NA"` where a rank is unknown. SeeDNAP normalizes `"NA"`/`""`/`"nan"` to a
genuine missing rank at the BLAST formatter, so neither LCA resolver treats
`"NA"` as a taxon; missing ranks surface as `Unassigned` in the output.

### `export`

```yaml
export:
  gbif:
    enabled: true                            # Generate GBIF-format output
    add_rank: true                           # Add taxonomic rank column
    add_taxon: true                          # Add lowest taxon column
```

### `metrics`

```yaml
metrics:
  generate_plots: true                       # Collect DADA2 QC metrics + plots (DADA2 path only)
```

### `report`

Per-step read/sequence tracking and the HTML run report. **Both are generated
on every `run-pipeline` by default** (the report is built automatically at the
end of the run). See [reporting.md](reporting.md) for full details.

```yaml
report:
  read_tracking: true                        # read_tracking.{csv,txt} + warnings (default: true)
  html_report: true                          # self-contained HTML report (default: true; set false to disable)
  output_dir: null                           # base dir for report artifacts; null -> "<output>/04_report"
  warn_below_retention_pct: 30.0             # warn for samples retaining < this % of raw reads
  warn_step_loss_pct: 70.0                   # warn when a single step drops more than this %
  # Optional dataset metadata for the report's "Dataset & provenance" section:
  sample_metadata: null                      # path to per-sample (field) metadata CSV; optional
  project_metadata: null                     # path to project metadata CSV; optional
```

By default, artifacts are written to `<paths.output>/04_report/<marker>/`. Set
`output_dir` to send them elsewhere; a per-marker subdirectory is created inside
it (so `output_dir: /data/reports` for marker `teleo` writes to
`/data/reports/teleo/`). `~` is expanded and relative paths are resolved. The
whole section is optional; defaults apply when omitted.

### `cleaning`

Control decontamination of the abundance table. **Off by default.** Control
identity (which samples are negative controls, and how they associate to
extraction/PCR batches) comes from the FAIRe manifest.

```yaml
cleaning:
  enabled: false                             # Run the cleaning step (default: false)
  mode: "flag"                               # "flag" or "subtract" (default: "flag")
```

`mode: "flag"` annotates OTUs/ASVs found in negative controls without changing
counts. `mode: "subtract"` removes those control reads from the associated
samples (extraction blanks clean their extraction batch; PCR blanks clean the
whole dataset). Can also be run standalone with
`seednap clean ABUNDANCE_CSV FIELD_METADATA OUTPUT --mode flag|subtract`.

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

Valid step names: `demultiplex`, `trim`, `dada2`, `swarm`, `taxonomy`, `export`. A `clean` step
is inserted automatically after `taxonomy` when `cleaning.enabled: true`, so you do not list it
here.

## Example Configs

See `config/markers/` for complete working examples:

- `teleo.yaml` -- Teleo 12S fish marker (Namibia dataset, ligation demux)
- `mifish.yaml` -- MiFish-U 12S fish marker (Argentina dataset)
- `mam07.yaml` -- MamP007 16S mammal marker (Switzerland dataset, SWARM path)
- `mam07_dada2.yaml` -- MamP007 16S mammal marker, DADA2 ASV path
- `teleo_rhone.yaml` -- Teleo 12S (Switzerland Rhone dataset)
