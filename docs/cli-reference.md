# CLI Reference

## Global Options

```
seednap --version           Show version
seednap --help              Show help
seednap -v <command>        Verbose output (DEBUG level)
seednap -q <command>        Quiet mode (errors only)
```

---

## `run-pipeline`

Run the complete pipeline end-to-end from a YAML config.

```
seednap run-pipeline CONFIG [OPTIONS]
```

| Option | Description |
|---|---|
| `--resume` | Resume from previous run (skip completed steps) |
| `--state-file PATH` | Custom state file path |
| `--stop-on-error / --continue-on-error` | Stop or continue on first error (default: stop) |

```bash
seednap run-pipeline config/markers/teleo.yaml
seednap run-pipeline config/markers/teleo.yaml --resume
seednap run-pipeline config/markers/teleo.yaml --continue-on-error
```

---

## `trim`

Two-pass primer trimming with Cutadapt.

```
seednap trim INPUT_DIR [OPTIONS]
```

| Option | Required | Description |
|---|---|---|
| `--forward-primer TEXT` | Yes | Forward primer sequence (5' to 3') |
| `--reverse-primer TEXT` | Yes | Reverse primer sequence (5' to 3') |
| `-o, --output-dir PATH` | Yes | Output directory for trimmed reads |
| `-c, --cores INTEGER` | No | Number of CPU cores |

```bash
seednap trim /path/to/raw/fastq \
  --forward-primer ACACCGCCCGTCACTCT \
  --reverse-primer CTTCCGGTACACTTACCATG \
  -o outputs/01_trim/teleo
```

---

## `swarm`

SWARM OTU clustering on trimmed reads.

```
seednap swarm MARKER TRIMMED_READS_DIR [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-o, --output-dir PATH` | `outputs/` | Base output directory |
| `-d, --distance INTEGER` | `1` | SWARM distance threshold |
| `-t, --threads INTEGER` | `4` | CPU threads |
| `--no-fastidious` | | Disable singleton refinement |
| `--no-chimera-filter` | | Skip chimera detection |

```bash
seednap swarm teleo /path/to/trimmed -o outputs -d 1 -t 8
```

---

## `dada2`

DADA2 ASV processing on trimmed reads.

```
seednap dada2 MARKER TRIMMED_READS_DIR [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-o, --output-dir PATH` | `outputs/` | Base output directory |
| `--max-ee FLOAT` | `2.0` | Maximum expected errors for filtering |
| `--trunc-q INTEGER` | `11` | Truncate at first base with quality below this |
| `--min-overlap INTEGER` | `20` | Minimum overlap for merging paired reads |
| `--assign-taxonomy` | | Run DADA2 taxonomic assignment (requires `--rdp-db` and `--species-db`) |
| `--rdp-db PATH` | | RDP-formatted taxonomy database |
| `--species-db PATH` | | Species-level database |

```bash
seednap dada2 teleo /path/to/trimmed -o outputs --max-ee 2.0 --trunc-q 11
```

> DADA2 can learn a separate error model per sequencing library (grouped
> from the FAIRe manifest's `seq_run_id`) and merge the per-library
> sequence tables, via the `dada2.per_library` config field (default
> `false` = single pooled error model). See
> [configuration.md](configuration.md). This is a YAML-only knob; it is
> not exposed on the `dada2` CLI command.

---

## `clean`

Decontaminate an abundance table against its negative controls. Control
identity (extraction vs PCR blanks, extraction batches) is derived from a
FAIRe manifest migrated from the field metadata.

```
seednap clean ABUNDANCE_CSV FIELD_METADATA OUTPUT [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `--mode {flag\|subtract}` | `flag` | `flag` annotates control-detected OTUs/ASVs without changing counts; `subtract` removes control reads (extraction blanks clean their `extraction_ID` batch, PCR blanks clean the whole dataset) |
| `--project-metadata PATH` | | Project metadata CSV (marker -> target_gene) |
| `--id-col TEXT` | first column | OTU/ASV identifier column in the abundance table |
| `--report PATH` | `<output stem>_report.csv` | Per-sample cleaning report CSV |

```bash
seednap clean outputs/02_swarm/teleo/otu_table.csv \
  metadata/metadata_field_my_dataset.csv \
  outputs/teleo_otu_clean.csv --mode subtract
```

In a full run this is the `cleaning:` config section
(`enabled`, default `false`; `mode`). See
[configuration.md](configuration.md).

---

## `blast`

BLAST taxonomic assignment with LCA resolution.

```
seednap blast QUERY_FASTA REF_FASTA ASV_COUNT [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-o, --output PATH` | auto | Output CSV file path |
| `--perc-identity FLOAT` | `80.0` | Minimum percent identity |
| `--qcov-hsp-perc FLOAT` | `80.0` | Minimum query coverage per HSP |
| `--evalue FLOAT` | `1e-25` | Maximum e-value |
| `--task {megablast\|blastn\|dc-megablast\|blastn-short}` | `megablast` | blastn task |
| `--threshold-species FLOAT` | `99.0` | Percent identity for species assignment |
| `--threshold-genus FLOAT` | `96.0` | Percent identity for genus assignment |
| `--threshold-family FLOAT` | `90.0` | Percent identity for family assignment |
| `--threshold-order FLOAT` | `80.0` | Percent identity for order assignment |
| `--threshold-class FLOAT` | `70.0` | Percent identity for class assignment |
| `--lca-algorithm {cascade\|collapsed_taxonomy}` | `cascade` | LCA algorithm (see below) |
| `--top-bitscore-pct FLOAT` | `10.0` | cascade LCA: include hits within this % of the best bitscore (MEGAN-LR band) |
| `--lca-pident-delta FLOAT` | `1.0` | cascade LCA: in-band hits must be within this %id of the best in-band hit |
| `--lca-pid FLOAT` | `90.0` | collapsed_taxonomy only: hard %identity floor for hits |
| `--lca-diff FLOAT` | `1.0` | collapsed_taxonomy only: identity-window width collapsed to the LCA |

**LCA algorithms.** `cascade` (default) applies the per-rank thresholds
above within a MEGAN-LR bitscore band (`--top-bitscore-pct`) above a
`--lca-pident-delta` floor. `collapsed_taxonomy` is the
eDNAFlow/OceanOmics %identity-window collapse-to-LCA: it ignores the
per-rank thresholds and instead collapses disagreeing hits within
`--lca-diff` %id of one another to their LCA, above the `--lca-pid` hard
floor. It is header-based (reads the CRABS lineage from the reference
FASTA headers), needs no NCBI taxids/taxdump, and runs offline.

```bash
seednap blast outputs/02_swarm/teleo/query.fasta \
  /path/to/reference.fasta \
  outputs/02_swarm/teleo/otu_table.csv \
  -o outputs/teleo_blast.csv \
  --evalue 1e-10 --threshold-species 100
```

Reference headers may carry the literal string `NA` for an unknown rank
(2025 CRABS DBs); SeeDNAP normalizes `NA`/empty/`nan` to a genuine
missing rank, which surfaces as `Unassigned` in the export rather than a
taxon named `NA`.

See [configuration.md](configuration.md#taxonomy) and
[taxonomy-methods.md](taxonomy-methods.md#blast--lca-recommended).

---

## `assign-taxonomy`

Generic taxonomic assignment supporting all methods.

```
seednap assign-taxonomy {blast|dada2|ecotag|decipher} MARKER QUERY_FASTA ASV_COUNT_CSV [OPTIONS]
```

Each method requires specific database options:

| Method | Required Options |
|---|---|
| `blast` | `--reference-fasta PATH` |
| `dada2` | `--rdp-db PATH`, `--species-db PATH` |
| `ecotag` | `--taxonomy-db PATH`, `--reference-db PATH` |
| `decipher` | `--trained-classifier PATH` |

Additional options:

| Option | Default | Description |
|---|---|---|
| `-o, --output-dir PATH` | `outputs/` | Base output directory |
| `--threshold-species FLOAT` | `99.0` | Species %ID threshold (BLAST) |
| `--threshold-genus FLOAT` | `96.0` | Genus %ID threshold (BLAST) |
| `--threshold-family FLOAT` | `90.0` | Family %ID threshold (BLAST) |
| `--threshold-order FLOAT` | `80.0` | Order %ID threshold (BLAST) |
| `--threshold-class FLOAT` | `70.0` | Class %ID threshold (BLAST) |
| `--lca-algorithm {cascade\|collapsed_taxonomy}` | `cascade` | BLAST LCA algorithm (see `blast` above) |
| `--top-bitscore-pct FLOAT` | `10.0` | cascade LCA bitscore band as % of best hit (BLAST) |
| `--lca-pident-delta FLOAT` | `1.0` | cascade LCA: in-band hits within this %id of the best in-band hit (BLAST) |
| `--lca-pid FLOAT` | `90.0` | collapsed_taxonomy: hard %identity floor (BLAST) |
| `--lca-diff FLOAT` | `1.0` | collapsed_taxonomy: identity-window width collapsed to the LCA (BLAST) |
| `--confidence-threshold INT` | `60` | Confidence threshold (DECIPHER) |
| `-c, --processors INTEGER` | `8` | CPU cores |

---

## `format-gbif`

Convert taxonomy results to GBIF long format.

```
seednap format-gbif INPUT_FILE [OPTIONS]
```

| Option | Required | Description |
|---|---|---|
| `-f, --format {dada2\|ecotag\|blast\|decipher}` | Yes | Input format type |
| `-o, --output PATH` | No | Output path (default: auto-generated) |

```bash
seednap format-gbif outputs/teleo_blast.csv -f blast -o outputs/teleo_gbif.csv
```

---

## `create-gbif`

Build a full DarwinCore-compliant GBIF occurrence CSV.

```
seednap create-gbif TAXONOMY_RESULTS SAMPLE_METADATA PROJECT_METADATA OUTPUT [OPTIONS]
```

| Option | Description |
|---|---|
| `--summarise-pcr / --no-summarise-pcr` | Aggregate PCR replicates per sample |
| `--skip-enrichment` | Skip NCBI/WORMS taxonomy enrichment |

Requires `NCBI_API_KEY` in `.env` for taxonomy enrichment. See `.env.example`.

```bash
seednap create-gbif outputs/teleo_gbif.csv metadata/samples.csv metadata/project.csv outputs/teleo_darwincore.csv
```

---

## `demultiplex`

Demultiplex ligation-based libraries.

```
seednap demultiplex RAW_READS_DIR LIBRARY_NAME METADATA_CSV [OPTIONS]
```

| Option | Required | Description |
|---|---|---|
| `-f, --forward-primer TEXT` | Yes | Forward primer sequence |
| `-r, --reverse-primer TEXT` | Yes | Reverse primer sequence |
| `-o, --output-dir PATH` | Yes | Output base directory |
| `-c, --cores INTEGER` | No | CPU cores |
| `--no-gunzip` | No | Keep output files gzipped |

---

## `manifest`

Build (and optionally validate) a canonical FAIRe sample manifest from the
lab's existing CSVs. The manifest is the source of truth for sample ->
library grouping (used by `dada2.per_library`) and control identity (used
by `clean`). Standalone and read-only on its inputs.

```
seednap manifest FIELD_METADATA [OPTIONS]
```

`FIELD_METADATA` is the per-sample field metadata CSV
(`metadata_field_*.csv`), or a legacy demux lab CSV
(`metadata_lab_*.csv`) carrying library/tag columns.

| Option | Default | Description |
|---|---|---|
| `--project-metadata PATH` | | Project metadata CSV (supplies the marker -> target_gene/assay_name) |
| `--lab-metadata PATH` | | Legacy demux metadata CSV (supplies seq_run_id/library and tag barcodes) |
| `--seq-run-id TEXT` | | Sequencing-run id for the whole dataset (overrides lab/derived value) |
| `--target-gene TEXT` | | Marker / target_gene (overrides the project metadata) |
| `--date-order {ymd\|dmy\|mdy}` | | Force the eventDate field order for genuinely-ambiguous dotted dates (otherwise such files raise rather than be guessed) |
| `-o, --output PATH` | | Write the canonical manifest CSV here |
| `--abundance PATH` | | Validate the manifest's eventIDs against this abundance/OTU table |
| `--strict` | | Raise if the abundance table has sample columns absent from the manifest (default: warn) |

```bash
seednap manifest metadata/metadata_field_my_dataset.csv \
  --project-metadata metadata/metadata_proj_my_dataset.csv \
  --abundance outputs/02_swarm/teleo/otu_table.csv \
  -o metadata/manifest_my_dataset.csv
```

---

## `init`

Generate an example configuration file.

```
seednap init [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-m, --marker TEXT` | `teleo` | Marker name |
| `-o, --output PATH` | `config/markers/example.yaml` | Output path |
| `--minimal / --full` | `--minimal` | Required-fields-only config (default) or the fully-annotated reference template |
| `-f, --force` | | Overwrite existing file |

---

## `validate`

Validate a YAML configuration file. Checks syntax, field types, and required fields (including
typos inside `taxonomy.databases.<method>`, which are rejected at load). The summary also reports
which database block is live for the selected method and flags any referenced database or
`raw_data` path that is missing on disk.

```
seednap validate CONFIG
```

---

## `report`

Build the read/sequence-tracking report (and optionally the HTML run report)
from an existing run's outputs. `run-pipeline` already generates both at the end
of every run by default (see the `report:` config block); this command is for
**regenerating** them from outputs that already exist.

```
seednap report MARKER [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-o, --output-dir PATH` | `outputs/` | Base output directory of the run |
| `--html` | | Also generate the self-contained HTML run report |
| `--warn-retention FLOAT` | `30.0` | Warn for samples retaining below this % of raw reads |
| `--warn-step-loss FLOAT` | `70.0` | Warn when a single step drops more than this % of a sample's reads |
| `--field-metadata PATH` | auto | Field metadata CSV (location, dates, sites) for the Dataset section |
| `--project-metadata PATH` | auto | Project metadata CSV (recorder, sequencing, reference DB) |
| `--log-file PATH` | auto | Run log to embed (colorized) in the HTML report's Run-log section; auto-located under `logs/` if omitted |

```bash
seednap report teleo -o outputs --html \
  --field-metadata metadata/metadata_field_my_dataset.csv \
  --project-metadata metadata/metadata_proj_my_dataset.csv
```

Writes `outputs/04_report/<marker>/read_tracking.{csv,txt}` and, with `--html`,
`report.html`. See [reporting.md](reporting.md) for details.

---

## `monitor`

Summarise a finished or in-progress run from its state JSON. Prints a
per-step status/duration table plus the read-tracking headline (raw ->
final reads, mean retention, warnings), and writes a
`monitoring_summary.csv` when per-sample counts are present. Standalone
and read-only; regenerable any time without a re-run.

```
seednap monitor MARKER [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-o, --output-dir PATH` | `outputs` | Base output directory of the run |
| `--state-file PATH` | `<output-dir>/.<marker>_state.json` | Run state JSON |

```bash
seednap monitor teleo -o outputs
```

---

## `version`

Show detailed version information.

```
seednap version
```
