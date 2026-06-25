# ⌨️ CLI Reference

<img src="../media/divider.svg" width="100%" alt="">

Complete reference for every `seednap` command, its arguments, and its options.

SeeDNAP exposes one command group with several subcommands. `run-pipeline` runs the whole pipeline from a YAML config; the remaining commands run or regenerate a single stage standalone. For per-stage behavior and algorithms see [pipeline-steps.md](pipeline-steps.md); for every config key see [configuration.md](configuration.md).

A few terms recur below. An **ASV** (amplicon sequence variant) is an exact denoised sequence produced by DADA2; an **OTU** (operational taxonomic unit) is a cluster of similar sequences produced by SWARM. Both are rows in the final abundance table (one feature per row, one sample per column). The two are alternative, mutually exclusive feature paths.

The full command surface, grouped by purpose:

<p align="center">
  <img src="../media/cli.svg" width="100%" alt="seednap commands grouped by purpose: run, individual stages, GBIF export, prepare data, inspect a run">
</p>

How a command reads, global options before the command, the argument and `--options` after:

<p align="center">
  <img src="../media/command-anatomy.svg" width="100%" alt="anatomy of a seednap command: program, global option, command, argument, option">
</p>

## 📂 Where commands write

`run-pipeline` writes under the canonical output tree rooted at `paths.output` (default `outputs/`):

| Path | Contents |
|---|---|
| `<output>/01_trim/<marker>/` | Trimmed reads (and `demux/` for the demultiplex step) |
| `<output>/02_dada2/<marker>/` or `<output>/02_swarm/<marker>/` | ASV or OTU tables |
| `<output>/03_taxo/<marker>/` | Taxonomy intermediates (BLAST TSV, etc.) |
| `<output>/04_report/<marker>/` | `read_tracking.{csv,txt}`, `step_summary.csv`, `report.html` |
| `<output>/<marker>_<method>.csv` | Final merged taxonomy + abundance table |
| `<output>/<marker>_<method>_cleaned.csv` | Cleaned table, written only when the `clean` step runs |
| `<output>/<marker>_<method>_gbif.csv` | Final GBIF long-format table |
| `<output>/.<marker>_state.json` | Pipeline state JSON (drives `--resume`) |
| `<output>/.<marker>_config.snapshot.yaml` | Effective merged config recorded for the run |

`<method>` is the `taxonomy.method` token, except the DADA2 taxonomy table uses `dada2RDP` (so `<marker>_dada2RDP.csv`; the GBIF and cleaned files still use the plain `dada2` token). The standalone commands write where you point `-o`/`--output`; the example paths below assume the default tree.

Every `run-pipeline` run is reproducible from its own outputs. The effective merged config (your YAML plus all defaults that were filled in) is snapshotted to `<output>/.<marker>_config.snapshot.yaml`, and the running SeeDNAP version is stamped into the state JSON. On `--resume`, a version mismatch between the recorded and running version emits a `[WARN]` rather than silently mixing outputs from two versions.

## 🌐 Global options

`-v`/`--verbose` and `-q`/`--quiet` are options on the top-level group.

| Option | Description |
|---|---|
| `--version` | Show the version and exit |
| `--help` | Show help and exit |
| `-v, --verbose` | Verbose logging (DEBUG level) |
| `-q, --quiet` | Quiet mode (errors only) |

```bash
seednap --version
seednap -v run-pipeline config/markers/teleo.yaml
```

> [!WARNING]
> `-v` and `-q` must come BEFORE the subcommand (`seednap -v run-pipeline ...`). Placed after the subcommand they are not recognized and the command fails.

## 🚀 `run-pipeline`

Run the complete pipeline end-to-end from a YAML config.

```
seednap run-pipeline CONFIG [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `--resume` | off | Skip steps already marked done in the state JSON |
| `--state-file PATH` | `<output>/.<marker>_state.json` | State JSON path |
| `--stop-on-error / --continue-on-error` | stop | Stop or continue on the first failed step |

```bash
seednap run-pipeline config/markers/teleo.yaml
seednap run-pipeline config/markers/teleo.yaml --resume
seednap run-pipeline config/markers/teleo.yaml --continue-on-error
```

`--resume` reads the state JSON (`--state-file`, defaulting to `<output>/.<marker>_state.json`) and skips steps recorded as completed, rerunning from the first incomplete or failed step.

`run-pipeline` runs a preflight before any compute: if a referenced database or `raw_data` path is missing on disk, it exits 1 instead of failing mid-run. Run `seednap validate CONFIG` first to catch this earlier.

> [!WARNING]
> If `pipeline.steps` lists `demultiplex` with any `demultiplex.protocol` other than `ligation` (including the default `none`), the config is REJECTED AT LOAD with a `ValueError`, before any step runs. Only the `ligation` protocol is implemented. Either set `demultiplex.protocol: ligation`, or, if your reads are already demultiplexed, drop `demultiplex` from `pipeline.steps` so the pipeline starts at trimming.

## ✂️ `trim`

Two-pass primer trimming with Cutadapt.

```
seednap trim INPUT_DIR [OPTIONS]
```

| Option | Required | Description |
|---|---|---|
| `--forward-primer TEXT` | Yes | Forward primer sequence (5' to 3') |
| `--reverse-primer TEXT` | Yes | Reverse primer sequence (5' to 3') |
| `-o, --output-dir PATH` | Yes | Output directory for trimmed reads |
| `-c, --cores INTEGER` | No | Number of CPU cores (default: 1) |

```bash
seednap trim /path/to/raw/fastq \
  --forward-primer ACACCGCCCGTCACTCT \
  --reverse-primer CTTCCGGTACACTTACCATG \
  -o outputs/01_trim/teleo
```

In a full run this is the `trimming:` config section. See [configuration.md](configuration.md).

## 🐝 `swarm`

SWARM OTU clustering on trimmed reads.

```
seednap swarm MARKER TRIMMED_READS_DIR [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-o, --output-dir PATH` | `outputs/` | Base output directory |
| `-d, --distance INTEGER` | `1` | SWARM distance threshold |
| `-t, --threads INTEGER` | `4` | CPU threads |
| `--no-fastidious` | off | Disable singleton refinement |
| `--no-chimera-filter` | off | Skip de novo chimera detection (chimeras are artificial sequences formed when two parent templates fuse during PCR) |

```bash
seednap swarm teleo /path/to/trimmed -o outputs -d 1 -t 8
```

In a full run this is the `swarm:` config section. See [configuration.md](configuration.md).

## 🧬 `dada2`

DADA2 ASV processing on trimmed reads.

```
seednap dada2 MARKER TRIMMED_READS_DIR [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-o, --output-dir PATH` | `outputs/` | Base output directory |
| `--max-ee FLOAT` | `2.0` | Maximum expected errors for filtering |
| `--trunc-q INTEGER` | `11` | Truncate at first base with quality at or below this |
| `--min-overlap INTEGER` | `20` | Minimum overlap for merging paired reads |
| `--assign-taxonomy` | off | Run DADA2 taxonomic assignment (requires `--rdp-db` and `--species-db`) |
| `--rdp-db PATH` | | RDP-formatted taxonomy database |
| `--species-db PATH` | | Species-level database |
| `--library-map PATH` | | Per-library error learning from a `sample,library` CSV |

```bash
seednap dada2 teleo /path/to/trimmed -o outputs --max-ee 2.0 --trunc-q 11
```

`--library-map` learns a separate DADA2 error model per sequencing library and merges the per-library tables (2+ libraries denoised separately, then identical ASVs collapsed; with a single library, or when the flag is omitted, the standard single-batch path runs). In `run-pipeline` this is the `dada2.per_library` config field (default `false`), with libraries grouped from the manifest's `seq_run_id`. See [configuration.md](configuration.md).

The standalone `dada2` command always collects ASV metrics (writes `metrics.json`/`csv`) and does not expose the chimera method or pooling flags. `dada2.chimera.method` and `dada2.pool` are configurable only via `run-pipeline`.

## 🧼 `clean`

Decontaminate an abundance table against its negative controls. A negative control (or "blank") is a sample carried through extraction and/or PCR with no biological template; any sequence it contains is contamination. SeeDNAP derives control identity (extraction vs PCR blanks, and which extraction batch each belongs to) from a FAIRe sample manifest migrated from the field metadata.

Cleaning is presence-based and operates at the OTU/ASV level: any OTU/ASV that has reads in an applicable control is removed from (`subtract`) or annotated in (`flag`) the samples that control covers. An extraction blank cleans only the biological samples that share its `extraction_ID` batch; a PCR blank cleans the whole dataset.

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

In a full run this is the `cleaning:` config section (`mode`, default `flag`); it runs only when `clean` is listed in `pipeline.steps`, after `taxonomy` and before `export`. There it cleans the taxonomy-annotated table (not the standalone `ABUNDANCE_CSV` argument) and takes control identity from `report.sample_metadata`; if that is not set, cleaning is skipped with a `[WARN]` and export uses the uncleaned table. See [configuration.md](configuration.md).

## 🎯 `blast`

BLAST taxonomic assignment with LCA resolution. When a query sequence has several good database hits that disagree on identity, the **LCA** (lowest common ancestor) is the most specific rank on which those hits still agree, so the sequence is assigned only as deeply as the evidence supports.

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

**LCA algorithms.** `cascade` (default) applies the per-rank thresholds within a MEGAN-LR bitscore band (`--top-bitscore-pct`) above a `--lca-pident-delta` floor. `collapsed_taxonomy` is the eDNAFlow/OceanOmics %identity-window collapse-to-LCA: it ignores the per-rank thresholds and collapses disagreeing hits within `--lca-diff` %id of one another, above the `--lca-pid` floor. It reads the lineage from the reference FASTA headers, needs no NCBI taxdump, and runs offline. See [taxonomy-methods.md](taxonomy-methods.md#blast--lca-recommended) for details.

```bash
seednap blast outputs/02_swarm/teleo/query.fasta \
  /path/to/reference.fasta \
  outputs/02_swarm/teleo/otu_table.csv \
  -o outputs/teleo_blast.csv \
  --evalue 1e-10 --threshold-species 100
```

Reference headers may carry the literal string `NA` for an unknown rank (2025 CRABS DBs); SeeDNAP normalizes `NA`/empty/`nan` to a genuine missing rank, which surfaces as `Unassigned` in the export rather than a taxon named `NA`.

See [configuration.md](configuration.md#taxonomy).

## 🔬 `assign-taxonomy`

Generic taxonomic assignment supporting all four methods.

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

<details>
<summary><b>Additional options (thresholds, LCA, processors)</b></summary>

| Option | Default | Description |
|---|---|---|
| `--config PATH` | | Marker YAML config: reuse its `taxonomy.databases.<method>` block (database path plus method parameters) so the standalone result matches `run-pipeline`. Any explicit option below overrides the config value |
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
| `--confidence-threshold INTEGER` | `60` | Confidence threshold (DECIPHER) |
| `-c, --processors INTEGER` | `8` | CPU cores |

</details>

Unlike the standalone `blast` command, the `assign-taxonomy` BLAST path does NOT expose `--task`, `--perc-identity`, `--qcov-hsp-perc`, or `--evalue` as flags; left on their own they stay at their defaults (`task=megablast`, `perc_identity=80`, `qcov_hsp_perc=80`, `evalue=1e-25`). To change them, pass `--config <marker.yaml>` (which carries those `taxonomy.databases.blast` values), or use the `blast` command.

When you pass `--config`, the config's `taxonomy.method` must equal the positional METHOD argument. Asking for one method with a config written for another is rejected with an error rather than silently applying the wrong parameter set.

## 📊 `format-gbif`

Convert taxonomy results to GBIF long format.

```
seednap format-gbif INPUT_FILE [OPTIONS]
```

| Option | Required | Description |
|---|---|---|
| `-f, --format {dada2\|ecotag\|blast\|decipher}` | Yes | Input format type |
| `-o, --output PATH` | No | Output path (default: `<input>_gbif_input.csv`) |

```bash
seednap format-gbif outputs/teleo_blast.csv -f blast -o outputs/teleo_gbif.csv
```

## 🌍 `create-gbif`

Build a full DarwinCore-compliant GBIF occurrence CSV. DarwinCore is the GBIF standard vocabulary for biodiversity records (each row is an occurrence with columns like `eventID`, `decimalLatitude`, `eventDate`, `scientificName`); this is the file you submit to GBIF.

```
seednap create-gbif TAXONOMY_RESULTS SAMPLE_METADATA PROJECT_METADATA OUTPUT [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `--summarise-pcr / --no-summarise-pcr` | off | Aggregate PCR replicates per sample (strips the `_XX` replicate suffix from the eventID and sums reads) |
| `--skip-enrichment` | off | Skip NCBI/WORMS taxonomy enrichment (kingdom/phylum lookup) |

```bash
seednap create-gbif outputs/teleo_gbif.csv metadata/samples.csv metadata/project.csv outputs/teleo_darwincore.csv
```

Sample metadata is joined onto the occurrences by `eventID`, on a normalized key so that a taxonomy table whose eventIDs were rewritten by R's `make.names()` (dots for dashes) still matches the canonical dashed eventIDs in the metadata; the output carries the canonical form. If some eventIDs do not match, a `[WARN]` names them and they get blank location/date fields. If ZERO eventIDs match, the build raises rather than emit a file with every coordinate and date blank.

> [!WARNING]
> Enrichment needs `NCBI_API_KEY` in `.env` (see `.env.example`). Without a key, NCBI throttles the step heavily; pass `--skip-enrichment` to skip it entirely.

See [gbif-export.md](gbif-export.md) for the metadata column requirements.

## 🗄️ `wis-metadata`

Generate the GBIF export's sample and project metadata CSVs from the **WIS database** (the PostgreSQL/PostGIS schema built by `wis_database_creator`), instead of hand-writing them. Reads each sample's `eventID`, date, coordinates (from the PostGIS point), environmental medium and size, and writes `<marker>_sample_metadata.csv` + `<marker>_project_metadata.csv`.

This is an optional feature: install the database extra first with `pip install 'seednap[wis]'` (adds SQLAlchemy + psycopg2; intentionally not in the core install).

```
seednap wis-metadata --database-url URL --marker MARKER --output-dir DIR \
  --recorded-by NAME --identification-remarks TEXT --identification-references TEXT [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `--database-url` | `$WIS_DATABASE_URL` | SQLAlchemy URL, e.g. `postgresql://user:pass@host:5432/wis` (required) |
| `--marker` | (required) | Marker name; used for the project row and the output filenames |
| `--output-dir` | (required) | Directory for the two CSVs |
| `--monitoring` | none | Restrict to one WIS `monitoring_id` (site / long-term project) |
| `--mission` | none | Restrict to one WIS `mission_id` (sampling campaign) |
| `--event-id-field` | `sample_id` | WIS identifier used as `eventID` (`sample_id` or `material_sample_id`); match your FASTQ naming |
| `--recorded-by` | (required) | DwC `recordedBy` for the project row |
| `--identification-remarks` | (required) | Identification-method note for the project row |
| `--identification-references` | (required) | Reference-DB / method citation for the project row |
| `--seq-meth` | empty | Sequencing-method description |
| `--otu-seq-comp-appr` | empty | OTU/ASV sequence-comparison approach |

```bash
seednap wis-metadata --database-url postgresql://user:pass@host:5432/wis \
  --marker teleo --monitoring fw_ch_rechy --output-dir metadata/ \
  --recorded-by "ELE Lab" --identification-remarks "BLAST + LCA" \
  --identification-references "10.1038/nmeth.3869"
```

See [gbif-export.md](gbif-export.md#sourcing-metadata-from-the-wis-database) for the full WIS-to-DarwinCore field mapping and the `env_medium` / `eventID` notes.

## 🔀 `demultiplex`

Demultiplex ligation-based libraries (Cutadapt under the hood).

```
seednap demultiplex RAW_READS_DIR LIBRARY_NAME METADATA_CSV [OPTIONS]
```

`METADATA_CSV` must contain `eventID`, `tag_demultiplex`, and `library` columns.

| Option | Required | Description |
|---|---|---|
| `-f, --forward-primer TEXT` | Yes | Forward primer sequence |
| `-r, --reverse-primer TEXT` | Yes | Reverse primer sequence |
| `-o, --output-dir PATH` | Yes | Output base directory |
| `-c, --cores INTEGER` | No | CPU cores (default: 1) |
| `--no-gunzip` | No | Keep output files gzipped (default: outputs are gunzipped) |

> [!WARNING]
> In `run-pipeline`, listing `demultiplex` in `pipeline.steps` with any `demultiplex.protocol` other than `ligation` is REJECTED AT CONFIG LOAD. Only the `ligation` protocol is implemented.

## 📋 `manifest`

Build (and optionally validate) a canonical FAIRe sample manifest from the lab's existing CSVs. The manifest is the source of truth for sample -> library grouping (used by `dada2.per_library`) and control identity (used by `clean`). Standalone and read-only on its inputs.

```
seednap manifest FIELD_METADATA [OPTIONS]
```

`FIELD_METADATA` is the per-sample field metadata CSV (`metadata_field_*.csv`), or a legacy demux lab CSV (`metadata_lab_*.csv`) carrying library/tag columns.

| Option | Default | Description |
|---|---|---|
| `--project-metadata PATH` | | Project metadata CSV (marker -> target_gene/assay_name) |
| `--lab-metadata PATH` | | Legacy demux metadata CSV (seq_run_id/library and tag barcodes) |
| `--seq-run-id TEXT` | | Sequencing-run id for the whole dataset (overrides lab/derived value) |
| `--target-gene TEXT` | | Marker / target_gene (overrides the project metadata) |
| `--date-order {ymd\|dmy\|mdy}` | | Force eventDate field order for ambiguous dotted dates (otherwise such files raise) |
| `-o, --output PATH` | | Write the canonical manifest CSV here |
| `--abundance PATH` | | Validate the manifest's eventIDs against this abundance/OTU table |
| `--strict` | off | Raise if the abundance table has sample columns absent from the manifest (default: warn) |

```bash
seednap manifest metadata/metadata_field_my_dataset.csv \
  --project-metadata metadata/metadata_proj_my_dataset.csv \
  --abundance outputs/02_swarm/teleo/otu_table.csv \
  -o metadata/manifest_my_dataset.csv
```

## 🆕 `init`

Generate an example configuration file.

```
seednap init [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-m, --marker TEXT` | `teleo` | Marker name |
| `-o, --output PATH` | `config/markers/example.yaml` | Output path |
| `--minimal / --full` | `--minimal` | Required-fields-only config (default) or the fully-annotated reference template |
| `-f, --force` | off | Overwrite existing file |

`--minimal` (the default) emits ONLY the required fields, leaving everything else on built-in defaults. Pass `--full` for the fully-annotated template that shows every knob.

## ✅ `validate`

Validate a YAML configuration file. Checks syntax, field types, and required fields (including typos inside `taxonomy.databases.<method>`, which are rejected at load). The summary reports which database block is live for the selected method and runs a preflight on referenced files.

```
seednap validate CONFIG
```

`validate` runs a preflight and EXITS 1 if a referenced database or `raw_data` path is missing on disk, or the taxonomy DB block is unresolved. A syntactically-valid config can still fail here. This is the recommended pre-run gate before `run-pipeline`.

## 📈 `report`

Build the read/sequence-tracking report (and optionally the HTML run report) from an existing run's outputs. When `report` is in `pipeline.steps` (it is by default), `run-pipeline` generates both at the end of the run (the HTML document is gated by `report.html_report`, default `true`); this command REGENERATES them from outputs that already exist.

```
seednap report MARKER [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-o, --output-dir PATH` | `outputs/` | Base output directory of the run |
| `--html` | off | Also generate the self-contained HTML run report |
| `--warn-retention FLOAT` | `30.0` | Warn for samples retaining below this % of raw reads |
| `--warn-step-loss FLOAT` | `70.0` | Warn when a single step drops more than this % of a sample's reads |
| `--field-metadata PATH` | auto | Field metadata CSV (location, dates, sites) for the Dataset section |
| `--project-metadata PATH` | auto | Project metadata CSV (recorder, sequencing, reference DB) |
| `--log-file PATH` | auto | Run log to embed in the HTML report; auto-located under `logs/` if omitted |

```bash
seednap report teleo -o outputs --html \
  --field-metadata metadata/metadata_field_my_dataset.csv \
  --project-metadata metadata/metadata_proj_my_dataset.csv
```

Writes `outputs/04_report/<marker>/read_tracking.{csv,txt}` and, with `--html`, `report.html`. See [reporting.md](reporting.md).

## 🔍 `monitor`

Summarise a finished or in-progress run from its state JSON. Prints a per-step status/duration table plus the read-tracking headline (raw -> final reads, mean retention, warnings), and writes `monitoring_summary.csv` when per-sample counts are present. Standalone, read-only, and regenerable any time without a re-run.

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

## 💡 `explain`

Explain a SeeDNAP error code in depth. Error messages may reference a code (e.g. `SDN-CFG-001`); pass it here for the full what/why/fix detail.

```
seednap explain [CODE]
```

With no `CODE`, lists every known error code with its title.

```bash
seednap explain SDN-CFG-001
seednap explain
```

## 🏷️ `version`

Show detailed version information.

```
seednap version
```

## 📖 See also

| Guide | What's inside |
|---|---|
| [configuration.md](configuration.md) | Every config key with type, default, and meaning |
| [pipeline-steps.md](pipeline-steps.md) | Per-stage behavior and algorithms |
| [taxonomy-methods.md](taxonomy-methods.md) | BLAST/DADA2/DECIPHER/ecotag details |
| [gbif-export.md](gbif-export.md) | DarwinCore export and metadata columns |
| [reporting.md](reporting.md) | Read tracking and the HTML run report |
</content>
</invoke>
