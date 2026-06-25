# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- WIS database bridge (`seednap wis-metadata`): generate the GBIF export's
  per-sample and project metadata CSVs straight from the WIS PostgreSQL/PostGIS
  database (the schema built by `wis_database_creator`) instead of hand-writing
  them. Reads each sample's `eventID`, date, coordinates (from the PostGIS
  point), environmental medium (mapped from the controlled `sample_type` code to
  the builder's ENVO vocabulary; an unmapped medium passes through with a `[WARN]`
  rather than being mislabelled) and size, and writes the two CSVs the DarwinCore
  export already consumes. The DarwinCore builder is unchanged. SQLAlchemy and a
  PostgreSQL driver are an optional dependency (`pip install 'seednap[wis]'`);
  the core pipeline stays dependency-light and the bridge fails with a clear
  install hint if they are absent.
- Error-explainability module with a `seednap explain` command: errors carry
  stable codes and actionable what / why / how-to-fix detail, and the codes can
  be looked up from the CLI.
- Sample discovery now finds per-sample FASTQs inside per-library / per-run
  subdirectories of `paths.raw_data`, not only at the top level, so
  already-demultiplexed data organized one folder per library is processed
  without flattening. Sample names must be unique across subfolders (an
  ambiguous name is rejected rather than guessed).
- DADA2-by-library (`dada2.per_library`) now derives the sample-to-library
  grouping from the per-library subfolders of `raw_data` when no metadata is
  configured, instead of silently falling back to a single-batch run.
- Early heavy-read-loss warning at the trim step: if primer trimming discards
  more than `report.warn_step_loss_pct` of the reads (summed across the run), a
  loud `[WARN]` is logged immediately, before the long downstream steps, naming
  the likely cause and fix. The textbook case is feeding already-primer-trimmed
  FASTQs into the default `trimming.discard_untrimmed: true` path (set it to
  `false`); a genuine low yield (off-target amplification, primer mismatch) is
  flagged too so the warning is not misread.
- New `darwincore` pipeline step: builds the GBIF-ready DarwinCore occurrence
  file in-pipeline (joining the long-format export to `report.sample_metadata` +
  `report.project_metadata`, with `export.darwincore` flags), rather than only
  via the standalone `create-gbif` command. Opt-in via `pipeline.steps`; required
  metadata is checked at config preflight. The reference-database (`otu_db`) and
  chimera-removal (`chimera_check`) provenance are filled automatically from the
  run config (a differing project-metadata value is warned and overridden). It also
  writes a deleted-entries report (`<output>_dropped.csv`) of the occurrences removed
  by the control and non-target filters, for QA.

### Changed

- Pipeline stage enable/disable now flows through a single `pipeline.steps`
  config model with dependency validation, replacing the previous scattered
  per-stage toggles.
- Standardized the documentation for accuracy and a consistent structure so the
  docs match the implementation.
- The shipped reference marker configs now use placeholder `raw_data` / metadata
  paths, so an unedited run fails the config preflight instead of silently
  processing a bundled example dataset.

### Fixed

- Correctness sweep across the pipeline focused on data integrity, removing
  silent fallbacks (fallbacks now warn or fail loudly), and catching
  wrong-environment misconfiguration earlier.
- The run report now reads the per-sample Cutadapt logs from the trim step's
  output directory, so raw and trimmed read counts (and `% retained`) are
  reported instead of `NA`.
- The trim step clears stale outputs from a previous run before writing, so a
  re-run that finds a different sample set cannot reuse the earlier run's
  trimmed reads downstream.
- Sample discovery fails loudly when `paths.raw_data` holds no FASTQ files,
  instead of returning an empty list and producing an empty run.
- Standalone CLI commands: `trim` no longer crashes on a shallow output path;
  `clean` no longer treats per-OTU annotation columns as biological samples;
  `create-gbif` normalizes eventID separators and warns or fails on a
  zero-match metadata join (rather than silently emitting blank dates and
  coordinates); `assign-taxonomy --config` uses the marker config's BLAST
  parameters; `format-gbif` logs the actual input format.

## [0.1.0]

- Initial alpha release of the SeeDNAP eDNA metabarcoding pipeline.
