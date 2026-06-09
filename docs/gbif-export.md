# GBIF and DarwinCore Export

How to turn a taxonomy table into a DarwinCore-compliant occurrence CSV for GBIF publishing.

Export is a two-stage process: reshape the wide taxonomy table into long format
(`format-gbif`), then merge it with sample and project metadata into a full
DarwinCore occurrence table (`create-gbif`). The first stage also runs
automatically as the pipeline `export` step.

> [!NOTE]
> The pipeline runs the same long-format conversion automatically when `export`
> is in `pipeline.steps`. It writes `<paths.output>/<marker>_<method>_gbif.csv`
> and honours the `export.gbif.add_rank` / `export.gbif.add_taxon` config keys
> (both default `true`). The `format-gbif` and `create-gbif` commands are the
> manual equivalents for working from existing files. See
> [configuration.md](configuration.md) for the `export` block and
> [pipeline-steps.md](pipeline-steps.md) for the step model.

## Step 1: Format for GBIF (`format-gbif`)

Converts the wide taxonomy table (OTUs as rows, samples as columns) into GBIF
long format (one row per sample-OTU observation).

```bash
seednap format-gbif outputs/teleo_blast.csv -f blast -o outputs/teleo_gbif.csv
```

| Flag | Type | Default | Meaning |
|---|---|---|---|
| `INPUT_FILE` (arg) | path | required | Wide taxonomy CSV from the taxonomy step |
| `-f` / `--format` | choice | required | Input parser: `dada2`, `ecotag`, `blast`, or `decipher` |
| `-o` / `--output` | path | `<input_stem>_gbif_input.csv` | Output path |

> [!NOTE]
> `blast` and `decipher` are parsed identically to `dada2` (same wide-table
> schema). `ecotag` differs: it renames `*_name` columns and drops ecotag
> metadata columns before reshaping.

### What it does

1. Reshapes wide format to long (one row per sample-OTU pair), excluding the
   per-OTU annotation columns `ASV_ID`, `pident`, and `is_contaminant_candidate`
   from the sample set.
2. Drops zero-count observations.
3. Adds a `rank` column (`species`, `genus`, `family`, or `higher`) when
   `add_rank` is set.
4. Adds a `taxon` column (lowest available taxonomic name) when `add_taxon` is set.
5. Normalises upstream schema differences: a capital-S `Sequence` column is
   renamed to `sequence`, and the literal `Unassigned` taxonomy value is mapped
   to empty before rank/taxon are computed.

> [!NOTE]
> On the manual `format-gbif` command, `rank` and `taxon` are always added. On
> the pipeline `export` step they are controlled by `export.gbif.add_rank` and
> `export.gbif.add_taxon` (both default `true`).

### Output columns

`kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `taxon`,
`rank`, `sequence`, `nb_reads`, `eventID`.

The `is_contaminant_candidate` column is also appended when the upstream
taxonomy table carried it (that is, when `taxonomy.contaminants` was set).
`create-gbif` reads this column to populate `contamination_flag`.

## Step 2: DarwinCore Publishing (`create-gbif`)

Merges the long-format taxonomy table with sample and project metadata to
produce a full DarwinCore occurrence CSV.

```bash
seednap create-gbif taxonomy_gbif.csv sample_metadata.csv project_metadata.csv output.csv
```

| Argument / flag | Type | Default | Meaning |
|---|---|---|---|
| `TAXONOMY_RESULTS` (arg) | path | required | Long-format output from `format-gbif` |
| `SAMPLE_METADATA` (arg) | path | required | Per-sample metadata CSV |
| `PROJECT_METADATA` (arg) | path | required | Single-row project metadata CSV |
| `OUTPUT` (arg) | path | required | Destination DarwinCore CSV |
| `--summarise-pcr` | flag | `false` | Collapse PCR replicates and sum reads |
| `--skip-enrichment` | flag | `false` | Skip NCBI/WORMS kingdom/phylum lookup |

### What it does

1. Validates inputs up front: coordinate ranges, `env_medium` values, date
   formats, and required metadata columns (see below). Invalid input raises a
   clear error instead of writing a corrupt submission.
2. Removes control samples by name (see the control-removal warning below).
3. Optionally summarises PCR replicates (`--summarise-pcr`).
4. Filters non-target taxa for the marker, then sums reads per sample.
5. Looks up primer and `target_gene` details from the bundled `primers_list.csv`.
6. Assigns a deterministic `occurrenceID`, maps `env_medium` to ENVO terms, and
   merges sample metadata (coordinates, dates, depth).
7. Enriches missing `kingdom`/`phylum` via NCBI and WORMS unless skipped.
8. Propagates the upstream `is_contaminant_candidate` flag into a
   `contamination_flag` boolean. Rows are never dropped on this flag; it is
   informational and downstream decides what to do.
9. Validates that every required DarwinCore output field is populated, then
   writes the CSV.

The output has 39 columns, including `scientificName` (the lowest assigned
name), six taxonomic ranks (`kingdom` through `genus`), occurrence/event
fields, location fields, sequencing fields, and `contamination_flag`.

> [!WARNING]
> Looking up an unknown marker hard-fails. If `project_metadata.marker` is not
> in the bundled `primers_list.csv`, the build aborts rather than emitting blank
> `target_gene` and primer columns.

> [!NOTE]
> The `occurrenceID` is `marker:eventID:sha256(sequence)[:8]` (literal `NOSEQ`
> in place of the hash when a row has no sequence). It is deterministic across
> re-runs of the same data, so resubmitting a dataset replaces records in GBIF
> rather than duplicating them.

### Options

| Option | Meaning |
|---|---|
| `--summarise-pcr` | Group replicates of the same sample and sum their reads. Replicates are detected by a trailing two-digit suffix on `eventID` matching `_NN` (for example `S1_01`, `S1_02` collapse to `S1`). |
| `--skip-enrichment` | Skip the NCBI/WORMS API calls; `kingdom`/`phylum` stay as supplied by upstream taxonomy. |

### Taxonomy enrichment and the NCBI API key

Enrichment fills missing `kingdom`/`phylum`/`class` by querying NCBI Entrez and
WORMS. Provide an NCBI API key in a `.env` file at the project root:

```text
NCBI_API_KEY=your_key_here
```

Get a key at https://www.ncbi.nlm.nih.gov/account/settings/.

> [!IMPORTANT]
> Without an NCBI API key the enrichment step is skipped entirely (it logs a
> `[WARN]` and returns the table unchanged), so `kingdom`/`phylum` may remain
> empty. Passing `--skip-enrichment` skips it regardless of the key. Enrichment
> only fills empty cells; existing upstream taxonomy is never overwritten.

### Input: Sample Metadata CSV

One row per sample. `eventID`, `eventDate`, and `env_medium` are required; the
rest are optional and validated only when present.

| Column | Required | Meaning | Example |
|---|---|---|---|
| `eventID` | yes | Sample identifier; must match the per-sample columns in the taxonomy table | `SPY221633_01` |
| `eventDate` | yes | Collection date, `yyyy`, `yyyy.mm`, or `yyyy.mm.dd` | `2023.06.15` |
| `decimalLatitude` | no | Latitude, validated to `[-90, 90]` | `-22.5678` |
| `decimalLongitude` | no | Longitude, validated to `[-180, 180]` | `17.1234` |
| `env_medium` | yes | Environment type, validated against the ENVO whitelist | `water` |
| `samp_size` | no | Sample volume/size (legacy `volume` is auto-renamed) | `1L` |
| `depth` | no | Depth in meters | `0.5` |
| `size_frac` | no | Filter size fraction | `0.22um` |

> [!WARNING]
> `env_medium` accepts only `water`, `soil`, `river`, `marine`, and `sediment`
> (case-insensitive). Any other value raises and aborts the build, both at input
> validation and at ENVO mapping. This prevents silently mislabelling samples.

> [!WARNING]
> `eventDate` month and day must be zero-padded to two digits. `2023.06.05` is
> valid; `2023.6.5` is rejected.

> [!NOTE]
> Legacy column names are auto-renamed on input: `filter_code` to `eventID` in
> the taxonomy table, and `volume` to `samp_size` in sample metadata.

### Input: Project Metadata CSV

A single-row table describing the run. Required columns are `marker`,
`recordedby`, `identificationRemarks`, and `identificationReferences`; the rest
are optional at load time.

| Column | Required | Meaning | Example |
|---|---|---|---|
| `marker` | yes | Marker name; must exist in `primers_list.csv` | `teleo` |
| `recordedby` | yes | Data recorder / contributor | `J. Smith` |
| `identificationRemarks` | yes | Method description | `BLAST + LCA` |
| `identificationReferences` | yes | Reference DOIs | `10.1038/nmeth.3869` |
| `seqmet` | no | Sequencing method | `MiSeq PE 2x150` |
| `otu_seq_comp_appr` | no | Sequence comparison approach | `SWARM d=1` |
| `otu_db` | no (load) | Reference database name | `CRABS MitoFish 2025` |
| `chimera_check` | no | Chimera method | `UCHIME de novo` |

> [!IMPORTANT]
> `otu_db` is optional at load time but is a required DarwinCore output field.
> If you leave it blank, the final required-field check fails and the write is
> aborted. Always supply an `otu_db` value (for example the reference database
> name) when building a real GBIF submission.

> [!WARNING]
> Headers are case-sensitive (for example `identificationReferences`, not
> `IdentificationReferences`), and the table must have exactly one data row.

> [!WARNING]
> Control removal in `create-gbif` is a name-based filter that drops only
> `eventID`s matching `blank`, `CNEG`, `CMET`, or `CEXT` (case-insensitive). It
> is independent of the manifest control classifier (`manifest.classify_control`)
> and recognises fewer patterns: it does not catch `CPCR`, `water`, or
> `EXT_NC`/`PCR_NC` controls. Such controls pass through into the GBIF output, so
> verify your control naming matches one of the four recognised patterns.

## See also

- [configuration.md](configuration.md) for the `export` block and `pipeline.steps`.
- [pipeline-steps.md](pipeline-steps.md) for where the `export` step runs.
- [taxonomy-methods.md](taxonomy-methods.md) for the taxonomy tables that feed export.
