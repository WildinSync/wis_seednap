# GBIF and DarwinCore Export

<img src="../media/divider.svg" width="100%" alt="">

How to turn a taxonomy table into a DarwinCore-compliant occurrence CSV for GBIF publishing.

GBIF (the Global Biodiversity Information Facility) is the public repository the lab submits occurrence records to. DarwinCore is the standardised biodiversity data vocabulary GBIF ingests: a fixed set of column names (`eventID`, `scientificName`, `decimalLatitude`, and so on) that make records comparable across datasets. The two commands here convert one of the pipeline's taxonomy tables into that vocabulary.

Export is a two-stage process: reshape the wide taxonomy table into long format (`format-gbif`), then merge it with sample and project metadata into a full DarwinCore occurrence table (`create-gbif`). The first stage also runs automatically as the pipeline `export` step.

<details>
<summary><b>Running export automatically as the pipeline <code>export</code> step</b></summary>

The pipeline runs the same long-format conversion automatically when `export` is in `pipeline.steps`. It writes `<paths.output>/<marker>_<method>_gbif.csv` and honours the `export.gbif.add_rank` / `export.gbif.add_taxon` config keys (both default `true`). If a `clean` step ran before `export` (decontamination against the blank/negative-control samples) and produced a cleaned table, the export step uses that decontaminated table instead of the raw taxonomy table. The `format-gbif` and `create-gbif` commands are the manual equivalents for working from existing files. See [configuration.md](configuration.md) for the `export` block and [pipeline-steps.md](pipeline-steps.md) for the step model.

</details>

## 🔄 Step 1: Format for GBIF (`format-gbif`)

Converts the wide taxonomy table into GBIF long format. A wide taxonomy table has one row per OTU or ASV (an OTU is a cluster of similar sequences; an ASV is a single denoised sequence variant; both stand in for a taxon) and one numeric column per sample holding that sequence's read count. Long format has one row per sample-OTU observation, with a single `nb_reads` count.

```bash
seednap format-gbif outputs/teleo_blast.csv -f blast -o outputs/teleo_gbif.csv
```

| Flag | Type | Default | Meaning |
|---|---|---|---|
| `INPUT_FILE` (arg) | path | required | Wide taxonomy CSV from the taxonomy step |
| `-f` / `--format` | choice | required | Input parser: `dada2`, `ecotag`, `blast`, or `decipher` |
| `-o` / `--output` | path | `<input_stem>_gbif_input.csv` | Output path |

`blast` and `decipher` are parsed identically to `dada2` (same wide-table schema). `ecotag` differs: it renames `*_name` columns and drops ecotag metadata columns before reshaping.

### What it does

1. Normalises upstream schema differences first: a capital-S `Sequence` column is renamed to `sequence`, and the literal `Unassigned` taxonomy value is mapped to empty so it is treated as a gap, not a real name.
2. Reshapes wide format to long (one row per sample-OTU pair), excluding the per-OTU annotation columns `ASV_ID`, `pident`, and `is_contaminant_candidate` from the sample set.
3. Drops zero-count observations.
4. Adds a `rank` column (`species`, `genus`, `family`, or `higher`) when `add_rank` is set: the finest rank that is confidently assigned. A species name containing `/` is an ambiguous tie between species and is treated as resolved only to genus.
5. Adds a `taxon` column (the lowest available taxonomic name) when `add_taxon` is set; this becomes `scientificName` in the final DarwinCore output.

On the manual `format-gbif` command, `rank` and `taxon` are always added. On the pipeline `export` step they are controlled by `export.gbif.add_rank` and `export.gbif.add_taxon` (both default `true`).

### Output columns

`kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `taxon`, `rank`, `sequence`, `nb_reads`, `eventID`.

The `is_contaminant_candidate` column is also appended when the upstream taxonomy table carried it (that is, when `taxonomy.contaminants` was set). `create-gbif` reads this column to populate `contamination_flag`.

## 🌍 Step 2: DarwinCore Publishing (`create-gbif`)

Merges the long-format taxonomy table with sample and project metadata to produce a full DarwinCore occurrence CSV.

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
| `--skip-enrichment` | flag | `false` | Skip NCBI/WORMS higher-taxonomy lookup |

### What it does

1. Validates inputs up front: coordinate ranges, `env_medium` values, date formats, and required metadata columns (see below). Invalid input raises a clear error instead of writing a corrupt submission.
2. Optionally summarises PCR replicates (`--summarise-pcr`).
3. Removes control samples by name (see [Control removal](#control-removal) below).
4. Looks up primer and `target_gene` details from the bundled `primers_list.csv`, then filters non-target taxa for the marker and sums reads per sample.
5. Maps `env_medium` to ENVO terms, merges sample metadata (coordinates, dates, depth) on a normalised `eventID` key (see [eventID matching](#eventid-matching) below), and assigns a deterministic `occurrenceID`.
6. Enriches missing `kingdom`/`phylum`/`class` via NCBI and WORMS unless skipped.
7. Propagates the upstream `is_contaminant_candidate` flag into a `contamination_flag` boolean. Rows are never dropped on this flag; it is informational and downstream decides what to do.
8. Validates that every required DarwinCore output field is populated, then writes the CSV.

The output has 39 columns, including `scientificName` (the lowest assigned name), six taxonomic ranks (`kingdom` through `genus`), occurrence/event fields, location fields, sequencing fields, and `contamination_flag`.

> [!WARNING]
> Looking up an unknown marker hard-fails. If `project_metadata.marker` is not in the bundled `primers_list.csv`, the build aborts rather than emitting blank `target_gene` and primer columns.

### Options

| Option | Meaning |
|---|---|
| `--summarise-pcr` | Group replicates of the same sample and sum their reads. Replicates are detected by a trailing two-digit suffix on `eventID` matching `_NN` (for example `S1_01`, `S1_02` collapse to `S1`). |
| `--skip-enrichment` | Skip the NCBI/WORMS API calls. `class` is carried from the upstream taxonomy, but `kingdom` and `phylum` are populated only by enrichment, so they stay blank when it is skipped. |

<details>
<summary><b><code>occurrenceID</code> construction</b></summary>

The `occurrenceID` is `marker:eventID:sha256(sequence)[:8]` (literal `NOSEQ` in place of the hash when a row has no sequence). The 8-character hash is the first eight hex digits of the SHA-256 digest of the uppercased sequence. The ID is deterministic across re-runs of the same data, so resubmitting a dataset replaces records in GBIF rather than duplicating them.

</details>

<details>
<summary><b><a name="eventid-matching"></a>eventID matching (normalisation and join behaviour)</b></summary>

Sample metadata is joined to the taxonomy table on a normalised `eventID`: any run of `.`, `_`, or `-` is collapsed to a single `-` and the result is upper-cased before matching. This is because the R taxonomy step rewrites a dashed identifier such as `DAR-2023-0025` to the dotted `DAR.2023.0025` when it becomes a column name, while the metadata sheet keeps the dashed form; the two forms must still join. The matched (canonical, dashed) `eventID` is what is written to GBIF.

If the normalised keys still collide (two metadata rows differing only by separators), the build aborts so reads are not fanned out. If the join matches ZERO rows the build aborts (every date and coordinate would be blank); if it matches some but not all rows, a `[WARN]` lists the unmatched `eventID`s and those rows ship with blank location/date fields.

</details>

### 🔑 Taxonomy enrichment and the NCBI API key

Taxonomic assignment usually resolves a read to a low rank (species or genus) but leaves the higher ranks blank, because a marker's reference database does not store the full lineage. GBIF expects those higher ranks. Enrichment fills the missing `kingdom`/`phylum`/`class` cells by looking each name up in NCBI Taxonomy (via the Entrez API) first and, when that returns nothing, in WORMS (the World Register of Marine Species). Provide an NCBI API key in a `.env` file at the project root:

```text
NCBI_API_KEY=your_key_here
```

Get a key at https://www.ncbi.nlm.nih.gov/account/settings/.

Without an NCBI API key the enrichment step is skipped entirely (it logs a `[WARN]` and returns the table unchanged), so the higher ranks may remain empty. Passing `--skip-enrichment` skips it regardless of the key. In the output the `class` rank is carried from the upstream taxonomy, but `kingdom` and `phylum` are filled only by enrichment, so when it is skipped (no key or `--skip-enrichment`) they remain blank. A name that neither NCBI nor WORMS resolves keeps its blank higher ranks, and the count of unresolved names is reported with a `[WARN]`.

### 📄 Input: Sample Metadata CSV

One row per sample. `eventID`, `eventDate`, and `env_medium` are required; the rest are optional and validated only when present.

| Column | Required | Meaning | Example |
|---|---|---|---|
| `eventID` | yes | Sample identifier; must match the per-sample columns in the taxonomy table | `SPY221633_01` |
| `eventDate` | yes | Collection date, `yyyy`, `yyyy.mm`, or `yyyy.mm.dd` | `2023.06.15` |
| `decimalLatitude` | no | Latitude, validated to `[-90, 90]` | `-22.5678` |
| `decimalLongitude` | no | Longitude, validated to `[-180, 180]` | `17.1234` |
| `env_medium` | yes | Environment type; must be one of the five values that map to an ENVO term (see warning) | `water` |
| `samp_size` | no | Sample volume/size (legacy `volume` is auto-renamed) | `1L` |
| `depth` | no | Depth in meters | `0.5` |
| `size_frac` | no | Filter size fraction | `0.22um` |

> [!WARNING]
> `env_medium` accepts only `water`, `soil`, `river`, `marine`, and `sediment` (case-insensitive). Each maps to a standard term from ENVO, the Environment Ontology that GBIF uses to describe sample media (for example `water` and `river` both become `liquid water [ENVO_00002006]`, `marine` becomes `sea water [ENVO_00002149]`). Any other value raises and aborts the build, both at input validation and at the ENVO mapping step. This prevents silently mislabelling samples.

> [!WARNING]
> `eventDate` month and day must be zero-padded to two digits. `2023.06.05` is valid; `2023.6.5` is rejected.

Legacy column names are auto-renamed on input: `filter_code` to `eventID` in the taxonomy table, and `volume` to `samp_size` in sample metadata.

### 📄 Input: Project Metadata CSV

A single-row table describing the run. Required columns are `marker`, `recordedby`, `identificationRemarks`, and `identificationReferences`; the rest are optional at load time.

| Column | Required | Meaning | Example |
|---|---|---|---|
| `marker` | yes | Marker name; must exist in `primers_list.csv` | `teleo` |
| `recordedby` | yes | Data recorder / contributor | `J. Smith` |
| `identificationRemarks` | yes | Free-text method description (for example BLAST plus LCA, lowest-common-ancestor assignment over the top hits) | `BLAST + LCA` |
| `identificationReferences` | yes | Reference DOIs | `10.1038/nmeth.3869` |
| `seqmet` | no | Sequencing method | `MiSeq PE 2x150` |
| `otu_seq_comp_appr` | no | Sequence comparison / clustering approach | `SWARM d=1` |
| `otu_db` | no (load) | Reference database name | `CRABS MitoFish 2025` |
| `chimera_check` | no | Chimera-detection method (a chimera is an artefactual sequence joining two parent templates during PCR) | `UCHIME de novo` |

> [!WARNING]
> `otu_db` is optional at load time but is a required DarwinCore output field. If you leave it blank, the final required-field check fails and the write is aborted. Always supply an `otu_db` value (for example the reference database name) when building a real GBIF submission.

> [!WARNING]
> Headers are case-sensitive (for example `identificationReferences`, not `IdentificationReferences`), and the table must have exactly one data row.

<details>
<summary><b><a name="control-removal"></a>Control removal</b></summary>

Controls (blanks and negative/positive controls: wells with no biological sample, used to detect contamination) must never be published as real biodiversity records, so `create-gbif` drops their rows. It uses `manifest.classify_control` as the single source of truth, so it removes the full set of recognised controls (`blank`, `CNEG`, `CMET`, `CEXT`, `CPCR`, `water`, `EXT_NC`/`PCR_NC`, and mock/positive forms), not just a few patterns. An `eventID` that looks control-like but cannot be classified is kept and reported with a `[WARN]` rather than silently dropped.

</details>

## 🗄️ Sourcing metadata from the WIS database

The two metadata CSVs above don't have to be hand-written. When your samples are registered in the **WIS database** (the normalized PostgreSQL/PostGIS schema built by `wis_database_creator`), `seednap wis-metadata` reads the per-sample field metadata straight from the database and writes the two CSVs the export consumes:

```bash
pip install 'seednap[wis]'   # one-time: adds SQLAlchemy + psycopg2 (optional, not in the core install)

seednap wis-metadata \
  --database-url postgresql://user:pass@host:5432/wis \
  --marker teleo \
  --monitoring fw_ch_rechy \
  --recorded-by "ELE Lab" \
  --identification-remarks "BLAST + LCA" \
  --identification-references "10.1038/nmeth.3869" \
  --output-dir metadata/
```

This writes `metadata/teleo_sample_metadata.csv` and `metadata/teleo_project_metadata.csv`; point `report.sample_metadata` / `report.project_metadata` (or the `create-gbif` arguments) at them. Filter to one site with `--monitoring <monitoring_id>` or one campaign with `--mission <mission_id>`.

What it maps, from the WIS schema to the export contract:

| Export field | WIS source |
|---|---|
| `eventID` | `sample_metadata.sample_id` (default) or `material_sample_id` (`--event-id-field`) |
| `eventDate` | `sample_metadata.event_date`, formatted `yyyy.mm.dd` |
| `decimalLatitude` / `decimalLongitude` | `ST_Y` / `ST_X` of the sample's `COORDINATE` point in `gis_point.geom` (SRID 4326) |
| `depth` / `size_frac` / `samp_size` | `sample_depth` / `sample_size_frac` / `sample_size` |
| `env_medium` | mapped from the controlled `sample_type` code (see below) |
| project row | `--marker` / `--recorded-by` / `--identification-*` (the WIS schema has no equivalent of these) |

> [!WARNING]
> **`env_medium` mapping.** The WIS database stores the environmental medium as a controlled `sample_type` code, not as an ENVO term. The bridge maps the aquatic / soil / sediment codes to the vocabulary the DarwinCore builder recognises: `FW`/`SU` → `water`, `MA` → `marine`, `SE` → `sediment`, `SO` → `soil` (the builder then maps those to ENVO). A code with no such analogue (`AI` air, `BL` blood, `HN` honey, `DW` deadwood) is passed through as its raw label with a `[WARN]`, so the builder fails loudly rather than mislabelling a published record. Extend `WIS_SAMPLE_TYPE_TO_ENV_MEDIUM` (and the builder's `_ENVO_TERMS`) to publish a new medium.

<details>
<summary><b><code>eventID</code> must match your sample names</b></summary>

The export joins on `eventID`, so the WIS identifier you choose has to correspond to the per-sample names in the taxonomy table (derived from the FASTQ filenames). The default is the short `sample_id`; switch to `material_sample_id` with `--event-id-field` if that is what your FASTQs are named after. The reference database (`otu_db`) and chimera-removal provenance are filled by the `darwincore` pipeline step from the run config, so they are not written by this command.

</details>

## 📖 See also

- [configuration.md](configuration.md) for the `export` block and `pipeline.steps`.
- [pipeline-steps.md](pipeline-steps.md) for where the `export` step runs.
- [taxonomy-methods.md](taxonomy-methods.md) for the taxonomy tables that feed export.
