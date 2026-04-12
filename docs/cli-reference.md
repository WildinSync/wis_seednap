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
| `--threshold-species FLOAT` | `98.0` | Percent identity for species assignment |
| `--threshold-genus FLOAT` | `96.0` | Percent identity for genus assignment |
| `--threshold-family FLOAT` | `86.5` | Percent identity for family assignment |

```bash
seednap blast outputs/02_swarm/teleo/query.fasta \
  /path/to/reference.fasta \
  outputs/02_swarm/teleo/otu_table.csv \
  -o outputs/teleo_blast.csv \
  --evalue 1e-10 --threshold-species 100
```

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
| `--threshold-species FLOAT` | `98.0` | Species %ID threshold (BLAST) |
| `--threshold-genus FLOAT` | `96.0` | Genus %ID threshold (BLAST) |
| `--threshold-family FLOAT` | `86.5` | Family %ID threshold (BLAST) |
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

## `init`

Generate an example configuration file.

```
seednap init [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-m, --marker TEXT` | `teleo` | Marker name |
| `-o, --output PATH` | `config/markers/example.yaml` | Output path |
| `-f, --force` | | Overwrite existing file |

---

## `validate`

Validate a YAML configuration file. Checks syntax, types, required fields.

```
seednap validate CONFIG
```

---

## `version`

Show detailed version information.

```
seednap version
```
