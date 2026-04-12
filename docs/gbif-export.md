# GBIF and DarwinCore Export

SeeDNAP provides two export steps for biodiversity data publishing.

## Step 1: Format for GBIF (`format-gbif`)

Converts the wide-format taxonomy table (OTUs as rows, samples as columns) into GBIF long format (one row per observation).

```bash
seednap format-gbif outputs/teleo_blast.csv outputs/teleo_gbif.csv --method blast
```

### What it does

1. Transforms wide format to long format (one row per sample-OTU pair)
2. Removes zero-count observations
3. Adds `rank` column (species, genus, family, or higher)
4. Adds `taxon` column (lowest available taxonomic name)

### Output columns

`kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `taxon`, `rank`, `sequence`, `nb_reads`, `eventID`

---

## Step 2: DarwinCore Publishing (`create-gbif`)

Merges taxonomy results with sample and project metadata to produce a full DarwinCore-compliant occurrence CSV.

```bash
seednap create-gbif taxonomy_gbif.csv sample_metadata.csv project_metadata.csv output.csv
```

### What it does

1. Loads taxonomy results from `format-gbif` output
2. Removes control samples (blanks, negative controls: CNEG, CMET, CEXT)
3. Optionally summarises PCR replicates (aggregates reads per sample)
4. Validates date formats (yyyy, yyyy.mm, or yyyy.mm.dd)
5. Filters non-target taxa based on marker type
6. Looks up primer/marker details from the bundled primer list
7. Computes total reads per sample and generates occurrence IDs
8. Maps environment medium to ENVO ontology terms
9. Merges sample metadata (coordinates, dates, depth)
10. Enriches missing kingdom/phylum via NCBI Entrez and WORMS APIs
11. Populates all DarwinCore columns and exports CSV

### Options

| Option | Description |
|---|---|
| `--summarise-pcr` | Aggregate PCR replicates by sample |
| `--skip-enrichment` | Skip NCBI/WORMS API calls |

### NCBI API Key

Taxonomy enrichment requires an NCBI API key. Create a `.env` file at the project root:

```
NCBI_API_KEY=your_key_here
```

Get your key at https://www.ncbi.nlm.nih.gov/account/settings/

Without a key, enrichment is skipped and kingdom/phylum columns may remain empty.

### Input: Sample Metadata CSV

One row per sample:

| Column | Description | Example |
|---|---|---|
| `eventID` | Sample identifier (must match taxonomy) | `SPY221633_01` |
| `decimalLatitude` | Latitude | `-22.5678` |
| `decimalLongitude` | Longitude | `17.1234` |
| `eventDate` | Date (yyyy.mm.dd) | `2023.06.15` |
| `env_medium` | Environment type | `water` |
| `samp_size` | Sample volume/size | `1L` |
| `depth` | Depth in meters | `0.5` |
| `size_frac` | Filter size fraction | `0.22um` |

### Input: Project Metadata CSV

One row per project:

| Column | Description | Example |
|---|---|---|
| `marker` | Marker name | `teleo` |
| `recordedby` | Data recorder | `J. Smith` |
| `seqmet` | Sequencing method | `MiSeq PE 2x150` |
| `identificationRemarks` | Method description | `BLAST + LCA` |
| `identificationReferences` | Reference DOIs | `10.1038/nmeth.3869` |
| `otu_seq_comp_appr` | Sequence comparison | `SWARM d=1` |
| `otu_db` | Reference database | `CRABS MitoFish 2025` |
| `chimera_check` | Chimera method | `UCHIME de novo` |
