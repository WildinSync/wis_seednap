# seednap

Modern eDNA metabarcoding pipeline with DADA2.

A Python-first pipeline for processing eDNA metabarcoding data with support for multiple taxonomic assignment methods (DADA2, DECIPHER, Ecotag, BLAST).

**Version:** 0.1.0 (Phase 0 - Infrastructure Complete)

---

## Installation

### Requirements

**Python:** >= 3.9

**System dependencies:**
- `cutadapt` (>= 4.0)
- `R` (>= 4.0)
- R packages: `tidyverse`, `dada2`, `Biostrings`

**Optional (depending on taxonomic assignment method):**
- **DECIPHER:** R package `decipher`
- **BLAST:** `ncbi-blast` (makeblastdb, blastn)
- **Ecotag:** `obitools` toolkit (v1)

### Install from Source

1. Clone the repository:
```bash
git clone https://github.com/eth-edna/seednap.git
cd seednap
```

2. Install the package in development mode:
```bash
pip install -e .
```

3. Verify installation:
```bash
seednap --version
```

### Install for Development

To install with development dependencies:

```bash
pip install -e ".[dev]"
```

This includes pytest, black, ruff, mypy, and other development tools.

### Conda Environment (Recommended)

A conda environment specification is coming soon. For now, install dependencies manually:

```bash
conda create -n seednap python=3.9 cutadapt r-base r-tidyverse r-dada2
conda activate seednap
pip install -e .
```

**For ETH eDNA server users:** The existing conda image `metabarcoding` contains all necessary dependencies. 

## Quick Start

### 1. Create a Configuration File

Generate an example configuration file:

```bash
seednap init --marker teleo --output config/markers/my_analysis.yaml
```

Or use the provided example:

```bash
cp config/markers/teleo.yaml config/markers/my_analysis.yaml
```

Edit the configuration file to match your analysis:
- Update `paths.raw_data` to point to your FASTQ files
- Choose `taxonomy.method` (dada2, blast, ecotag, or decipher)
- Update database paths for your chosen method
- Adjust computational resources (`resources.max_cores`)

### 2. Validate Your Configuration

```bash
seednap validate config/markers/my_analysis.yaml
```

This checks that:
- YAML syntax is valid
- All required fields are present
- Values are within valid ranges
- Configuration is internally consistent

### 3. Run the Pipeline

**Note:** Full pipeline execution will be available after Phase 6. For now, use the legacy `main.sh` script.

```bash
# Coming soon in Phase 6:
seednap run config/markers/my_analysis.yaml

# For now, use the legacy script:
bash main.sh config/config_teleo.sh
```

## CLI Commands

### `seednap init`

Create an example configuration file:

```bash
seednap init --marker teleo --output config/markers/example.yaml
```

Options:
- `--marker, -m`: Marker name (default: teleo)
- `--output, -o`: Output path (default: config/markers/example.yaml)
- `--force, -f`: Overwrite existing file

### `seednap validate`

Validate a configuration file:

```bash
seednap validate config/markers/my_analysis.yaml
```

### `seednap version`

Show version information:

```bash
seednap version
```

### `seednap run` (Coming in Phase 6)

Run the full pipeline:

```bash
seednap run config/markers/my_analysis.yaml
```

Options:
- `--resume-from`: Resume from a specific step (trim, dada2, taxonomy, export)
- `--dry-run`: Show what would be run without executing

## Configuration

Configuration files use YAML format with the following main sections:

- **marker**: Marker name, primers, and description
- **paths**: Input/output directory paths
- **demultiplex**: Demultiplexing settings (if applicable)
- **trimming**: Cutadapt primer trimming parameters
- **dada2**: DADA2 filtering, merging, and chimera removal settings
- **taxonomy**: Taxonomic assignment method and database paths
- **export**: Output formats (CSV, GBIF)
- **metrics**: Quality control metrics configuration
- **logging**: Logging level and format
- **resources**: CPU cores, memory limits
- **pipeline**: Steps to execute

See [config/markers/teleo.yaml](config/markers/teleo.yaml) for a complete example with detailed comments.

---

## Legacy Usage (Shell Scripts)

The original shell-script based pipeline is still available during the migration. To process one marker, edit its config file (e.g., `config/config_teleo.sh`) and run:

```bash
bash main.sh path/to/config.sh
```

The content of the config file is as follows: it contains infos on paths and marker name

``` sh
# Config parameters
marker="teleo" # No uppercase here
raw_path="/home/shared/edna/raw/ma_ga_akand_2024/"
primer_F="ACACCGCCCGTCACTCT"
primer_R="CTTCCGGTACACTTACCATG"
method_demultiplex="" # among: primer_trim / ligation_trim - not ready yet
convert_for_GBIF="TRUE" # values= TRUE or FALSE

# Assignment parameters
method_assignment="dada2" # among: "dada2 / decipher / ecotag / blast" # Note: blast is in development

# For blast below
path_blast_fasta="/home/shared/edna/reference_database/2024/teleo/blast_db/refdb_all_fish_teleo.fasta" # If not already created, it will create the accompanying files in the same folder (ndb;nhr;nin etc)

# For ecotag below
path_ecotag_tree="/home/shared/edna/reference_database/2023_06/teleo_custom_embl/customtaxonomy/" # Directory not file # Caution: here it needs to be the path to the tree (NCBI tree ; -t ecotag option not -d)
path_ecotag_fasta="/home/shared/edna/reference_database/2023_06/teleo_custom_embl/db_teleo_custom_and_embl.fasta"

# For dada2 below
path_dada_all="utils/teleo_crabs_dada2_all.fasta"
path_dada_species="utils/teleo_crabs_dada2_species.fasta"

# For decipher below
path_decipher_trained="utils/teleo_trained.rds"

# General
CORES=12
```

**marker** is the marker name (lowercase)

**raw_path** is the location of the raw data  

**convert_for_GBIF** TRUE/FALSE: create an alternative output file. This would be the correct input format for creating the GBIF processed file. 

**CORES** indicates the number of cores to be used for a single cutadapt command

Some parameters are still in development. 
Note that to use the `ecotag` assignment option, you need to have a conda image containing the obitools named `obitools` or handle this part manually.  
This is because the obitools are not compatible with cutadapt within the same conda due to python version conflicts. 

## Content 

The script first cuts the primers from the sequences. It always cuts on the 5' end, and will primers on the 3' if they are present (but will not discard reads if there are no primers present in the 3' end).

Reads will be discarded is **no** ends are trimmed. If only trimmed on the 5', only trimmed on the 3' end, or trimmed on both, the reads will be kept. 
 
See the cutadapt documentation for more details: https://cutadapt.readthedocs.io/en/stable/guide.html

No direction precaution is taken as the samples are not prepared via a ligation protocol (checks were done and not reads were in the opposite direction)
Samples are already demultiplexed so we directly cut the primers (not anchorage since parts of the tags remain)

Next part of the script is launching dada2 in R. For now it processes all together, later dev will handle it on a library basis. 

Taxonomic assignment is done either with dada2, decipher, ecotag or blast (in dev). 

## Ligation-based compatibility

If necessary, the script demultiplex_ligation.sh can now be used to demultiplex and re-order reads from a pooled library created using the ligation method ("SPYGEN" type). Final output is saved in the original raw path, by creating a sub-directory `raw_for_dada2` in fastq format (not zipped). 

This script can be executed like this: 

```
bash scripts/demultiplex_ligation.sh [path_to_folder_with_libraries] [path_to_lab_file] [path_to_config.sh]
```

For for example: 

```
bash scripts/demultiplex_ligation.sh /home/shared/edna/raw/MA_Calanques_2023/ /home/shared/edna/raw/MA_Calanques_2023/metadata/metadata_lab_ma_fr_calanq_2023.csv config_teleo.sh
```

Config is the same as for the normal pipeline execution, and the metadata lab file is the one from our internal standards. Expected structure is: 


```
| `eventID`   | `tag_demultiplex` | `library`             | `pcr_primer_name_forward` | `pcr_primer_name_reverse` | `pcr_primer_forward`   | `pcr_primer_reverse`    |
|-------------|------------------|------------------------|----------------------------|----------------------------|------------------------|--------------------------|
| CNEG03_01   | aattgccg         | VL331___MB0523B1__     | teleoF                     | teleoR                     | ACACCGCCCGTCACTCT      | CTTCCGGTACACTTACCATG     |
| CNEG03_10   | ttaggcac         | VL331___MB0523B1__     | teleoF                     | teleoR                     | ACACCGCCCGTCACTCT      | CTTCCGGTACACTTACCATG     |
| CNEG03_11   | gtgttgga         | VL331___MB0523B1__     | teleoF                     | teleoR                     | ACACCGCCCGTCACTCT      | CTTCCGGTACACTTACCATG     |
| CNEG03_12   | aacgcgat         | VL331___MB0523B1__     | teleoF                     | teleoR                     | ACACCGCCCGTCACTCT      | CTTCCGGTACACTTACCATG     |
| CNEG03_02   | atgcttgg         | VL331___MB0523B1__     | teleoF                     | teleoR                     | ACACCGCCCGTCACTCT      | CTTCCGGTACACTTACCATG     |
| CNEG03_03   | atggaggt         | VL331___MB0523B1__     | teleoF                     | teleoR                     | ACACCGCCCGTCACTCT      | CTTCCGGTACACTTACCATG     |
| CNEG03_04   | tgaggaca         | VL331___MB0523B1__     | teleoF                     | teleoR                     | ACACCGCCCGTCACTCT      | CTTCCGGTACACTTACCATG     |
| CNEG03_05   | acaagacc         | VL331___MB0523B1__     | teleoF                     | teleoR                     | ACACCGCCCGTCACTCT      | CTTCCGGTACACTTACCATG     |
| CNEG03_06   | agttccac         | VL331___MB0523B1__     | teleoF                     | teleoR                     | ACACCGCCCGTCACTCT      | CTTCCGGTACACTTACCATG     |
```

## Issues with pattern recognition 

We have now implemented that the match of the sample name is strict, meaning that the script will only work if your input files are named [sample_name]_R1.fastq (or R2). Characters between sample name and R1/2 are no longer allowed, as this was too flexible and causing issues with diluted samples (with a pattern of sample_named). If you wish to use the older more flexible pattern, this is the code to replace in the main script

```
  r1_in=$(ls "${raw_path}" | grep "^${s}.*_R1\.fastq\.gz$")
  r2_in=$(ls "${raw_path}" | grep "^${s}.*_R2\.fastq\.gz$")
```

