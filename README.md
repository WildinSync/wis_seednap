# seednap

**Modern eDNA metabarcoding pipeline with DADA2**

[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-passing-brightgreen.svg)](tests/)

A production-ready Python pipeline for processing environmental DNA (eDNA) metabarcoding data with support for multiple taxonomic assignment methods (DADA2, BLAST, DECIPHER, ecotag).

**Version:** 0.1.0

---

## Table of Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
  - [Complete Pipeline](#complete-pipeline)
  - [Individual Steps](#individual-steps)
- [CLI Reference](#cli-reference)
- [Configuration](#configuration)
- [Pipeline Steps](#pipeline-steps)
- [Taxonomic Assignment Methods](#taxonomic-assignment-methods)
- [Testing](#testing)
- [Architecture](#architecture)
- [Contributing](#contributing)
- [Citation](#citation)

---

## Pipeline Steps
1. **Trimming** (cutadapt): Remove primers and adapters
2. **DADA2**: Quality filtering, denoising, ASV detection, chimera removal
3. **Taxonomic Assignment**: BLAST, DADA2, DECIPHER, or ecotag
4. **Export**: GBIF-compatible format with rank determination

---

## Quick Start

### 1. Installation

```bash
# Clone repository
git clone https://gitlab.ethz.ch/ele-projects/edna/edna-app/seednap.git

cd seednap

# Install with pip
pip install -e .

# Verify installation
seednap --version
```

### 2. Create Configuration

```bash
# Generate example config
seednap init --marker teleo --output my_analysis.yaml

# Validate configuration
seednap validate my_analysis.yaml
```

### 3. Run Pipeline

```bash
# Run complete pipeline
seednap run-pipeline my_analysis.yaml

# Resume from failed step
seednap run-pipeline my_analysis.yaml --resume

# Stop on first error (default) or continue
seednap run-pipeline my_analysis.yaml --continue-on-error
```

---

## Installation

### Requirements

**Python:** >= 3.9

**System dependencies:**
- `cutadapt` (>= 4.0)
- `R` (>= 4.0)
- R packages: `tidyverse`, `dada2`, `Biostrings`

**Optional (for specific taxonomic methods):**
- **BLAST:** `ncbi-blast+` (makeblastdb, blastn)
- **DECIPHER:** R package `DECIPHER`
- **Ecotag:** `obitools` v1 (requires separate conda environment)

### Install from Source

```bash
# Clone repository
git clone https://gitlab.ethz.ch/ele-projects/edna/edna-app/seednap.git
cd seednap

# Install in development mode
pip install -e .

# Install with development dependencies
pip install -e ".[dev]"
```

### Conda Environment (Recommended)

```bash
# Create conda environment
conda env create -f environment.yml
conda activate seednap

# Install seednap
pip install -e .
```

**NOTE**: **For ETH eDNA server users:** Use the existing `metabarcoding` conda environment which contains all dependencies.

---

## Usage

### Complete Pipeline

Run the entire pipeline end-to-end:

```bash
# Basic usage
seednap run-pipeline config.yaml

# Resume from previous run
seednap run-pipeline config.yaml --resume

# Specify custom state file
seednap run-pipeline config.yaml --state-file pipeline_state.json

# Continue on errors (don't stop pipeline if one step fails)
seednap run-pipeline config.yaml --continue-on-error
```

**Output:**
- Trimmed FASTQ files
- DADA2 ASV table and sequences
- Taxonomic assignments
- GBIF-formatted export
- Quality plots and metrics
- Pipeline state file for resumability

### Individual Steps

Run pipeline steps separately:

#### 1. Primer Trimming

```bash
seednap trim config.yaml
```

Options:
- `--cores`: Number of CPU cores (default: from config)
- `--output-dir`: Output directory (default: from config)

#### 2. DADA2 Processing

```bash
seednap dada2 config.yaml
```

Performs:
- Quality filtering (maxEE, truncQ, maxN)
- Denoising and error learning
- Read merging
- Chimera removal
- ASV table generation

#### 3. Taxonomic Assignment

**BLAST:**
```bash
seednap blast config.yaml
```

**DADA2 RDP Classifier:**
```bash
seednap assign-taxonomy config.yaml --method dada2
```

**DECIPHER:**
```bash
seednap assign-taxonomy config.yaml --method decipher
```

**Ecotag:**
```bash
seednap assign-taxonomy config.yaml --method ecotag
```

#### 4. GBIF Formatting

```bash
seednap format-gbif input.csv output.csv --method dada2
```

Options:
- `--method`: Source method (dada2, blast, decipher, ecotag)
- `--add-rank/--no-add-rank`: Add rank column
- `--add-taxon/--no-add-taxon`: Add taxon column

---

## CLI Reference

### Global Options

```bash
seednap --help              # Show help
seednap --version           # Show version
seednap -v [command]        # Verbose output (DEBUG level)
seednap -q [command]        # Quiet mode (errors only)
```

### Commands

| Command | Description |
|---------|-------------|
| `init` | Create example configuration file |
| `validate` | Validate configuration file |
| `run-pipeline` | Run complete pipeline (trim → dada2 → taxonomy → export) |
| `trim` | Run primer trimming with cutadapt |
| `dada2` | Run DADA2 quality filtering and ASV detection |
| `blast` | Run BLAST taxonomic assignment with LCA resolution |
| `assign-taxonomy` | Run taxonomic assignment (dada2/decipher/ecotag) |
| `format-gbif` | Convert outputs to GBIF format |
| `demultiplex` | Demultiplex pooled libraries |
| `version` | Show detailed version information |

### Command Details

#### `seednap init`

Create an example configuration file:

```bash
seednap init [OPTIONS]

Options:
  -m, --marker TEXT     Marker name (default: teleo)
  -o, --output PATH     Output path (default: config/markers/example.yaml)
  -f, --force           Overwrite existing file
```

#### `seednap validate`

Validate a configuration file:

```bash
seednap validate CONFIG_FILE

Checks:
  - YAML syntax validity
  - Required fields present
  - Field types and values correct
  - Referenced paths exist
```

#### `seednap run-pipeline`

Run the complete pipeline:

```bash
seednap run-pipeline CONFIG_FILE [OPTIONS]

Options:
  --resume                    Resume from previous run
  --state-file PATH           Custom state file path
  --stop-on-error            Stop on first error (default)
  --continue-on-error        Continue pipeline if step fails
```

#### `seednap trim`

Run primer trimming:

```bash
seednap trim CONFIG_FILE [OPTIONS]

Options:
  --cores INTEGER       Number of CPU cores
  --output-dir PATH     Output directory
```

#### `seednap dada2`

Run DADA2 processing:

```bash
seednap dada2 CONFIG_FILE [OPTIONS]

Options:
  --input-dir PATH      Input directory (default: from config)
  --output-dir PATH     Output directory (default: from config)
```

#### `seednap blast`

Run BLAST taxonomic assignment:

```bash
seednap blast CONFIG_FILE [OPTIONS]

Options:
  --input PATH          ASV sequences file
  --database PATH       BLAST database FASTA
  --output PATH         Output CSV file
```

#### `seednap format-gbif`

Convert to GBIF format:

```bash
seednap format-gbif INPUT OUTPUT [OPTIONS]

Options:
  --method TEXT               Source method (dada2/blast/ecotag/decipher)
  --add-rank/--no-add-rank   Add rank column (default: yes)
  --add-taxon/--no-add-taxon Add taxon column (default: yes)
```

---

## Configuration

Configuration files use YAML format with the following structure:

### Example Configuration

```yaml
# Pipeline version
version: "0.1.0"

# Marker information
marker:
  name: "teleo"
  description: "Teleost fish eDNA metabarcoding using 12S rRNA"
  primers:
    forward: "ACACCGCCCGTCACTCT"
    reverse: "CTTCCGGTACACTTACCATG"
    name: "Teleo"
    target: "12S rRNA"
    amplicon_length: [100, 200]

# File paths
paths:
  raw_data: "/path/to/raw/fastq/files"
  output: "outputs"
  logs: "logs"
  references: "/path/to/reference/databases"

# Trimming configuration
trimming:
  tool: "cutadapt"
  min_length: 20
  max_error_rate: 0.1
  cores: 12
  discard_untrimmed: true

# DADA2 configuration
dada2:
  filter:
    max_ee: 2.0
    trunc_q: 11
    max_n: 0
  merge:
    min_overlap: 20
    max_mismatch: 0
  chimera:
    method: "consensus"

# Taxonomic assignment
taxonomy:
  method: "blast"  # Options: blast, dada2, decipher, ecotag

  databases:
    blast:
      fasta: "/path/to/blast/database.fasta"
      threshold_species: 98.0
      threshold_genus: 96.0
      threshold_family: 86.5

# Export configuration
export:
  formats: ["csv"]
  gbif:
    enabled: true
    add_rank: true
    add_taxon: true

# Resources
resources:
  max_cores: 12
  memory_limit: "32G"

# Pipeline steps
pipeline:
  steps: ["trim", "dada2", "taxonomy", "export"]
  skip: []
```

### Configuration Sections

| Section | Description |
|---------|-------------|
| `version` | Configuration format version |
| `marker` | Primer sequences and marker metadata |
| `paths` | Input/output directory paths |
| `demultiplex` | Demultiplexing settings (if needed) |
| `trimming` | Cutadapt primer trimming parameters |
| `dada2` | DADA2 filtering, merging, chimera settings |
| `taxonomy` | Taxonomic assignment method and databases |
| `export` | Output formats and GBIF settings |
| `metrics` | Quality control metrics configuration |
| `logging` | Logging level and format |
| `resources` | CPU cores, memory limits |
| `pipeline` | Steps to execute and skip |

See [config/markers/teleo.yaml](config/markers/teleo.yaml) for a complete annotated example.

---

## Pipeline Steps

### 1. Trimming (cutadapt)

**Input:** Raw FASTQ files (paired-end)
**Output:** Trimmed FASTQ files

**Process:**
- Remove forward and reverse primers
- Discard reads where primers not found (optional)
- Quality filtering by minimum length
- Generate trimming statistics

**Configuration:**
```yaml
trimming:
  min_length: 20          # Minimum read length after trimming
  max_error_rate: 0.1     # Maximum error rate in primer matching
  cores: 12               # CPU cores for parallel processing
  discard_untrimmed: true # Discard reads without primers
  overlap: 3              # Minimum overlap for primer detection
```

### 2. DADA2 Processing

**Input:** Trimmed FASTQ files
**Output:** ASV table, ASV sequences, quality plots

**Process:**
1. **Quality Filtering:**
   - Filter by expected errors (maxEE)
   - Truncate by quality score (truncQ)
   - Remove reads with N bases
   - Remove PhiX contamination

2. **Denoising:**
   - Learn error rates from data
   - Denoise reads to infer ASVs
   - Pool samples or process independently

3. **Merging:**
   - Merge paired-end reads
   - Require minimum overlap
   - Allow maximum mismatches

4. **Chimera Removal:**
   - Detect chimeric sequences
   - Methods: consensus, pooled, or none

**Configuration:**
```yaml
dada2:
  filter:
    max_ee: 2.0        # Maximum expected errors
    trunc_q: 11        # Truncate reads at first base with Q <= truncQ
    max_n: 0           # Maximum N bases allowed
    rm_phix: true      # Remove PhiX contamination
  merge:
    min_overlap: 20    # Minimum overlap for merging
    max_mismatch: 0    # Maximum mismatches in overlap
  chimera:
    method: "consensus" # Chimera detection method
  pool: false          # Pool samples for error learning
  multithread: true    # Use multiple threads
```

### 3. Taxonomic Assignment

Choose from four methods:

#### A. BLAST + LCA (Recommended)

**Process:**
1. Create BLAST database from reference FASTA
2. Run blastn search for each ASV
3. Extract phylogeny from database headers
4. Filter hits by percent identity thresholds
5. Resolve ambiguous assignments using Lowest Common Ancestor (LCA)

**Configuration:**
```yaml
taxonomy:
  method: "blast"
  databases:
    blast:
      fasta: "/path/to/refdb.fasta"
      perc_identity: 80.0       # Minimum percent identity
      qcov_hsp_perc: 80.0       # Minimum query coverage
      evalue: 1.0e-25           # Maximum e-value
      max_target_seqs: 5        # Maximum hits to return
      threshold_species: 98.0   # %ID threshold for species assignment
      threshold_genus: 96.0     # %ID threshold for genus assignment
      threshold_family: 86.5    # %ID threshold for family assignment
```

**Features:**
- Configurable percent identity thresholds by rank
- LCA resolution for ambiguous hits
- Phylogeny extraction from FASTA headers
- Handles multiple database formats

#### B. DADA2 RDP Classifier

**Process:**
1. Assign taxonomy using naive Bayesian classifier
2. Requires reference FASTA with taxonomy in headers

**Configuration:**
```yaml
taxonomy:
  method: "dada2"
  databases:
    dada2:
      all: "/path/to/dada2_all.fasta"
      species: "/path/to/dada2_species.fasta"
```

#### C. DECIPHER

**Process:**
1. Use trained DECIPHER model for assignment
2. Confidence-based taxonomic assignment

**Configuration:**
```yaml
taxonomy:
  method: "decipher"
  databases:
    decipher:
      trained: "/path/to/trained_model.rds"
      threshold: 60           # Confidence threshold
      processors: 8           # CPU cores
```

#### D. Ecotag (OBITools)

**Process:**
1. Use OBITools ecotag for assignment
2. Requires NCBI taxonomy tree

**Configuration:**
```yaml
taxonomy:
  method: "ecotag"
  databases:
    ecotag:
      tree: "/path/to/taxonomy/tree/"
      fasta: "/path/to/reference.fasta"
```

**Note:** Ecotag requires a separate conda environment (`obitools`) due to Python version conflicts with cutadapt.

### 4. Export (GBIF Formatting)

**Input:** Taxonomic assignment CSV
**Output:** GBIF-compatible CSV

**Process:**
1. Transform from wide format (samples as columns) to long format (samples as rows)
2. Determine taxonomic rank (species/genus/family/higher)
3. Extract lowest available taxon
4. Filter zero-count observations
5. Export to GBIF format

**Configuration:**
```yaml
export:
  formats: ["csv"]
  gbif:
    enabled: true
    add_rank: true   # Add 'rank' column
    add_taxon: true  # Add 'taxon' column
```

**Output Columns:**
- `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`
- `taxon`: Lowest available taxonomic assignment
- `rank`: Taxonomic rank of assignment (species/genus/family/higher)
- `sequence`: ASV sequence
- `nb_reads`: Read count
- `eventID`: Sample identifier

---

## Taxonomic Assignment Methods

### Method Comparison

| Method | Speed | Accuracy | Database | Best For |
|--------|-------|----------|----------|----------|
| **BLAST + LCA** | Moderate | High | Custom FASTA | Flexible, multiple hits |
| **DADA2 RDP** | Fast | Good | DADA2 format | Quick assignments |
| **DECIPHER** | Fast | Good | Trained model | Pre-trained databases |
| **Ecotag** | Slow | High | OBITools format | Legacy workflows |

### Choosing a Method

**Use BLAST + LCA when:**
- You need configurable percent identity thresholds
- You want LCA resolution for ambiguous hits
- You have a custom reference database
- You need detailed control over assignment logic

**Use DADA2 when:**
- You want fast, straightforward assignments
- Your database is in DADA2 format
- You're already using DADA2 for ASV detection

**Use DECIPHER when:**
- You have a pre-trained DECIPHER model
- You want confidence scores
- Speed is important

**Use ecotag when:**
- You're migrating from OBITools workflows
- You have existing ecotag databases

---

## Architecture

### Project Structure

```
seednap/
├── src/seednap/
│   ├── cli.py                    # Command-line interface
│   ├── config/                   # Configuration management
│   │   ├── models.py            # Pydantic config models
│   │   ├── loader.py            # YAML loading and validation
│   │   └── __init__.py
│   ├── pipeline/                 # Pipeline orchestration
│   │   ├── state.py             # State tracking and persistence
│   │   ├── orchestrator.py      # Pipeline execution
│   │   └── __init__.py
│   ├── steps/                    # Pipeline step implementations
│   │   ├── trimming/            # Cutadapt integration
│   │   ├── dada2/               # DADA2 processing
│   │   ├── taxonomic_assignment/ # BLAST, DECIPHER, ecotag
│   │   └── formatting/          # GBIF formatting
│   └── utils/                    # Shared utilities
│       ├── logging.py           # Rich console logging
│       └── sequences.py         # Sequence manipulation
├── tests/                        # Test suite
├── config/                       # Example configurations
│   └── markers/
│       └── teleo.yaml           # Teleo marker example
├── docs/                         # Documentation (coming soon)
├── environment.yml               # Conda environment
├── pyproject.toml               # Python package metadata
└── README.md                     # This file
```

### State Management

The pipeline tracks state in JSON format:

```json
{
  "marker": "teleo",
  "started_at": "2025-01-15T10:30:00",
  "current_step": "taxonomy",
  "steps": {
    "trim": {
      "status": "completed",
      "started_at": "2025-01-15T10:30:00",
      "completed_at": "2025-01-15T10:45:00",
      "duration_seconds": 900,
      "outputs": {
        "trimmed_dir": "/path/to/trimmed"
      }
    },
    "dada2": {
      "status": "completed",
      ...
    },
    "taxonomy": {
      "status": "running",
      ...
    }
  }
}
```

This allows:
- Resume from any failed step
- Track step timing and duration
- Pass outputs between steps
- Debug pipeline issues

---

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

```bash
# Install in development mode with dev dependencies
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install

# Run tests
pytest

# Run linters
ruff check .
black --check .
mypy src/
```

---

## Acknowledgments

- DADA2 pipeline and algorithm: Callahan et al. (2016)
- cutadapt: Martin (2011)
- BLAST: Altschul et al. (1990)
- DECIPHER: Wright (2016)
- OBITools: Boyer et al. (2016)