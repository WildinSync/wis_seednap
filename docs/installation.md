# Installation

## Conda Environment (Recommended)

The simplest way to install SeeDNAP with all dependencies:

```bash
git clone https://gitlab.ethz.ch/ele-projects/edna/edna-app/seednap.git
cd seednap
conda env create -f environment.yml
conda activate seednap
pip install -e .
```

Verify:

```bash
seednap --version
```

## Manual Installation

If you prefer to manage dependencies yourself:

```bash
pip install -e .
```

You must separately install the external tools:

```bash
# Primer trimming
pip install cutadapt

# SWARM clustering
conda install -c bioconda vsearch swarm

# BLAST taxonomy
conda install -c bioconda blast

# DADA2 / DECIPHER (R packages)
conda install -c bioconda bioconductor-dada2 bioconductor-decipher
```

## Requirements

### Python Dependencies

Defined in `pyproject.toml`:

- `pydantic >= 2.0` -- Configuration validation
- `pyyaml >= 6.0` -- YAML config parsing
- `click >= 8.0` -- CLI framework
- `pandas >= 2.0` -- Data manipulation
- `numpy >= 1.24` -- Numerical operations
- `biopython >= 1.80` -- Sequence handling
- `rich >= 13.0` -- Console output formatting
- `jinja2 >= 3.0` -- Template rendering
- `python-dotenv >= 1.0` -- Environment variable loading

### External Tools

| Tool | Min. Version | Required For | Install |
|---|---|---|---|
| Cutadapt | 4.0 | Primer trimming | `pip install cutadapt` |
| VSEARCH | 2.0 | Read merging, dereplication, chimera detection | `conda install -c bioconda vsearch` |
| SWARM | 3.0 | OTU clustering | `conda install -c bioconda swarm` |
| BLAST+ | 2.12 | BLAST taxonomic assignment | `conda install -c bioconda blast` |
| R | 4.0 | DADA2 and DECIPHER methods | `conda install r-base=4.2` |

### R Packages (for DADA2/DECIPHER)

- `dada2` (Bioconductor)
- `DECIPHER` (Bioconductor)
- `Biostrings` (Bioconductor)
- `tidyverse`
- `patchwork`

All R packages are included in `environment.yml`.

## ETH eDNA Server Users

On the ETH ELE eDNA server, use the pre-configured conda environment:

```bash
conda activate seednap
```

All tools (cutadapt, vsearch, swarm, blast, R + packages) are already installed.
