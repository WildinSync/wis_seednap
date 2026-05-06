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

You must separately install the external tools (versions match
`environment.yml`):

```bash
# Primer trimming
pip install cutadapt==5.2

# SWARM clustering
conda install -c bioconda vsearch=2.30.5 swarm=3.1.6

# BLAST taxonomy
conda install -c bioconda blast=2.17.0

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

The conda environment pins each tool to the version we validate against on
the ETH ELE eDNA server. If you install manually, target these versions
unless you have a reason not to.

| Tool | Pinned Version | Required For | Install |
|---|---|---|---|
| Cutadapt | 5.2 | Primer trimming | `pip install cutadapt==5.2` |
| VSEARCH | 2.30.5 | Read merging, dereplication, chimera detection | `conda install -c bioconda vsearch=2.30.5` |
| SWARM | 3.1.6 | OTU clustering | `conda install -c bioconda swarm=3.1.6` |
| BLAST+ | 2.17.0 | BLAST taxonomic assignment | `conda install -c bioconda blast=2.17.0` |
| R | 4.2 | DADA2 and DECIPHER methods | `conda install r-base=4.2` |

**ecotag (optional).** OBITools v1 has Python 2 dependencies that conflict
with the rest of the environment, so it lives in its own conda env. The
runner auto-discovers it via `SEEDNAP_OBITOOLS_BIN`, `PATH`, or a set of
well-known install paths. See [ecotag-setup.md](ecotag-setup.md).

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
OBITools (for the optional `ecotag` method) lives in
`/opt/anaconda3/envs/obitools` and is auto-discovered; no extra activation
needed.

## Running the Test Suite

The repository ships with a pytest suite that exercises the post-processing
logic without invoking any external bioinformatics tools, so it runs
locally in under a second.

```bash
pytest
```

Tests cover BLAST LCA / cascade-null behavior, the DarwinCore builder,
SWARM OTU non-zero invariants, the shared taxonomy post-processor,
demultiplex robustness, ecotag discovery, runner signatures, utility
coverage, and YAML round-tripping. Add new tests under `tests/unit/` or
`tests/integration/` when you change load-bearing logic.
