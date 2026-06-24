# Installation

How to install SeeDNAP and the external bioinformatics tools it drives.

This page covers the recommended conda setup, manual installation, and the requirements you need for each pipeline step. On the ETH ELE eDNA server, installation is already done; see [On the ETH eDNA server](#on-the-eth-edna-server).

## On the ETH eDNA server

> [!IMPORTANT]
> On the ETH ELE eDNA server you do not need to install anything. SeeDNAP is already installed in a shared conda env that every user can activate:
> ```bash
> conda activate /home/shared/edna/envs/seednap
> seednap --version
> ```
> Always activate by full path. A bare `conda activate seednap` only works if you have a personal env of that name, and silently picks the wrong one if you do. All tools (cutadapt, vsearch, swarm, blast, R + packages) are present. OBITools for the optional `ecotag` method lives in its own env and is auto-discovered; no extra activation needed.
>
> The shared env is editable-installed from the canonical checkout, so it is always the current version. You do **not** clone the repo to run on the server, and you should **not** clone it per dataset (a private clone only goes stale and ships example configs pointing at other datasets). All you provide per dataset is a config file, and each run writes into the output folder named in it (`paths.output`). See the [server quickstart in the README](../README.md#quick-start) for the copy-a-config-and-run recipe.

The rest of this page is only for a fresh setup elsewhere (local development or a new machine).

## Conda environment (recommended)

The simplest way to install SeeDNAP with every dependency, including the external tools and R packages, is the bundled conda environment:

```bash
git clone https://github.com/WildinSync/wis_seednap.git
cd wis_seednap
conda env create -f environment.yml
conda activate seednap
pip install -e .
```

Verify the install:

```bash
seednap --version
```

> [!TIP]
> After the version check, scaffold and validate a config as a post-install sanity check:
> ```bash
> seednap init --marker teleo -o my_config.yaml   # write a starter config
> seednap validate my_config.yaml                 # schema + preflight checks
> ```
> `seednap validate` runs preflight: it fails if the config references files that are missing on disk or a taxonomy database block that cannot be resolved, so problems surface before a run starts. Run `seednap --help` to see all 17 commands; see [cli-reference.md](cli-reference.md).

## Manual installation

If you manage dependencies yourself, install the package and then the external tools separately. Target the versions pinned in `environment.yml` unless you have a reason not to.

```bash
pip install -e .
```

External tools and R packages (these are the exact channels and pins used in `environment.yml`):

```bash
# Primer trimming
conda install -c bioconda cutadapt=5.2

# SWARM OTU clustering + read merging, dereplication, chimera detection
# (chimeras are artefactual sequences formed when two real templates fuse during PCR)
conda install -c bioconda vsearch=2.30.5 swarm=3.1.6

# BLAST taxonomy
conda install -c bioconda blast=2.17.0

# DADA2 / DECIPHER R stack. r-base, r-tidyverse and r-patchwork come from
# conda-forge; the bioconductor-* packages come from bioconda, so give both channels.
conda install -c conda-forge -c bioconda r-base=4.2 r-tidyverse=2.0.0 r-patchwork=1.2.0 \
  bioconductor-biostrings=2.66.0 bioconductor-dada2=1.26.0 bioconductor-decipher=2.26.0
```

> [!WARNING]
> Install the full R stack, not just dada2 and DECIPHER. The R scripts also call `library(Biostrings)`, `library(dplyr)`, `library(ggplot2)` (the last two ship with `r-tidyverse`), and `library(patchwork)`; a partial install passes `pip install -e .` and Python import checks but fails at runtime the first time an R step loads a missing package. Note also that `environment.yml` does not install OBITools (the optional `ecotag` taxonomy method); set it up separately if you need it, per [ecotag-setup.md](ecotag-setup.md).

## Requirements

### Python dependencies

Defined in `pyproject.toml`:

| Package | Version | Used for |
|---|---|---|
| `pydantic` | >= 2.0 | Configuration validation |
| `pyyaml` | >= 6.0 | YAML config parsing |
| `click` | >= 8.0 | CLI framework |
| `pandas` | >= 2.2 | Data manipulation |
| `numpy` | >= 1.24 | Numerical operations |
| `biopython` | >= 1.80 | Sequence handling |
| `rich` | >= 13.0 | Console output formatting |
| `jinja2` | >= 3.0 | Template rendering |
| `matplotlib` | >= 3.5 | Charts for the optional HTML run report |
| `python-dotenv` | >= 1.0 | Loads `.env` (NCBI Entrez key for GBIF enrichment) |

> [!NOTE]
> The NCBI Entrez API key in `.env` is only needed for taxonomy enrichment in the `create-gbif` command. Copy `.env.example` to `.env` and fill in `NCBI_API_KEY` before running that command. No key is required for the core pipeline or for BLAST assignment. See [gbif-export.md](gbif-export.md).

### External tools

The conda environment pins each tool to the version validated on the ETH ELE eDNA server.

| Tool | Pinned version | Required for | Install |
|---|---|---|---|
| Cutadapt | 5.2 | Primer trimming | `conda install -c bioconda cutadapt=5.2` |
| VSEARCH | 2.30.5 | Read merging, dereplication, chimera detection | `conda install -c bioconda vsearch=2.30.5` |
| SWARM | 3.1.6 | OTU clustering (groups near-identical reads into operational taxonomic units) | `conda install -c bioconda swarm=3.1.6` |
| BLAST+ | 2.17.0 | BLAST taxonomic assignment | `conda install -c bioconda blast=2.17.0` |
| R | 4.2 | DADA2 (ASV inference) and DECIPHER taxonomy methods | `conda install -c conda-forge r-base=4.2` |

### R packages (DADA2 / DECIPHER)

Loaded by the R scripts and pinned in `environment.yml`. ASV inference and merging come from DADA2; DECIPHER is the alternative taxonomy classifier; the rest support sequence I/O and the diagnostic plots:

| conda package | Pinned version | Source |
|---|---|---|
| `bioconductor-dada2` | 1.26.0 | Bioconductor |
| `bioconductor-decipher` | 2.26.0 | Bioconductor |
| `bioconductor-biostrings` | 2.66.0 | Bioconductor |
| `r-tidyverse` (provides `dplyr`, `ggplot2`) | 2.0.0 | conda-forge |
| `r-patchwork` | 1.2.0 | conda-forge |

### ecotag (optional)

OBITools v1 has Python 2 dependencies that conflict with the rest of the environment, so it lives in its own conda env and is not installed by `environment.yml`. The `ecotag` runner discovers it in this order:

1. `PATH` -- an activated OBITools env wins, if `ecotag`, `obiannotate`, and `obitab` all resolve there.
2. `SEEDNAP_OBITOOLS_BIN` -- a bin directory used when the tools are not all on `PATH`.
3. A set of well-known install locations.

See [ecotag-setup.md](ecotag-setup.md) for setup.

## Running the test suite

The repository ships a pytest suite that exercises the post-processing logic without invoking any external bioinformatics tools, so it runs locally in a few seconds:

```bash
pytest
```

Tests cover, among other things:

- BLAST lowest-common-ancestor (LCA) behavior: when several reference sequences match a feature equally well, the assignment is collapsed to the lowest taxonomic rank they all agree on (e.g. genus rather than a guessed species), including the case where that resolves to nothing.
- The DarwinCore builder (DarwinCore is the GBIF biodiversity record standard that the `create-gbif` export targets).
- SWARM OTU non-zero invariants: a guard that the OTU (operational taxonomic unit, a cluster of near-identical reads treated as one feature) by sample table is not silently all zeros, which would mean reads were dropped during clustering.
- The shared taxonomy post-processor, demultiplex robustness, ecotag tool discovery, runner signatures, config-snapshot reproducibility, R-script packaging, and YAML round-tripping.

Add new tests under `tests/unit/` or `tests/integration/` when you change load-bearing logic.

## See also

- [cli-reference.md](cli-reference.md) -- all commands and options
- [configuration.md](configuration.md) -- YAML config reference
- [ecotag-setup.md](ecotag-setup.md) -- installing OBITools for the `ecotag` method
- [gbif-export.md](gbif-export.md) -- `create-gbif` and the NCBI Entrez key
