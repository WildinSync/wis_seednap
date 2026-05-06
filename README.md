<p align="center">
  <img src="media/teaser.png" alt="SeeDNAP">
</p>

**Modern eDNA metabarcoding pipeline with DADA2 and SWARM**

[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## What is SeeDNAP?

SeeDNAP is an end-to-end Python pipeline for processing environmental DNA (eDNA) metabarcoding data. It takes raw paired-end FASTQ files and produces taxonomically assigned OTU/ASV tables ready for biodiversity analysis or GBIF submission.

```mermaid
flowchart LR
    raw[Raw FASTQ<br/>Paired-end R1 / R2]:::io
    trim[Trim<br/>Cutadapt 2-pass]:::pass

    subgraph s2 [STEP 2 : Cluster - pick one]
        direction TB
        swarm[SWARM<br/>VSEARCH + SWARM]:::recommended
        dada2[DADA2<br/>R / Bioconductor]:::alt
    end

    subgraph s3 [STEP 3 : Taxonomy - pick one]
        direction TB
        blast[BLAST + LCA<br/>Recommended]:::recommended
        rdp[DADA2 RDP<br/>Naive Bayesian]:::alt
        decipher[DECIPHER<br/>IdTaxa classifier]:::alt
        ecotag[ecotag<br/>OBITools]:::alt
    end

    export[Export<br/>GBIF / DarwinCore]:::pass
    final[Final CSV<br/>Taxonomy + Abundances]:::io

    raw --> trim --> s2 --> s3 --> export --> final

    classDef io fill:#ffffff,stroke:#bbbbbb,stroke-width:1px,color:#000
    classDef pass fill:#d4edda,stroke:#28a745,stroke-width:2px,color:#000
    classDef recommended fill:#d4edda,stroke:#28a745,stroke-width:2px,color:#000
    classDef alt fill:#ffffff,stroke:#bbbbbb,stroke-width:1px,color:#000
```

## Quick Start

```bash
# Install
git clone https://github.com/WildinSync/wis_seednap.git
cd wis_seednap
conda env create -f environment.yml
conda activate seednap
pip install -e .

# Create and edit a config
seednap init --marker teleo --output config/markers/my_marker.yaml

# Run the pipeline
seednap run-pipeline config/markers/my_marker.yaml
```

That's it. See [docs/](docs/) for configuration details, step-by-step guides, and CLI reference.

## Requirements

| Tool | Pinned Version | Purpose |
|---|---|---|
| Python | 3.9 | Pipeline runtime |
| Cutadapt | 5.2 | Primer trimming |
| VSEARCH | 2.30.5 | Read merging, dereplication, chimera detection |
| SWARM | 3.1.6 | OTU clustering |
| BLAST+ | 2.17.0 | Taxonomic assignment |
| R | 4.2 | DADA2 / DECIPHER (optional) |

External tool versions are pinned in `environment.yml` to the set we validate against. OBITools (for the optional `ecotag` method) lives in a separate env -- see [docs/ecotag-setup.md](docs/ecotag-setup.md).

## Pipeline Steps

| Step | Tool | Description |
|---|---|---|
| **Demultiplex** *(optional)* | Built-in | Ligation-tag demultiplexing; `skip: true` for pre-demultiplexed inputs |
| **Trim** | Cutadapt | Two-pass primer removal (5' then 3') |
| **Cluster** | SWARM or DADA2 | OTU clustering or ASV denoising |
| **Taxonomy** | BLAST, DADA2, DECIPHER, or ecotag | Taxonomic assignment with cascade-null per-rank thresholds and MEGAN-LR top-bitscore LCA (BLAST) or RDP bootstrap (DADA2) |
| **Export** | Built-in | GBIF long format and DarwinCore occurrence CSV with deterministic `occurrenceID` and `contamination_flag` |

## CLI Commands

| Command | Description |
|---|---|
| `run-pipeline CONFIG` | Run the full pipeline from a YAML config |
| `init` | Generate an example config file |
| `validate CONFIG` | Validate a config file |
| `trim INPUT_DIR` | Primer trimming with Cutadapt |
| `swarm MARKER READS_DIR` | SWARM OTU clustering |
| `dada2 MARKER READS_DIR` | DADA2 ASV processing |
| `blast QUERY REF COUNTS` | BLAST taxonomic assignment with LCA |
| `assign-taxonomy METHOD MARKER QUERY COUNTS` | Generic taxonomy (blast/dada2/decipher/ecotag) |
| `format-gbif INPUT` | Convert results to GBIF long format |
| `create-gbif TAXO SAMPLE_META PROJECT_META OUTPUT` | Build DarwinCore GBIF occurrence CSV |
| `demultiplex READS LIB META` | Demultiplex ligation-based libraries |

Run `seednap --help` or `seednap <command> --help` for full options.

## Configuration

Everything is controlled by a single YAML file per marker. Example configs are in [config/markers/](config/markers/). Key sections:

```yaml
marker:
  name: "teleo"
  primers:
    forward: "ACACCGCCCGTCACTCT"
    reverse: "CTTCCGGTACACTTACCATG"

paths:
  raw_data: "/path/to/fastq/files"
  output: "outputs"

pipeline:
  steps: ["trim", "swarm", "taxonomy"]
```

Full configuration reference: [docs/configuration.md](docs/configuration.md)

## Documentation

| Document | Description |
|---|---|
| [docs/installation.md](docs/installation.md) | Installation and environment setup |
| [docs/configuration.md](docs/configuration.md) | Complete YAML configuration reference |
| [docs/pipeline-steps.md](docs/pipeline-steps.md) | Detailed description of each pipeline step |
| [docs/cli-reference.md](docs/cli-reference.md) | Full CLI command reference |
| [docs/taxonomy-methods.md](docs/taxonomy-methods.md) | Taxonomy assignment methods comparison |
| [docs/gbif-export.md](docs/gbif-export.md) | GBIF and DarwinCore export guide |
| [docs/ecotag-setup.md](docs/ecotag-setup.md) | OBITools / ecotag installation and discovery |

## Project Structure

```
seednap/
  src/seednap/
    cli.py                          # CLI entry point
    config/                         # Pydantic config models + YAML loader
    pipeline/                       # Orchestrator + state management
    steps/
      trimming/                     # Cutadapt integration
      dada2/                        # DADA2 R wrapper
      swarm/                        # VSEARCH + SWARM clustering
      taxonomic_assignment/         # BLAST, DADA2, DECIPHER, ecotag
      formatting/                   # GBIF + DarwinCore export
    utils/                          # Subprocess, logging, sequence tools
  config/markers/                   # Example YAML configs
  scripts/                          # R scripts (DADA2, DECIPHER)
```

## Acknowledgments

SeeDNAP builds on: [Cutadapt](https://cutadapt.readthedocs.io/) (Martin, 2011), [VSEARCH](https://github.com/torognes/vsearch) (Rognes et al., 2016), [SWARM](https://github.com/torognes/swarm) (Mahe et al., 2015), [BLAST+](https://blast.ncbi.nlm.nih.gov/) (Camacho et al., 2009), [DADA2](https://benjjneb.github.io/dada2/) (Callahan et al., 2016).

## License

MIT. See [LICENSE](LICENSE).
