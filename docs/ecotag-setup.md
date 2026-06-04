# ecotag (OBITools) setup

The ecotag taxonomy method uses OBITools v1, which has Python 2 dependencies
that conflict with Python 3 tools (Cutadapt 5.x, DADA2 R packages, etc.).
For that reason OBITools is **not** installed in the main `seednap` conda
environment -- it lives in its own env.

The seednap pipeline auto-discovers OBITools in three ways, in this order:

1. **`SEEDNAP_OBITOOLS_BIN` environment variable** -- explicit override.
2. **`PATH`** -- whatever's already activated (`conda activate obitools`).
3. **Well-known install locations** -- `/opt/anaconda3/envs/obitools/bin`,
   `/opt/conda/envs/obitools/bin`, `~/miniconda3/envs/obitools/bin`,
   `~/.conda/envs/obitools/bin`, `~/anaconda3/envs/obitools/bin`.

If none of those work, the runner emits a clear error pointing here.

## On the ETH ELE eDNA server (recommended)

OBITools is already installed at `/opt/anaconda3/envs/obitools`. No setup
needed; the runner will find it automatically. Just run:

```bash
seednap run-pipeline config/markers/mymarker.yaml
```

## On a fresh machine

Install OBITools v1 in a separate conda env:

```bash
conda create -n obitools -c bioconda obitools -y
```

Then either activate it or point seednap at it:

```bash
# Option 1: activate before running seednap
conda activate obitools

# Option 2: point seednap at it without activating
export SEEDNAP_OBITOOLS_BIN=$(conda info --base)/envs/obitools/bin
```

Verify the runner can find it:

```python
from seednap.steps.taxonomic_assignment.ecotag_runner import EcotagRunner
EcotagRunner()  # raises EcotagError with install instructions if not found
```

## Why a separate env?

OBITools v1's installation chain:

- `obitools` 1.x package -> Python 2.7
- Cutadapt 5.x (used in seednap's trim step) -> Python 3.10+
- These can't coexist in the same conda environment.

OBITools v3 / OBITools4 are the maintained Python 3 / Go rewrites, but they
use a different command set (`obitag` instead of `ecotag`) and have not been
benchmarked against the same reference databases as v1 in this lab. For now
seednap targets OBITools v1; migration to v4 is a separate work stream.

## Reference databases

ecotag needs an NCBI-format taxonomy tree (a directory of `.tdx`/`.adx` files
from `obitaxonomy --download-ncbi-taxdump`) and a reference sequence FASTA in
OBITools format. On the server these live at

```
/home/shared/edna/reference_database/2023_06/teleo_custom_embl/customtaxonomy/
/home/shared/edna/reference_database/2023_06/teleo_custom_embl/db_teleo_custom_and_embl.fasta
```

and are referenced from `config/markers/teleo.yaml` under
`taxonomy.databases.ecotag`.
