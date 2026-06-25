# ecotag (OBITools) setup

<img src="../media/divider.svg" width="100%" alt="">

How to install and configure OBITools v1 so seednap can run the `ecotag` taxonomy method.

ecotag is the taxonomic-assignment program from OBITools, a long-standing DNA-metabarcoding toolkit. It compares each query sequence (an ASV or OTU, the per-marker sequence variants the pipeline produces) against a reference database and assigns it to the lowest taxonomic rank still consistent with its closest matches. When several references tie, it falls back to their last (lowest) common ancestor (LCA), so an ambiguous sequence is reported at, say, genus or family rather than forced to a single species.

OBITools v1 has Python 2 dependencies and therefore lives in its own conda environment, separate from the main `seednap` env. This doc covers the install, how seednap discovers it, and the config keys it needs.

ecotag only runs when a marker config sets `taxonomy.method: "ecotag"` **and** lists `"taxonomy"` in `pipeline.steps`. The shipped configs default to `method: "blast"` (see `config/markers/teleo.yaml`), so `seednap run-pipeline <config>` alone does not exercise ecotag. Switch the method first.

## 🔍 How seednap finds OBITools

The runner auto-discovers the OBITools bin directory by probing these sources in order, requiring all three of `ecotag`, `obiannotate`, and `obitab`:

| Order | Source | When used |
| :--: | --- | --- |
| 1 | `PATH` | All three tools resolve on PATH (e.g. after `conda activate obitools`). The directory of the first tool is used. |
| 2 | `SEEDNAP_OBITOOLS_BIN` env var | A bin directory used when the tools are not all on PATH. |
| 3 | Well-known install locations | `/opt/anaconda3/envs/obitools/bin`, `/opt/conda/envs/obitools/bin`, `~/miniconda3/envs/obitools/bin`, `~/.conda/envs/obitools/bin`, `~/anaconda3/envs/obitools/bin`. |

PATH wins over `SEEDNAP_OBITOOLS_BIN` only when `ecotag`, `obiannotate`, and `obitab` all resolve on PATH; otherwise the env-var fallback and well-known locations are tried in turn. If none contain all three binaries, the runner raises `EcotagError` with inline install instructions (the `conda activate` and `SEEDNAP_OBITOOLS_BIN` options) and the list of probed locations.

## 🖥️ On the ETH ELE eDNA server

OBITools is already installed at `/opt/anaconda3/envs/obitools`. No setup is needed; the runner finds it automatically. With a config that selects ecotag (see above):

```bash
seednap run-pipeline config/markers/mymarker.yaml
```

## 🔧 On a fresh machine

Install OBITools v1 in its own conda env:

```bash
conda create -n obitools -c bioconda obitools -y
```

Then either activate it, or point seednap at it without activating:

```bash
# Option 1: activate before running seednap
conda activate obitools

# Option 2: point seednap at the env's bin directory
export SEEDNAP_OBITOOLS_BIN=$(conda info --base)/envs/obitools/bin
```

> [!WARNING]
> OBITools v1 cannot share an env with the main seednap tools: its Python 2.7 dependencies conflict with seednap's Python 3 stack (Cutadapt 5.x in the trim step). It MUST be a separate conda env. OBITools v3/v4 are the maintained Python 3 / Go rewrites but use a different command set (`obitag`, not `ecotag`) and have not been benchmarked against the v1 reference databases in this lab, so seednap targets v1.

Verify the runner can find it:

```python
from seednap.steps.taxonomic_assignment.ecotag_runner import EcotagRunner

EcotagRunner()                          # auto-discovers; raises EcotagError if not found
EcotagRunner(bin_dir="/path/to/bin")    # or pass the bin directory explicitly
```

`EcotagRunner` also accepts `timeout` (seconds per OBITools command, default `3600` = 1 hour).

## 📂 Reference databases

ecotag needs two inputs, configured under `taxonomy.databases.ecotag` in the marker YAML. The first is the taxonomy backbone (the rank hierarchy that lets ecotag walk up to a common ancestor), the second is the labelled reference sequences themselves:

| Key | Type | Default | Meaning |
| --- | --- | :--: | --- |
| `tree` | Path | required | Directory holding the NCBI taxonomy tree in OBITools binary form (the `.tdx`/`.adx`/`.ndx` files built by `obitaxonomy` from an NCBI taxdump). This encodes which taxa nest inside which, so ecotag can resolve an LCA. Passed to `ecotag -t`. |
| `fasta` | Path | required | Reference sequence FASTA in OBITools format: each record's header carries the NCBI taxid of the organism it came from, which is how a match becomes a taxonomic assignment. Passed to `ecotag -R`. |

On the server these are, for the teleo marker:

```yaml
taxonomy:
  method: "ecotag"
  databases:
    ecotag:
      tree: "/home/shared/edna/reference_database/2023_06/teleo_custom_embl/customtaxonomy/"
      fasta: "/home/shared/edna/reference_database/2023_06/teleo_custom_embl/db_teleo_custom_and_embl.fasta"
```

<details>
<summary><b>Validate database paths before a full run</b></summary>

Run `seednap validate <config>` before a full run. The Pydantic model only path-expands `tree` and `fasta` at load time; the model check itself does not confirm the paths exist. Existence is checked separately by a read-only preflight step that both `seednap validate` and the start of `seednap run-pipeline` run: it lists each database path as found or MISSING in the summary table and, if any is missing, prints a what/why/fix error and exits with status 1 before trimming or clustering begins. The `FileNotFoundError` inside the ecotag runner is only a last-resort safety net for paths that vanish mid-run.

</details>

## 📖 See also

- [taxonomy-methods.md](taxonomy-methods.md): all taxonomy methods, output schema, and the `is_contaminant_candidate` flag from `taxonomy.contaminants`.
- [configuration.md](configuration.md): full config reference.
