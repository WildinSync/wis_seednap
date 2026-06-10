# Contributing to SeeDNAP

SeeDNAP is the eDNA metabarcoding pipeline run on the ETH ELE eDNA server to
process real biodiversity samples. Outputs are submitted to GBIF and shared
with collaborators, so the bar is correctness, reproducibility of any specific
run, and surviving environment drift on the server. Please keep changes
surgical and well-tested.

## Environment setup

SeeDNAP targets Python 3.9 and wraps external bioinformatics tools (cutadapt,
vsearch, swarm, blastn) and R packages (DADA2, DECIPHER). Those tools are
pinned in `environment.yml`, so a conda environment is the supported way to get
a working install.

### On the ETH ELE eDNA server

A shared, prebuilt environment already exists. Activate it; you do not need to
create your own or `pip install` anything:

```bash
conda activate /home/shared/edna/envs/seednap
```

### Fresh setup elsewhere (local development)

```bash
conda env create -f environment.yml
conda activate seednap
pip install -e .
```

`environment.yml` installs the package in editable mode (`-e .`) and the pinned
tool versions. For development tooling not pulled in by the conda recipe,
install the `dev` extras:

```bash
pip install -e ".[dev]"
```

## Running tests

Tests live under `tests/` (split into `tests/unit/` and `tests/integration/`)
and are run with pytest from the repo root:

```bash
pytest                 # whole suite
pytest tests/unit      # unit tests only
pytest -m integration  # integration tests only
```

Test discovery, markers (`unit`, `integration`), and default options are
configured in `pyproject.toml` under `[tool.pytest.ini_options]`.

Some integration tests touch real validation datasets on the server. If a
dataset is not present in your environment those tests will be skipped rather
than failing; do not work around a skip by faking inputs.

"It compiled" and "the import succeeded" are not verification. For load-bearing
changes (anything in `swarm/`, `dada2/`, `taxonomic_assignment/`, or the
orchestrator), run the pipeline end-to-end against one of the small validation
datasets and confirm the results stay consistent before opening a PR.

## Lint and type-check

```bash
ruff check .   # lint
mypy src       # type-check
```

Ruff and mypy are both configured in `pyproject.toml`. The package ships a
`py.typed` marker and is expected to stay mypy-clean (typed defs are required;
see `[tool.mypy]`). Formatting follows `black` (line length 100), also
configured in `pyproject.toml`.

> Note: a `pre-commit` package is listed in the `dev` extras, but there is no
> `.pre-commit-config.yaml` checked in, so there is no pre-commit hook step to
> run. Run `ruff check`, `mypy`, and `pytest` manually.

## Commit message convention

Commit titles carry a bracketed prefix that names the kind of change. Use one
of:

- `[FIX]` -- bug fix or correctness change
- `[FEAT]` -- new feature
- `[REFACTOR]` -- internal restructuring with no behavior change
- `[DOCS]` -- documentation only
- `[CONFIG]` -- changes to marker/config files or build configuration
- `[TEST]` -- tests only

Write the message to describe the change, not the author. Keep messages
professional; this repo is shared inside ELE and may be made public. Avoid
embedding result numbers or other specifics that rot when regenerated.

## Branch-then-PR flow

`main` is the default branch and the deployed branch. Do not commit directly to
`main`:

1. Branch off `main` for your change.
2. Make the change, keeping it surgical (every changed line should trace to the
   task), and add or update tests.
3. Run `ruff check .`, `mypy src`, and `pytest` locally.
4. Open a pull request against `main` on GitHub
   (`https://github.com/WildinSync/wis_seednap`).
5. Get review, then merge.

Do not push branches or merge proactively; confirm intent first.

## Deploy model

There is no separate release process. Production deployment is simply checking
out `main` and reinstalling the package in editable mode inside the conda
environment:

```bash
conda activate /home/shared/edna/envs/seednap   # or: conda activate seednap
git checkout main
git pull
pip install -e .
```

## No silent fallbacks

The single most important rule in this codebase: any fallback path must say so.
A silent default, a swallowed exception, or a quietly disabled feature can
produce a biodiversity dataset that looks valid and is wrong, and the samples
cannot be re-collected. Every fallback must emit a `[WARN]` line describing what
was expected, what happened, and what fallback was chosen, or `raise` with a
descriptive error naming the offending file or config key. See the project
guidelines for the exact log format.

## License

This project is licensed under the MIT License; see [LICENSE](LICENSE) at the
repo root.
