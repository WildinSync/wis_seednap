"""Guards for the packaging metadata: pyproject.toml, environment.yml, lockfile.

These are cheap structural checks that the version is sourced dynamically from
``seednap.__version__``, the load-bearing dependency floors are pinned (pandas
needs >=2.2 for ``DataFrameGroupBy.apply(include_groups=...)`` in the BLAST LCA),
the R scripts ship in the wheel, the project URLs point at the live repo, and
the validated server environment is captured as a reproducible lockfile. They
catch silent drift between environment.yml and pyproject.toml.
"""

from pathlib import Path

import pytest
import yaml

try:  # Python 3.11+ ships tomllib; 3.9/3.10 use the tomli backport.
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised only on <3.11
    tomli = pytest.importorskip("tomli")
    tomllib = tomli

from seednap.__version__ import __version__

REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_pyproject() -> dict:
    path = REPO_ROOT / "pyproject.toml"
    assert path.is_file(), "pyproject.toml must exist at the repo root"
    with path.open("rb") as fh:
        return tomllib.load(fh)


def _conda_deps() -> list[str]:
    """Flattened list of conda + pip dependency strings from environment.yml."""
    path = REPO_ROOT / "environment.yml"
    assert path.is_file(), "environment.yml must exist at the repo root"
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    flat: list[str] = []
    for entry in data["dependencies"]:
        if isinstance(entry, dict):
            flat.extend(entry.get("pip", []))
        else:
            flat.append(str(entry))
    return flat


def test_version_is_dynamic_from_version_module() -> None:
    """pyproject must source the version from seednap.__version__, not hardcode it."""
    data = _load_pyproject()
    project = data["project"]
    assert "version" in project.get("dynamic", []), "version must be declared dynamic"
    assert "version" not in project, "version must not be hardcoded in [project]"
    attr = data["tool"]["setuptools"]["dynamic"]["version"]["attr"]
    assert attr == "seednap.__version__.__version__"


def test_dynamic_version_resolves_to_module_value() -> None:
    """The setuptools attr resolver must yield the same string the module exposes."""
    from setuptools.config.pyprojecttoml import read_configuration

    cfg = read_configuration(str(REPO_ROOT / "pyproject.toml"))
    assert cfg["project"]["version"] == __version__


def test_pandas_floor_supports_include_groups() -> None:
    """pandas floor must be >=2.2 (BLAST LCA uses apply(include_groups=...))."""
    data = _load_pyproject()
    pandas_dep = next(
        d for d in data["project"]["dependencies"] if d.replace(" ", "").startswith("pandas")
    )
    assert ">=2.2" in pandas_dep, f"pandas floor must be >=2.2, got {pandas_dep!r}"


def test_python_dotenv_is_a_runtime_dependency() -> None:
    """python-dotenv loads the NCBI key from .env and must be a declared runtime dep."""
    data = _load_pyproject()
    assert any(d.startswith("python-dotenv") for d in data["project"]["dependencies"])
    assert any(d.startswith("python-dotenv") for d in _conda_deps())


def test_r_scripts_ship_as_package_data() -> None:
    """The DADA2/DECIPHER R scripts must be packaged so they ship in the wheel."""
    data = _load_pyproject()
    patterns = data["tool"]["setuptools"]["package-data"]["seednap"]
    assert "scripts/*.R" in patterns


def test_project_urls_point_at_live_repo() -> None:
    """[project.urls] must reference the live WildinSync/wis_seednap repository."""
    data = _load_pyproject()
    urls = data["project"]["urls"]
    assert urls["Repository"] == "https://github.com/WildinSync/wis_seednap"
    for value in urls.values():
        assert "WildinSync/wis_seednap" in value
        assert "eth-edna/seednap" not in value


def test_pinned_versions_consistent_between_files() -> None:
    """pandas and python-dotenv floors must agree across environment.yml and pyproject."""
    py = _load_pyproject()
    conda = _conda_deps()

    def floor(deps, name):
        match = next(d for d in deps if d.replace(" ", "").startswith(name))
        return match.split("#", 1)[0].strip()

    assert floor(py["project"]["dependencies"], "pandas") == floor(conda, "pandas") == "pandas>=2.2"


@pytest.mark.parametrize(
    "pin",
    [
        "bioconductor-dada2=1.26.0",
        "bioconductor-decipher=2.26.0",
        "bioconductor-biostrings=2.66.0",
        "r-tidyverse=2.0.0",
        "r-patchwork=1.2.0",
        "r-base=4.2",
    ],
)
def test_r_stack_is_pinned_to_validated_versions(pin: str) -> None:
    """The R/Bioconductor stack that produces the biology results must stay pinned."""
    assert pin in _conda_deps(), f"environment.yml must pin {pin}"


def test_lockfile_is_present_and_explicit() -> None:
    """environment.lock.txt must be a header-commented, explicit conda lockfile."""
    path = REPO_ROOT / "environment.lock.txt"
    assert path.is_file(), "environment.lock.txt must exist at the repo root"
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    assert lines[0].startswith("#"), "lockfile must open with a header comment"
    assert "validated" in lines[0].lower()
    assert "@EXPLICIT" in text, "must be an explicit (URL-pinned) conda lockfile"
    assert any(line.startswith("https://") for line in lines)
