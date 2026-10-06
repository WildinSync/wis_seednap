"""Unit tests for the ``seednap init`` templates (minimal / full).

The templates are hand-written YAML strings in config/loader.py, so nothing ties them to the
Pydantic config models. These tests do: when a parameter is added, removed or renamed, or a
default changes, they fail and name the template line that needs updating.
"""

from __future__ import annotations

import re
import typing
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner
from pydantic import BaseModel

from seednap.cli import main
from seednap.config import create_example_config, load_config
from seednap.config.loader import INIT_TEMPLATES, _lookup_primers
from seednap.config.models.pipeline import PipelineConfig
from seednap.config.models.taxonomy import (
    BlastDatabaseConfig,
    Dada2DatabaseConfig,
    EcotagDatabaseConfig,
)

# Template values that deliberately differ from the model default. If the model default
# changes to match, or the template value changes, update this table.
INTENTIONAL_OVERRIDES = {
    "demultiplex.protocol": "ligation",  # lab protocol; model default is "none"
    "trimming.cores": 4,  # model default is 1
    "pipeline.steps": ["trim", "dada2", "taxonomy", "clean", "export", "report"],
}

# Filled with example values or /path/to placeholders rather than the model default.
# (Required fields have no default and are skipped automatically.)
EXAMPLE_VALUES = {
    "marker.description",
    "paths.raw_data",
    "demultiplex.metadata",
    "taxonomy.databases",
    "taxonomy.contaminants",
}

DB_MODELS = {
    "blast": BlastDatabaseConfig,
    "dada2": Dada2DatabaseConfig,
    "ecotag": EcotagDatabaseConfig,
}


def _render(tmp_path: Path, template: str, marker: str = "teleo") -> tuple[str, dict]:
    path = tmp_path / f"{template}.yaml"
    create_example_config(path, marker=marker, template=template)
    text = path.read_text(encoding="utf-8")
    return text, yaml.safe_load(text)


def _submodel(annotation) -> type[BaseModel] | None:
    for arg in typing.get_args(annotation) or (annotation,):
        if isinstance(arg, type) and issubclass(arg, BaseModel):
            return arg
    return None


def _walk(model: type[BaseModel], data: dict, text: str, prefix: str = "") -> list[str]:
    """Compare a template mapping against a model. Returns human-readable problems."""
    problems = []
    fields = {(f.alias or name): f for name, f in model.model_fields.items()}

    for key in data:
        if key not in fields:
            problems.append(
                f"{prefix}{key}: in the template but not in {model.__name__} (stale key)"
            )

    for key, field in fields.items():
        dotted = f"{prefix}{key}"
        if key not in data:
            # Optional fields without a default may be shown commented out ("# key: ...").
            commented = re.search(rf"^\s*#\s*{re.escape(key)}:", text, re.MULTILINE)
            if field.is_required() or field.get_default() is not None or not commented:
                problems.append(f"{dotted}: missing from the full template")
            continue

        sub = _submodel(field.annotation)
        if sub is not None and isinstance(data[key], dict):
            problems += _walk(sub, data[key], text, dotted + ".")
            continue
        if field.is_required() or dotted in EXAMPLE_VALUES:
            continue

        default = field.get_default(call_default_factory=True)
        value = data[key]
        if isinstance(default, Path):
            default = str(default)
        if dotted in INTENTIONAL_OVERRIDES:
            if value != INTENTIONAL_OVERRIDES[dotted]:
                problems.append(
                    f"{dotted}: template has {value!r}, INTENTIONAL_OVERRIDES says "
                    f"{INTENTIONAL_OVERRIDES[dotted]!r}"
                )
            if value == default:
                problems.append(
                    f"{dotted}: now equals the model default; drop it from INTENTIONAL_OVERRIDES"
                )
        elif value != default:
            problems.append(f"{dotted}: template has {value!r}, model default is {default!r}")
    return problems


def test_full_template_matches_models(tmp_path: Path) -> None:
    text, data = _render(tmp_path, "full")
    problems = _walk(PipelineConfig, data, text)
    assert not problems, "Update the full template in config/loader.py:\n  " + "\n  ".join(
        problems
    )


def test_full_template_database_blocks_match_models(tmp_path: Path) -> None:
    text, data = _render(tmp_path, "full")
    blast = data["taxonomy"]["databases"]["blast"]
    problems = _walk(BlastDatabaseConfig, blast, text, "taxonomy.databases.blast.")
    # The other methods are shown commented out: every field must still be listed.
    for method in ("dada2", "ecotag"):
        for key in DB_MODELS[method].model_fields:
            if not re.search(rf"^\s*#\s+{key}:", text, re.MULTILINE):
                problems.append(
                    f"taxonomy.databases.{method}.{key}: missing from the commented block"
                )
    assert not problems, "Update the full template in config/loader.py:\n  " + "\n  ".join(
        problems
    )


def test_minimal_template_keys_exist_in_full(tmp_path: Path) -> None:
    _, minimal = _render(tmp_path, "minimal")
    _, full = _render(tmp_path, "full")

    def keys(d: dict, prefix: str = "") -> set[str]:
        out = set()
        for k, v in d.items():
            out.add(prefix + k)
            if isinstance(v, dict):
                out |= keys(v, f"{prefix}{k}.")
        return out

    assert keys(minimal) <= keys(full)


@pytest.mark.parametrize("template", INIT_TEMPLATES)
def test_steps_come_first_without_demultiplex(tmp_path: Path, template: str) -> None:
    _, data = _render(tmp_path, template)
    assert next(iter(data)) == "pipeline"
    assert "demultiplex" not in data["pipeline"]["steps"]


@pytest.mark.parametrize("template", INIT_TEMPLATES)
def test_template_loads_once_paths_are_set(tmp_path: Path, template: str) -> None:
    text, _ = _render(tmp_path, template)
    raw = tmp_path / "raw"
    raw.mkdir()
    ref = tmp_path / "ref.fasta"
    ref.write_text(">a\nACGT\n")
    text = text.replace("/path/to/raw_fastq", str(raw)).replace(
        "/path/to/teleo_reference.fasta", str(ref)
    )
    text = re.sub(
        r'^(\s*(output|logs):\s*)"(\w+)"', rf'\1"{tmp_path}/\3"', text, flags=re.MULTILINE
    )
    cfg_path = tmp_path / "filled.yaml"
    cfg_path.write_text(text)

    cfg = load_config(cfg_path, validate=True)
    assert cfg.marker.name == "teleo"
    assert cfg.taxonomy.method == "blast"


def test_known_marker_gets_its_primers(tmp_path: Path) -> None:
    forward, reverse, found = _lookup_primers("mifish")
    assert found
    text, data = _render(tmp_path, "minimal", marker="mifish")
    assert data["marker"]["primers"] == {"forward": forward, "reverse": reverse}
    assert "placeholder" not in text


def test_unknown_marker_gets_placeholder_primers(tmp_path: Path) -> None:
    text, data = _render(tmp_path, "minimal", marker="not_a_marker")
    assert data["marker"]["name"] == "not_a_marker"
    assert text.count("# placeholder") == 2


def test_cli_init_default_output_and_force(tmp_path: Path) -> None:
    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path):
        assert runner.invoke(main, ["init"]).exit_code == 0
        assert Path("teleo.yaml").exists()

        again = runner.invoke(main, ["init"])
        assert again.exit_code == 1
        assert "already exists" in again.output

        forced = runner.invoke(main, ["init", "--full", "--force"])
        assert forced.exit_code == 0
        assert "demultiplex:" in Path("teleo.yaml").read_text()


def test_cli_init_rejects_template_argument() -> None:
    result = CliRunner().invoke(main, ["init", "complete"])
    assert result.exit_code != 0


@pytest.mark.parametrize(
    "filename, marker, template",
    [
        ("teleo.yaml", "teleo", "full"),
        ("mifish.yaml", "mifish", "full"),
        ("mam07.yaml", "mam07", "full"),
        ("minimal.example.yaml", "teleo", "minimal"),
    ],
)
def test_shipped_marker_configs_match_init(
    tmp_path: Path, filename: str, marker: str, template: str
) -> None:
    shipped = Path(__file__).resolve().parents[2] / "config" / "markers" / filename
    text, _ = _render(tmp_path, template, marker=marker)
    assert shipped.read_text(encoding="utf-8") == text, (
        f"config/markers/{filename} is out of date. Regenerate it with:\n"
        f"  seednap init --{template} -m {marker} -o config/markers/{filename} --force"
    )
