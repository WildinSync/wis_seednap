"""Regression tests for two standards-fix behaviors.

1. Quiet mode (-q) must still surface WARNING/ERROR on the console -- the no-silent-fallback
   policy means data-loss [WARN]s stay visible even when INFO chatter is suppressed.
2. Run state stamps the seednap version on new runs; a state file written before version
   stamping existed (no field) loads as None so resume can flag unknown provenance.
"""

import json
import logging

from rich.logging import RichHandler

from seednap.__version__ import __version__
from seednap.pipeline.state import PipelineState
from seednap.utils.logging import setup_logging


def test_quiet_mode_keeps_warning_console_handler(tmp_path):
    setup_logging(level="WARNING", log_file=tmp_path / "q.log", console_output=False)
    root = logging.getLogger()
    rich = [h for h in root.handlers if isinstance(h, RichHandler)]
    assert rich, "quiet mode must still attach a console handler so warnings/errors stay visible"
    assert all(h.level <= logging.WARNING for h in rich)
    setup_logging(level="INFO", console_output=True)  # restore default config


def test_state_stamps_version_on_new_runs(tmp_path):
    assert PipelineState.from_config(marker="t").seednap_version == __version__


def test_legacy_state_without_version_loads_as_none(tmp_path):
    s = PipelineState.from_config(marker="t")
    d = json.loads(s.model_dump_json())
    d.pop("seednap_version", None)  # simulate a state file written before version stamping
    p = tmp_path / "state.json"
    p.write_text(json.dumps(d))
    assert PipelineState.load(p).seednap_version is None
