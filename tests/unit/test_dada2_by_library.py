"""Unit tests for the DADA2-by-library wiring (item C, Python side).

The R per-library branch and the no-op byte-identical guarantee are validated by end-to-end
DADA2 runs; here we pin the config and the runner's arg threading (the library_map must be
the 15th positional arg, matching scripts/dada2_process.R's args[15]).
"""

from pathlib import Path

from seednap.config.models import Dada2Config
from seednap.steps.dada2.dada2_runner import Dada2Runner


def test_per_library_config_default_off():
    assert Dada2Config().per_library is False
    assert Dada2Config(per_library=True).per_library is True


def test_runner_threads_library_map_as_arg15(monkeypatch, tmp_path):
    captured = {}

    def fake_run(self, script_path, args, log_file=None):
        captured["args"] = args
        return ""

    monkeypatch.setattr(Dada2Runner, "_run_r_script", fake_run)
    runner = Dada2Runner()

    # With a library map: it must be the 15th positional arg.
    runner.run_dada2_process(
        marker="m", trimmed_reads_dir=tmp_path, output_dir=tmp_path,
        library_map=tmp_path / "library_map.csv",
    )
    args = captured["args"]
    assert len(args) == 15, f"expected 15 args, got {len(args)}: {args}"
    assert args[14] == str(tmp_path / "library_map.csv")

    # Without a library map: arg 15 is empty (R falls back to the single-batch path).
    runner.run_dada2_process(marker="m", trimmed_reads_dir=tmp_path, output_dir=tmp_path)
    assert captured["args"][14] == ""
