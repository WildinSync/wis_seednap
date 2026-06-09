"""Unrunnable demultiplexing protocols must fail at config load, not mid-run.

Only 'ligation' demultiplexing is implemented. If 'demultiplex' is in pipeline.steps, the
protocol must be 'ligation'; anything else ('standard' -> NotImplementedError, or the default
'none' -> "Unknown protocol") would otherwise crash mid-run after trimming has already run.
"""

import pytest
from pydantic import ValidationError

from seednap.config.models import PipelineConfig


def _cfg(tmp_path, steps, protocol):
    return {
        "marker": {
            "name": "t",
            "primers": {"forward": "ACACCGCCCGTCACTCT", "reverse": "CTTCCGGTACACTTACCATG"},
        },
        "paths": {
            "raw_data": str(tmp_path),
            "output": str(tmp_path / "o"),
            "logs": str(tmp_path / "o" / "l"),
        },
        "taxonomy": {"method": "blast", "databases": {"blast": {"fasta": str(tmp_path / "r.fasta")}}},
        "demultiplex": {"protocol": protocol},
        "pipeline": {"steps": steps},
    }


@pytest.mark.parametrize("protocol", ["standard", "none"])
def test_non_ligation_protocol_rejected_when_demux_runs(tmp_path, protocol):
    with pytest.raises(ValidationError) as exc:
        PipelineConfig(**_cfg(tmp_path, ["demultiplex", "trim", "swarm", "taxonomy"], protocol))
    msg = str(exc.value)
    assert "ligation" in msg and "demultiplex.protocol" in msg


def test_ligation_protocol_with_demux_is_accepted(tmp_path):
    cfg = PipelineConfig(**_cfg(tmp_path, ["demultiplex", "trim", "swarm", "taxonomy"], "ligation"))
    assert cfg.demultiplex.protocol == "ligation"


@pytest.mark.parametrize("protocol", ["standard", "none"])
def test_protocol_ignored_when_demux_not_in_steps(tmp_path, protocol):
    # protocol is irrelevant if demultiplexing does not run; do not error on a leftover value.
    cfg = PipelineConfig(**_cfg(tmp_path, ["trim", "swarm", "taxonomy"], protocol))
    assert cfg.demultiplex.protocol == protocol
