"""Base Pydantic model for SeeDNAP config (strict: rejects unknown fields)."""

from pydantic import BaseModel, ConfigDict


class StrictModel(BaseModel):
    """Base config model that rejects unknown fields (``extra="forbid"``).

    Every SeeDNAP config model inherits from this so an unrecognised YAML key (typically a
    typo, e.g. ``primers.forwrd``) errors at load time rather than being silently ignored and
    producing wrong-but-valid-looking behaviour on real biodiversity samples.
    """

    model_config = ConfigDict(extra="forbid")
