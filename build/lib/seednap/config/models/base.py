"""Base Pydantic model for SeeDNAP config (strict: rejects unknown fields)."""

from pydantic import BaseModel, ConfigDict


class StrictModel(BaseModel):
    """Base model that rejects unknown fields to catch config typos."""

    model_config = ConfigDict(extra="forbid")
