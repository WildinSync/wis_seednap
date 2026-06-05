"""Control decontamination (cleaning) for SeeDNAP abundance tables.

Removes (or flags) OTU/ASV reads that appear in negative controls, using the FAIRe
manifest's control taxonomy and extraction-batch associations. See
:mod:`seednap.steps.cleaning.processor`.
"""

from seednap.steps.cleaning.processor import CleaningProcessor, CleaningResult

__all__ = ["CleaningProcessor", "CleaningResult"]
