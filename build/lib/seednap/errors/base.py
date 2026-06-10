"""Structured, user-facing errors for SeeDNAP.

The design follows the consensus of the rustc / Elm compilers and Google's technical-writing
guidance (and the no-silent-fallbacks policy): a good error states WHAT failed, WHY, and HOW to fix it,
never fails silently, and surfaces the root cause. Rather than ad-hoc strings, an error carries
those parts as fields plus an optional stable code and docs pointer, and renders them
consistently. The original exception is preserved via ``raise ... from`` so ``-v`` can still
show a developer traceback.
"""

from typing import Optional


class SeednapError(Exception):
    """A user-facing error that explains itself.

    Args:
        summary: One line: what went wrong (no trailing period needed).
        why: Why it happened / the root cause, in plain language.
        fix: The concrete, declarative remedy (an instruction, not a question).
        code: Optional stable error code (e.g. ``SDN-CFG-001``) for `seednap explain`.
        docs: Optional pointer to a docs section.
    """

    def __init__(
        self,
        summary: str,
        *,
        why: Optional[str] = None,
        fix: Optional[str] = None,
        code: Optional[str] = None,
        docs: Optional[str] = None,
    ) -> None:
        self.summary = summary
        self.why = why
        self.fix = fix
        self.code = code
        self.docs = docs
        super().__init__(self.render())

    def render(self) -> str:
        """Render the what / why / fix triad as a single multi-line message."""
        lines = [self.summary]
        if self.why:
            lines.append(f"  Why: {self.why}")
        if self.fix:
            lines.append(f"  Fix: {self.fix}")
        if self.docs:
            lines.append(f"  See: {self.docs}")
        if self.code:
            lines.append(f"  [{self.code}]  (run `seednap explain {self.code}` for more)")
        return "\n".join(lines)
