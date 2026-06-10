"""Version information for seednap."""

__version__ = "0.1.0"
# isdigit() drops any non-numeric segment, so a pre-release like "1.2.0-rc1"
# yields (1, 2): the "0-rc1" patch segment is skipped, truncating the tuple.
__version_info__ = tuple(int(i) for i in __version__.split(".") if i.isdigit())