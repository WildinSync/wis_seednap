# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Error-explainability module with a `seednap explain` command: errors carry
  stable codes and actionable what / why / how-to-fix detail, and the codes can
  be looked up from the CLI.

### Changed

- Pipeline stage enable/disable now flows through a single `pipeline.steps`
  config model with dependency validation, replacing the previous scattered
  per-stage toggles.
- Standardized the documentation for accuracy and a consistent structure so the
  docs match the implementation.

### Fixed

- Correctness sweep across the pipeline focused on data integrity, removing
  silent fallbacks (fallbacks now warn or fail loudly), and catching
  wrong-environment misconfiguration earlier.

## [0.1.0]

- Initial alpha release of the SeeDNAP eDNA metabarcoding pipeline.
