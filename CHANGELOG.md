# Changelog

All notable changes to NGS Core are documented here.

The project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.1.0] - 2026-08-20

### Added

- Streaming FASTQ and FASTQ.GZ parsing with structural validation.
- Paired-end name and record-count synchronization checks.
- QC metrics covering bases, GC, N, Q20/Q30, length, duplication estimate, cycles,
  composition, and common adapter signatures.
- JSON, TSV, and self-contained HTML QC reports.
- Fixed, adapter, and sliding-window quality trimming.
- Single- and paired-end filtering with gzip output.
- Typed Python API, command-line interface, tests, CI, container definition, and docs.

[Unreleased]: https://github.com/nishenthaharan/ngs-core/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/nishenthaharan/ngs-core/releases/tag/v0.1.0
