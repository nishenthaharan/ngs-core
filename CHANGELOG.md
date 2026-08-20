# Changelog

All notable changes to NGS Core are documented here. The project follows Semantic
Versioning.

## [Unreleased]

## [0.3.0] - 2026-08-20

### Added

- Fixed 5′/3′ trimming and exact adapter clipping.
- Sliding-window 3′ quality trimming.
- Minimum length, mean quality, and ambiguous-base filtering.
- Synchronized paired-end retention with gzip output.
- Atomic FASTQ output files and machine-readable filtering statistics.
- Collision protection for FASTQ, report, filtered-read, and statistics paths.

## [0.2.0] - 2026-08-20

### Added

- Streaming base, quality, length, cycle, duplication, and adapter metrics.
- JSON schema 1.0, TSV export, and self-contained HTML QC reports.
- Single- and paired-end QC with optional read sampling.

## [0.1.0] - 2026-08-20

### Added

- Streaming FASTQ and FASTQ.GZ parsing.
- Structural and paired-end synchronization validation.
- JSON validation summaries and typed Python iterators.

[Unreleased]: https://github.com/nishenthaharan/ngs-core/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/nishenthaharan/ngs-core/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/nishenthaharan/ngs-core/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/nishenthaharan/ngs-core/releases/tag/v0.1.0
