# NGS Core

[![CI](https://github.com/nishenthaharan/ngs-core/actions/workflows/ci.yml/badge.svg)](https://github.com/nishenthaharan/ngs-core/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-3776AB.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

**NGS Core v0.1** provides a strict, streaming FASTQ foundation for reproducible
next-generation sequencing workflows. It validates single- and paired-end files without
loading complete datasets into memory.

## Capabilities

- Reads plain FASTQ and gzip-compressed FASTQ detected from file bytes.
- Validates headers, separators, record completeness, sequence/quality lengths, and bases.
- Detects paired files with different read counts or unsynchronised identifiers.
- Supports common `/1`, `/2`, and CASAVA-style read naming.
- Produces a machine-readable validation summary.
- Exposes typed Python iterators for use in future NGS modules.

## Installation

```bash
git clone https://github.com/nishenthaharan/ngs-core.git
cd ngs-core
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

NGS Core requires Python 3.10 or later and has no runtime dependencies.

## Command-line use

Validate one file:

```bash
ngs-core validate sample.fastq.gz
```

Validate paired-end reads and confirm that the mates are synchronised:

```bash
ngs-core validate sample_R1.fastq.gz --read2 sample_R2.fastq.gz
```

Example output:

```json
{
  "valid": true,
  "files": 2,
  "reads": 5000000,
  "bases": 750000000
}
```

Malformed records return exit code `2` with the affected file and record number.

## Python API

```python
from ngs_core import read_fastq, read_paired_fastq

for record in read_fastq("sample.fastq.gz"):
    print(record.identifier, len(record.sequence))

for read1, read2 in read_paired_fastq("R1.fastq.gz", "R2.fastq.gz"):
    assert read1.pair_key == read2.pair_key
```

## Why this is the first version

Every later QC, trimming, alignment, or visualisation component depends on correct FASTQ
iteration and pair integrity. Version 0.1 therefore establishes and tests that boundary
before analytical features are added.

## Performance model

Parsing is `O(n)` in the number of bases. Memory usage is normally `O(r)`, where `r` is
the current read length. No full-file record list is created by the command-line tool.

## Development

```bash
python -m pip install -e ".[dev]"
ruff check .
ruff format --check .
pytest
python -m build
```

Synthetic test fixtures cover plain and gzip input, truncated records, length mismatches,
paired-name mismatches, unequal read counts, CLI summaries, and error handling.

## Next version

Version 0.2 adds streaming QC metrics and JSON, TSV, and self-contained HTML reports on
top of this validated FASTQ layer.

## License

Copyright © 2026 Nishenthaharan Balachandran. Licensed under the
[Apache License 2.0](LICENSE).
