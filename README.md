# NGS Core

[![CI](https://github.com/nishenthaharan/ngs-core/actions/workflows/ci.yml/badge.svg)](https://github.com/nishenthaharan/ngs-core/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-3776AB.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

**NGS Core v0.2** adds memory-efficient FASTQ quality control and visual reporting to the
validated single- and paired-end parsing foundation introduced in v0.1.

## Capabilities

- Plain and gzip FASTQ parsing with structural and paired-end validation.
- Read count, total bases, GC%, N%, Q20, Q30, and mean quality.
- Minimum, maximum, mean, median, distribution, and read-length N50.
- Per-cycle mean quality, Q20/Q30, and A/C/G/T/N composition.
- Bounded duplication estimation using up to 100,000 read prefixes.
- Exact signatures for common Illumina, Nextera, and small-RNA adapters.
- JSON, TSV, and self-contained offline HTML reports.
- Protection against accidentally overwriting FASTQ inputs with reports.

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

## Validate input

```bash
ngs-core validate sample_R1.fastq.gz --read2 sample_R2.fastq.gz
```

## Generate QC reports

Create a self-contained HTML report:

```bash
ngs-core qc sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --sample tumour-replicate-1 \
  --output reports/tumour-replicate-1.html
```

Create analysis-friendly formats:

```bash
ngs-core qc sample.fastq.gz --output reports/sample.json
ngs-core qc sample.fastq.gz --output reports/sample.tsv
```

Sample only the first one million reads when rapid screening is appropriate:

```bash
ngs-core qc sample.fastq.gz --max-reads 1000000 --output quick-qc.json
```

## Stable report schema

Reports declare `schema_version: "1.0"`. Paired input produces `read1` and `read2`
objects with identical fields, allowing dashboards and workflow engines to consume either
layout predictably.

```json
{
  "schema_version": "1.0",
  "tool": "ngs-core",
  "tool_version": "0.2.0",
  "sample": "tumour-replicate-1",
  "results": {
    "read1": {
      "reads": 1250000,
      "q30_percent": 91.72,
      "gc_percent": 48.13
    }
  }
}
```

The example is abbreviated; generated reports include cycle and length arrays.

## Python API

```python
from ngs_core import calculate_qc

result = calculate_qc("sample.fastq.gz")
print(result.q30_percent, result.gc_percent)
```

## Metric interpretation

The duplication value is a screening estimate, not an optical-duplicate or UMI-aware
calculation. Adapter observations use exact sequence signatures. Assay-specific decisions
should consider the platform, library preparation, expected insert size, and downstream
analysis rather than applying universal thresholds.

## Development

```bash
python -m pip install -e ".[dev]"
ruff check .
ruff format --check .
pytest
python -m build
```

## Next version

Version 0.3 adds deterministic trimming, filtering, gzip output, paired-read retention,
atomic output handling, and machine-readable preprocessing statistics.

## License

Copyright © 2026 Nishenthaharan Balachandran. Licensed under the
[Apache License 2.0](LICENSE).
