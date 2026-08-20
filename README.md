# NGS Core

[![CI](https://github.com/nishenthaharan/ngs-core/actions/workflows/ci.yml/badge.svg)](https://github.com/nishenthaharan/ngs-core/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-3776AB.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

**NGS Core** is a dependency-free Python toolkit for validating, assessing, and
preprocessing FASTQ data. It is designed as the reliable foundation for a broader
collection of command-line NGS analysis and visualisation tools.

The implementation is streaming: reads are processed one record at a time, so
routine operations remain memory-efficient on multi-gigabyte sequencing files.

## What it does

- Reads plain or gzip-compressed FASTQ, detected from the file content.
- Validates four-line FASTQ structure, sequence/quality lengths, and quality encoding.
- Detects unsynchronised paired-end reads before downstream analysis.
- Calculates read count, bases, GC%, N%, Q20, Q30, length distribution, read-length
  N50, estimated duplication, per-cycle quality, base composition, and adapter hits.
- Produces reproducible JSON, analysis-friendly TSV, or self-contained HTML reports.
- Performs fixed trimming, exact adapter clipping, 3′ sliding-window quality trimming,
  and length/quality/ambiguous-base filtering.
- Preserves paired-end synchronisation by discarding both mates when either fails.
- Runs without uploading sequence data or requiring a third-party service.

## Installation

NGS Core requires Python 3.10 or later.

```bash
git clone https://github.com/nishenthaharan/ngs-core.git
cd ngs-core
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

For development tools:

```bash
python -m pip install -e ".[dev]"
```

## Quick start

Validate a single-end file:

```bash
ngs-core validate sample.fastq.gz
```

Validate paired-end structure and read-name synchronisation:

```bash
ngs-core validate sample_R1.fastq.gz --read2 sample_R2.fastq.gz
```

Generate a self-contained QC report:

```bash
ngs-core qc sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --sample tumour-replicate-1 \
  --output reports/tumour-replicate-1.html
```

Write machine-readable metrics instead:

```bash
ngs-core qc sample.fastq.gz --output reports/sample.json
ngs-core qc sample.fastq.gz --output reports/sample.tsv
```

Filter paired-end reads while keeping the mates synchronized:

```bash
ngs-core filter sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --output clean_R1.fastq.gz \
  --output2 clean_R2.fastq.gz \
  --min-length 50 \
  --min-mean-quality 20 \
  --quality-trim 20 \
  --quality-window 4 \
  --max-n-fraction 0.05 \
  --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --stats reports/filtering.json
```

Use `ngs-core <command> --help` for every option.

## Commands

| Command | Purpose | Principal output |
|---|---|---|
| `validate` | Validate FASTQ records and optional R1/R2 synchronisation | JSON summary |
| `qc` | Calculate streaming sequence and quality metrics | JSON, TSV, or HTML |
| `filter` | Trim and remove failing reads or pairs | FASTQ/FASTQ.GZ and JSON stats |

## Quality-control schema

JSON reports declare `schema_version: "1.0"`. Each mate has the same metric fields,
including a `per_cycle` array and `length_distribution` array. This stable structure is
intended for workflow engines, dashboards, notebooks, and future NGS Core modules.

```json
{
  "schema_version": "1.0",
  "tool": "ngs-core",
  "tool_version": "0.1.0",
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

The example is abbreviated; generated reports include the complete schema.

## Python API

The parsing and QC components can also be embedded in a Python workflow:

```python
from ngs_core import calculate_qc, read_fastq

result = calculate_qc("sample.fastq.gz")
print(result.q30_percent)

for record in read_fastq("sample.fastq.gz"):
    print(record.identifier, len(record.sequence))
```

## Design principles

1. **Streaming by default** — memory use depends mainly on read length, not file size.
2. **Fail clearly** — malformed records and broken pairs stop with record-level diagnostics.
3. **Reproducible outputs** — explicit parameters and versioned machine-readable reports.
4. **Composable components** — CLI functions are also typed Python APIs.
5. **Local-first privacy** — raw reads stay in the user's analysis environment.

See [Architecture](docs/ARCHITECTURE.md) for module boundaries and implementation notes.

## Scope and scientific use

NGS Core v0.1 focuses on FASTQ integrity, descriptive QC, and deterministic
preprocessing. It does not replace assay-specific validation, FastQC/MultiQC comparisons,
alignment, contamination screening, or clinical-grade pipeline verification. Parameters
must be selected for the sequencing platform, library preparation, read length, and assay.

## Roadmap

- Multi-file sample manifests and batch QC aggregation
- Approximate adapter matching and poly-G/poly-X trimming
- K-mer spectra and contamination signatures
- FASTA support and sequence-format auto-detection
- MultiQC-compatible exports
- Workflow wrappers for Nextflow and Snakemake
- Alignment and coverage modules in separate focused repositories

## Development

```bash
ruff check .
ruff format --check .
pytest
python -m build
```

The tests include malformed FASTQ, gzip detection, paired-end synchronisation, exact QC
metrics, filtering behaviour, report rendering, and CLI error handling. See
[Contributing](CONTRIBUTING.md) before opening a change.

## Citation

If NGS Core supports published work, cite the software using the metadata in
[`CITATION.cff`](CITATION.cff) and record the exact release and command parameters.

## License

Copyright © 2026 Nishenthaharan Balachandran. Licensed under the
[Apache License 2.0](LICENSE).
