# Architecture — v0.1

Version 0.1 intentionally establishes one narrow boundary: trusted streaming FASTQ input.

```text
FASTQ / FASTQ.GZ
      │
      ▼
open_fastq_text()
      │
      ▼
read_fastq() ──► FastqRecord
      │
      └──► read_paired_fastq() ──► synchronized R1/R2 tuples
```

## `fastq.py`

- Detects gzip using magic bytes instead of relying only on a filename suffix.
- Reads exactly four lines for each record and reports the failing record number.
- Normalises sequence bases to uppercase while preserving quality characters.
- Validates sequence/quality length equality and decodes Phred scores explicitly.
- Iterates paired files together and stops on count or identifier mismatches.

The parser rejects wrapped FASTQ. Silent acceptance of ambiguous wrapping can shift record
boundaries and corrupt paired-read interpretation.

## `cli.py`

The `validate` command consumes the same public iterators used by the Python API. Expected
data and filesystem failures become concise diagnostics with exit code `2`; unexpected
programming failures are not hidden.

## Extension contract

Later modules must consume `FastqRecord` and the public readers rather than implementing
new FASTQ parsing. This keeps validation behaviour consistent across QC and preprocessing.
