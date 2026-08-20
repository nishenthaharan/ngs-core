# Architecture

NGS Core separates FASTQ I/O, analysis, transformation, presentation, and command-line
concerns so each component can be tested or embedded independently.

## Data flow

```text
FASTQ / FASTQ.GZ
      │
      ▼
fastq.py ── validation and paired-read synchronization
      │
      ├──► qc.py ── streaming accumulators ──► report.py ──► JSON / TSV / HTML
      │
      └──► filtering.py ── trimming and predicates ──► FASTQ / FASTQ.GZ
```

## Modules

### `fastq.py`

Owns the `FastqRecord` data model, content-based gzip detection, FASTQ parsing, Phred
decoding, and paired-end iteration. It is the only module that interprets FASTQ layout.
The parser is deliberately strict and rejects wrapped records because silently accepting
ambiguous structure can desynchronise downstream reads.

### `qc.py`

Uses `QCAccumulator` to update exact counts for every read and base. Per-cycle arrays grow
only to the longest observed read. Read-length storage is a histogram rather than one
entry per record. Duplication is explicitly an estimate based on at most the first 100,000
read prefixes, bounding memory use while keeping the metric useful for screening.

### `filtering.py`

Applies operations in a stable order:

1. fixed 5′/3′ trimming;
2. exact adapter clipping;
3. 3′ sliding-window quality trimming;
4. minimum length;
5. maximum N fraction;
6. minimum mean quality.

Paired-end filtering uses an all-or-none policy. If either mate fails, neither is written.

### `report.py`

Serializes the stable report payload. HTML output has no JavaScript, web font, CDN, or
external asset dependency; it can be opened offline and archived with an analysis run.

### `cli.py`

Defines user-facing commands, validates option combinations, and translates expected
domain failures into exit code `2`. Unexpected programming failures are not hidden.

## Complexity

For `n` total bases and maximum read length `r`:

- Parsing, QC, and filtering time: `O(n)`
- Core per-cycle memory: `O(r)`
- Length histogram memory: `O(number of distinct read lengths)`
- Duplication estimate: bounded to 100,000 sampled prefixes

## Extension rules

- Keep FASTQ parsing and validation centralized in `fastq.py`.
- Add new metrics to `QCResult` and every report format together.
- Increment `schema_version` for breaking report changes.
- Add representative unit tests before exposing a new CLI option.
- Keep the core dependency-free; optional integrations should be isolated extras or
  separate repositories.
