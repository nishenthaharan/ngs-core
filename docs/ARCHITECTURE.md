# Architecture — v0.2

Version 0.2 retains the v0.1 streaming parser and adds bounded-memory QC accumulation and
three report renderers.

```text
FASTQ / FASTQ.GZ
      │
      ▼
fastq.py ──► validated FastqRecord stream
      │
      ▼
qc.py ──► QCAccumulator ──► QCResult
                              │
                              ▼
                    report.py ──► JSON / TSV / HTML
```

## QC accumulation

Every read updates exact base, quality, length, and per-cycle counts. The length
distribution is stored as a histogram rather than one element per read. Per-cycle arrays
grow only to the longest read observed.

Duplication estimation stores at most the first 100,000 sequence prefixes, explicitly
bounding that component. Adapter screening counts exact signatures and does not claim
approximate alignment.

## Reporting

JSON is the stable machine-readable representation. TSV flattens summary metrics for
statistical tools. HTML embeds its CSS and SVG chart, requires no JavaScript or network
connection, and can be archived with an analysis run.

## Safety boundary

The CLI resolves paths before writing a report and rejects a destination that aliases an
input FASTQ. The analytical result is fully rendered in memory before the destination is
written.

## Complexity

For `n` bases and maximum read length `r`, QC time is `O(n)` and core cycle memory is
`O(r)`. The number of distinct read lengths and bounded duplication sample add small,
explicit state.
