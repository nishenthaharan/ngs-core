# Benchmarking NGS Core

This document defines the protocol for measuring NGS Core performance. It is the first
deliverable for [issue #5](https://github.com/nishenthaharan/ngs-core/issues/5); benchmark
automation and published baseline results will be added separately.

## Goals

The benchmark suite should answer four questions:

1. How many reads and bases can each command process per second?
2. Does peak memory remain bounded as input files grow?
3. What overhead is introduced by gzip decompression and paired-end validation?
4. Do later changes cause statistically meaningful performance regressions?

The protocol measures engineering performance only. It does not establish biological or
clinical validity.

## Dataset profiles

All generated reads must be synthetic, deterministic, and free of identifiable sequence
data. The generator seed and every generation parameter must be stored with the results.

| Profile | Records | Read length | Layout | Compression | Purpose |
|---|---:|---:|---|---|---|
| Smoke | 10,000 | 150 bp | Single-end | Plain and gzip | Fast command and CI checks |
| Standard | 1,000,000 | 150 bp | Single-end | Plain and gzip | Routine throughput comparison |
| Paired | 1,000,000 pairs | 150 bp | Paired-end | Plain and gzip | Synchronisation overhead |
| Variable | 1,000,000 | 50–250 bp | Single-end | gzip | Per-cycle array growth and mixed lengths |
| Filtering | 1,000,000 pairs | 150 bp | Paired-end | gzip | Trimming, rejection, and output compression |

Only the smoke fixture belongs in the repository. Larger inputs should be generated in a
temporary directory and removed after measurement.

## Commands under test

Run the installed command-line entry point rather than importing private functions:

```bash
ngs-core validate reads.fastq.gz
ngs-core qc reads.fastq.gz --output qc.json
ngs-core filter reads.fastq.gz \
  --output filtered.fastq.gz \
  --stats filtering.json \
  --min-length 50 \
  --min-mean-quality 20 \
  --quality-trim 20 \
  --quality-window 4
```

Paired profiles must add `--read2` and `--output2`. Benchmark outputs should be written to
the same filesystem used for the input unless the storage experiment explicitly varies
that condition.

## Measurements

Record these values for every run:

- wall-clock and CPU time
- peak resident memory
- reads and bases processed per second
- input and output byte counts
- command, complete arguments, exit code, and NGS Core version
- Python version, operating system, CPU model and count, available memory, and storage type
- gzip implementation and compression level where applicable
- dataset profile, generator seed, and input checksum

Machine-readable results should use JSON or TSV. Raw timing output must remain available so
summary calculations can be audited.

## Execution procedure

1. Use an otherwise idle machine with a fixed power profile.
2. Build a clean virtual environment and install the exact commit under test.
3. Generate inputs once and verify their checksums.
4. Run one unmeasured warm-up for each command and profile.
5. Run five measured repetitions in alternating baseline/candidate order.
6. Report the median and range; retain every individual observation.
7. Compare revisions only on the same machine and storage path.

File-system cache state can materially change results. Reports must state whether they
measure warm-cache or cold-cache performance; the two modes must never be combined.

## Result record

Each observation should contain at least:

```json
{
  "schema_version": "1.0",
  "revision": "git-commit-sha",
  "command": "validate",
  "profile": "standard-gzip",
  "repetition": 1,
  "reads": 1000000,
  "bases": 150000000,
  "wall_seconds": 0.0,
  "cpu_seconds": 0.0,
  "peak_rss_bytes": 0,
  "input_bytes": 0,
  "environment": {
    "python": "3.x",
    "operating_system": "recorded by runner",
    "cpu": "recorded by runner",
    "storage": "recorded by runner"
  }
}
```

Zero values are placeholders in this schema example, not measured performance claims.

## Regression interpretation

Do not compare isolated best runs. Compare medians from the same environment and inspect
the complete distribution. A regression threshold should be introduced only after enough
baseline runs exist to quantify normal variance. Correctness tests must pass before any
performance result is accepted.

## Reporting checklist

- [ ] Exact commit and command arguments recorded
- [ ] Dataset seed, parameters, and checksums recorded
- [ ] Hardware and software environment recorded
- [ ] Warm-up and five measured repetitions completed
- [ ] Median, range, throughput, and peak memory reported
- [ ] Raw observations attached or committed when reasonably small
- [ ] No identifiable or licensed sequence dataset redistributed
