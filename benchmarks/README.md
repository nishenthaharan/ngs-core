# NGS Core benchmarks

These tools generate synthetic FASTQ inputs, measure installed NGS Core commands, and
compare median results between revisions. They implement the protocol in
[`docs/BENCHMARKING.md`](../docs/BENCHMARKING.md) and track
[issue #5](https://github.com/nishenthaharan/ngs-core/issues/5).

Benchmark inputs contain no biological samples. The generator records its seed and file
checksums so a result can be audited or reproduced.

## 1. Generate an input

```bash
python -m benchmarks.generate_fastq \
  --output /tmp/ngs-core-R1.fastq.gz \
  --output2 /tmp/ngs-core-R2.fastq.gz \
  --records 10000 \
  --read-length 150 \
  --seed 20260821 \
  --manifest /tmp/ngs-core-manifest.json
```

Use `profiles.json` for the versioned smoke, standard, paired, and variable-length workload
definitions. Large profiles are generated locally and are never committed.

## 2. Measure a command

Install NGS Core first, then run the harness against its command-line entry point:

```bash
python -m benchmarks.run_benchmarks \
  --operation validate \
  --input /tmp/ngs-core-R1.fastq.gz \
  --read2 /tmp/ngs-core-R2.fastq.gz \
  --profile paired \
  --records 10000 \
  --bases 3000000 \
  --repetitions 5 \
  --output /tmp/baseline.json
```

The harness records wall and child CPU time, throughput, command arguments, revision,
environment, exit status, and Linux peak RSS when `/proc` is available. Shared CI runners
are useful only for smoke testing and must not be published as stable performance baselines.

## 3. Compare revisions

Produce a candidate file using the same machine, dataset, arguments, and cache state, then
compare matching profile/operation groups:

```bash
python -m benchmarks.compare_results \
  --baseline /tmp/baseline.json \
  --candidate /tmp/candidate.json \
  --output /tmp/comparison.json \
  --max-wall-regression-percent 10
```

Exit code `1` means at least one median wall-time regression exceeded the threshold. Exit
code `2` means the files had no comparable successful groups. Always inspect the complete
observation distribution before treating a threshold breach as a product regression.
