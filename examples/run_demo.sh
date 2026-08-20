#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
report_dir="${project_root}/reports"
mkdir -p "${report_dir}"

ngs-core validate \
  "${project_root}/tests/data/reads_R1.fastq" \
  --read2 "${project_root}/tests/data/reads_R2.fastq"

ngs-core qc \
  "${project_root}/tests/data/reads_R1.fastq" \
  --read2 "${project_root}/tests/data/reads_R2.fastq" \
  --sample synthetic-paired-demo \
  --output "${report_dir}/synthetic-paired-demo.html"

echo "Demo report: ${report_dir}/synthetic-paired-demo.html"
