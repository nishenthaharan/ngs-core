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

ngs-core filter \
  "${project_root}/tests/data/reads_R1.fastq" \
  --read2 "${project_root}/tests/data/reads_R2.fastq" \
  --output "${report_dir}/clean_R1.fastq.gz" \
  --output2 "${report_dir}/clean_R2.fastq.gz" \
  --min-length 1 \
  --max-n-fraction 1 \
  --quality-trim=-1 \
  --stats "${report_dir}/filtering.json"

ngs-core validate \
  "${report_dir}/clean_R1.fastq.gz" \
  --read2 "${report_dir}/clean_R2.fastq.gz"

echo "Demo report: ${report_dir}/synthetic-paired-demo.html"
echo "Filtering stats: ${report_dir}/filtering.json"
