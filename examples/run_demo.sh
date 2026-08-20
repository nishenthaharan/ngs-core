#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

ngs-core validate \
  "${project_root}/tests/data/reads_R1.fastq" \
  --read2 "${project_root}/tests/data/reads_R2.fastq"
