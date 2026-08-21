from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from benchmarks.generate_fastq import generate_fastq, main
from ngs_core.fastq import read_fastq, read_paired_fastq


class BenchmarkGeneratorTests(unittest.TestCase):
    def test_plain_output_is_byte_reproducible(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first.fastq"
            second = Path(directory) / "second.fastq"
            generate_fastq(first, records=12, read_length=25, seed=7)
            generate_fastq(second, records=12, read_length=25, seed=7)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(len(list(read_fastq(first))), 12)

    def test_gzip_output_is_byte_reproducible(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first.fastq.gz"
            second = Path(directory) / "second.fastq.gz"
            generate_fastq(first, records=8, read_length=20, seed=11)
            generate_fastq(second, records=8, read_length=20, seed=11)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_paired_outputs_have_matching_identifiers(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            read1 = Path(directory) / "R1.fastq"
            read2 = Path(directory) / "R2.fastq"
            generate_fastq(read1, records=10, read_length=30, seed=19, mate=1)
            generate_fastq(read2, records=10, read_length=30, seed=19, mate=2)
            self.assertEqual(len(list(read_paired_fastq(read1, read2))), 10)

    def test_cli_writes_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "reads.fastq"
            manifest = Path(directory) / "manifest.json"
            exit_code = main(
                [
                    "--output",
                    str(output),
                    "--records",
                    "5",
                    "--read-length",
                    "18",
                    "--manifest",
                    str(manifest),
                ]
            )
            payload = json.loads(manifest.read_text(encoding="utf-8"))
            self.assertEqual(exit_code, 0)
            self.assertEqual(payload["files"][0]["records"], 5)
            self.assertEqual(len(payload["files"][0]["sha256"]), 64)

    def test_rejects_invalid_lengths(self) -> None:
        with (
            tempfile.TemporaryDirectory() as directory,
            self.assertRaisesRegex(ValueError, "minimum"),
        ):
            generate_fastq(
                Path(directory) / "reads.fastq",
                records=1,
                read_length=10,
                min_read_length=11,
                seed=1,
            )


if __name__ == "__main__":
    unittest.main()
