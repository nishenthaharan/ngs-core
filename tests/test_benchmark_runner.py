from __future__ import annotations

import sys
import unittest
from pathlib import Path

from benchmarks.run_benchmarks import _build_command, measure_command


class BenchmarkRunnerTests(unittest.TestCase):
    def test_measures_successful_process(self) -> None:
        result = measure_command(
            [sys.executable, "-c", "print('benchmark-probe')"],
            operation="probe",
            profile="test",
            repetition=1,
            reads=10,
            bases=100,
            poll_interval=0.001,
        )
        self.assertEqual(result["return_code"], 0)
        self.assertGreater(result["wall_seconds"], 0)
        self.assertGreater(result["reads_per_second"], 0)
        self.assertEqual(len(result["stdout_sha256"]), 64)

    def test_records_failed_process_stderr(self) -> None:
        result = measure_command(
            [sys.executable, "-c", "import sys; print('failed', file=sys.stderr); sys.exit(3)"],
            operation="probe",
            profile="test",
            repetition=1,
            reads=1,
            bases=1,
            poll_interval=0.001,
        )
        self.assertEqual(result["return_code"], 3)
        self.assertIn("failed", result["stderr"])

    def test_builds_paired_filter_command(self) -> None:
        command = _build_command(
            "ngs-core",
            "filter",
            Path("R1.fastq"),
            Path("R2.fastq"),
            Path("results"),
            2,
        )
        self.assertIn("--read2", command)
        self.assertIn("--output2", command)
        self.assertIn("filtered-R1-2.fastq", " ".join(command))


if __name__ == "__main__":
    unittest.main()
