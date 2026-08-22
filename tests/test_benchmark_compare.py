from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from benchmarks.compare_results import compare_payloads, main, summarize


def _payload(wall_times: list[float], throughputs: list[float]) -> dict[str, object]:
    return {
        "observations": [
            {
                "profile": "smoke",
                "operation": "validate",
                "return_code": 0,
                "wall_seconds": wall,
                "reads_per_second": throughput,
                "bases_per_second": throughput * 150,
            }
            for wall, throughput in zip(wall_times, throughputs, strict=True)
        ]
    }


class BenchmarkComparisonTests(unittest.TestCase):
    def test_summarizes_medians(self) -> None:
        summary = summarize(_payload([3.0, 1.0, 2.0], [30.0, 10.0, 20.0]))
        row = summary[("smoke", "validate")]
        self.assertEqual(row["observations"], 3)
        self.assertEqual(row["median_wall_seconds"], 2.0)
        self.assertEqual(row["median_reads_per_second"], 20.0)

    def test_compares_matching_groups(self) -> None:
        rows = compare_payloads(_payload([2.0], [50.0]), _payload([2.2], [45.0]))
        self.assertEqual(rows[0]["wall_regression_percent"], 10.0)
        self.assertEqual(rows[0]["throughput_change_percent"], -10.0)

    def test_cli_returns_failure_above_threshold(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline = root / "baseline.json"
            candidate = root / "candidate.json"
            output = root / "comparison.json"
            baseline.write_text(json.dumps(_payload([1.0], [100.0])), encoding="utf-8")
            candidate.write_text(json.dumps(_payload([1.2], [80.0])), encoding="utf-8")
            exit_code = main(
                [
                    "--baseline",
                    str(baseline),
                    "--candidate",
                    str(candidate),
                    "--output",
                    str(output),
                    "--max-wall-regression-percent",
                    "10",
                ]
            )
            self.assertEqual(exit_code, 1)
            self.assertEqual(len(json.loads(output.read_text())["comparisons"]), 1)


if __name__ == "__main__":
    unittest.main()
