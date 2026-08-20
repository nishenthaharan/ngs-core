from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ngs_core.fastq import FastqRecord
from ngs_core.qc import QCAccumulator, calculate_qc
from ngs_core.report import infer_report_format, render_report


DATA = Path(__file__).parent / "data"


class QCTests(unittest.TestCase):
    def test_calculates_expected_summary(self) -> None:
        result = calculate_qc(DATA / "reads_R1.fastq")
        self.assertEqual(result.reads, 3)
        self.assertEqual(result.bases, 30)
        self.assertEqual(result.min_length, 10)
        self.assertEqual(result.max_length, 10)
        self.assertEqual(result.read_length_n50, 10)
        self.assertEqual(result.median_length, 10.0)
        self.assertAlmostEqual(result.gc_percent, 40.0)
        self.assertAlmostEqual(result.n_percent, 6.6667)
        self.assertAlmostEqual(result.q20_percent, 66.6667)
        self.assertAlmostEqual(result.q30_percent, 33.3333)
        self.assertEqual(len(result.per_cycle), 10)
        self.assertAlmostEqual(result.per_cycle[0]["mean_quality"], 20.0)

    def test_handles_empty_fastq(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "empty.fastq"
            path.write_text("", encoding="utf-8")
            result = calculate_qc(path)
            self.assertEqual(result.reads, 0)
            self.assertEqual(result.mean_quality, 0.0)
            self.assertEqual(result.median_length, 0.0)

    def test_counts_adapter_signatures(self) -> None:
        accumulator = QCAccumulator("memory")
        sequence = "ACGTAGATCGGAAGAGTT"
        accumulator.add(FastqRecord("@r1", sequence, "+", "I" * len(sequence)))
        self.assertEqual(accumulator.finish().adapter_hits["illumina_universal"], 1)

    def test_renders_all_report_formats(self) -> None:
        result = calculate_qc(DATA / "reads_R1.fastq").as_dict()
        payload = {
            "sample": "sample-1",
            "generated_at": "2026-08-20T10:00:00+00:00",
            "tool_version": "0.1.0",
            "results": {"read1": result},
        }
        self.assertIn('"sample": "sample-1"', render_report(payload, "json"))
        self.assertIn("mate\tsource\tmetric\tvalue", render_report(payload, "tsv"))
        rendered_html = render_report(payload, "html")
        self.assertIn("<!doctype html>", rendered_html)
        self.assertIn("Mean quality by cycle", rendered_html)
        self.assertIn("<svg", rendered_html)

    def test_infers_report_format(self) -> None:
        self.assertEqual(infer_report_format("report.html", "auto"), "html")
        self.assertEqual(infer_report_format("report.tsv", "auto"), "tsv")
        self.assertEqual(infer_report_format("-", "auto"), "json")
        self.assertEqual(infer_report_format("report.tsv", "json"), "json")

    def test_rejects_invalid_max_reads(self) -> None:
        with self.assertRaises(ValueError):
            calculate_qc(DATA / "reads_R1.fastq", max_reads=0)


if __name__ == "__main__":
    unittest.main()
