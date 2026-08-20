from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from ngs_core.cli import main

DATA = Path(__file__).parent / "data"


class CLITests(unittest.TestCase):
    def test_validate_outputs_machine_readable_summary(self) -> None:
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            exit_code = main(["validate", str(DATA / "reads_R1.fastq")])
        self.assertEqual(exit_code, 0)
        self.assertEqual(json.loads(stdout.getvalue())["reads"], 3)

    def test_qc_writes_html_report(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "qc.html"
            with contextlib.redirect_stderr(io.StringIO()):
                exit_code = main(["qc", str(DATA / "reads_R1.fastq"), "--output", str(output)])
            self.assertEqual(exit_code, 0)
            self.assertIn("<!doctype html>", output.read_text(encoding="utf-8"))

    def test_paired_qc_outputs_two_results(self) -> None:
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            exit_code = main(
                [
                    "qc",
                    str(DATA / "reads_R1.fastq"),
                    "--read2",
                    str(DATA / "reads_R2.fastq"),
                    "--max-reads",
                    "2",
                ]
            )
        payload = json.loads(stdout.getvalue())
        self.assertEqual(exit_code, 0)
        self.assertEqual(set(payload["results"]), {"read1", "read2"})
        self.assertEqual(payload["results"]["read1"]["reads"], 2)

    def test_bad_fastq_has_concise_exit_code(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.fastq"
            path.write_text("not-fastq\n", encoding="utf-8")
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                exit_code = main(["validate", str(path)])
            self.assertEqual(exit_code, 2)
            self.assertIn("ngs-core: error:", stderr.getvalue())

    def test_qc_refuses_to_overwrite_input(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "reads.fastq"
            original = "@r1\nACGT\n+\nIIII\n"
            input_path.write_text(original, encoding="utf-8")
            with contextlib.redirect_stderr(io.StringIO()):
                exit_code = main(["qc", str(input_path), "--output", str(input_path)])
            self.assertEqual(exit_code, 2)
            self.assertEqual(input_path.read_text(encoding="utf-8"), original)


if __name__ == "__main__":
    unittest.main()
