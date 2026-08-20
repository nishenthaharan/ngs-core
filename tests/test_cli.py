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

    def test_validate_paired_reads(self) -> None:
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            exit_code = main(
                [
                    "validate",
                    str(DATA / "reads_R1.fastq"),
                    "--read2",
                    str(DATA / "reads_R2.fastq"),
                ]
            )
        payload = json.loads(stdout.getvalue())
        self.assertEqual(exit_code, 0)
        self.assertEqual(payload["files"], 2)
        self.assertEqual(payload["reads"], 6)

    def test_bad_fastq_has_concise_exit_code(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.fastq"
            path.write_text("not-fastq\n", encoding="utf-8")
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                exit_code = main(["validate", str(path)])
            self.assertEqual(exit_code, 2)
            self.assertIn("ngs-core: error:", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
