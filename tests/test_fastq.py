from __future__ import annotations

import gzip
import io
import tempfile
import unittest
from pathlib import Path

from ngs_core.exceptions import FastqFormatError, PairingError
from ngs_core.fastq import FastqRecord, open_fastq_text, read_fastq, read_paired_fastq


DATA = Path(__file__).parent / "data"


class FastqReaderTests(unittest.TestCase):
    def test_reads_valid_records(self) -> None:
        records = list(read_fastq(DATA / "reads_R1.fastq"))
        self.assertEqual(len(records), 3)
        self.assertEqual(records[0].identifier, "read1/1")
        self.assertEqual(records[0].pair_key, "read1")
        self.assertEqual(records[0].sequence, "ACGTACGTNN")

    def test_serializes_record(self) -> None:
        record = FastqRecord("@r1", "ACGT", "+", "IIII")
        self.assertEqual(record.to_fastq(), "@r1\nACGT\n+\nIIII\n")

    def test_detects_gzip_without_gz_suffix(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.bin"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("@r1\nACGT\n+\nIIII\n")
            self.assertEqual(len(list(read_fastq(path))), 1)

    def test_writes_gzip_by_suffix(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "out.fastq.gz"
            with open_fastq_text(path, "wt") as handle:
                handle.write("@r1\nACGT\n+\nIIII\n")
            self.assertEqual(path.read_bytes()[:2], b"\x1f\x8b")
            self.assertEqual(len(list(read_fastq(path))), 1)

    def test_rejects_truncated_record(self) -> None:
        with self.assertRaisesRegex(FastqFormatError, "truncated"):
            list(read_fastq(io.StringIO("@r1\nACGT\n+\n")))

    def test_rejects_length_mismatch(self) -> None:
        with self.assertRaisesRegex(FastqFormatError, "sequence length 4"):
            list(read_fastq(io.StringIO("@r1\nACGT\n+\nIII\n")))

    def test_reads_paired_records(self) -> None:
        pairs = list(
            read_paired_fastq(DATA / "reads_R1.fastq", DATA / "reads_R2.fastq")
        )
        self.assertEqual(len(pairs), 3)
        self.assertEqual(pairs[1][0].pair_key, pairs[1][1].pair_key)

    def test_rejects_mismatched_pair_names(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            read2 = Path(directory) / "R2.fastq"
            read2.write_text("@different/2\nACGT\n+\nIIII\n", encoding="utf-8")
            read1 = Path(directory) / "R1.fastq"
            read1.write_text("@expected/1\nACGT\n+\nIIII\n", encoding="utf-8")
            with self.assertRaisesRegex(PairingError, "out of sync"):
                list(read_paired_fastq(read1, read2))

    def test_rejects_unequal_pair_counts(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            read1 = Path(directory) / "R1.fastq"
            read2 = Path(directory) / "R2.fastq"
            read1.write_text("@r1/1\nACGT\n+\nIIII\n", encoding="utf-8")
            read2.write_text("", encoding="utf-8")
            with self.assertRaisesRegex(PairingError, "different read counts"):
                list(read_paired_fastq(read1, read2))


if __name__ == "__main__":
    unittest.main()
