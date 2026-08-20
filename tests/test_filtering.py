from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ngs_core.exceptions import ConfigurationError, PairingError
from ngs_core.fastq import FastqRecord, read_fastq
from ngs_core.filtering import FilterConfig, filter_fastq, rejection_reason, trim_record

DATA = Path(__file__).parent / "data"


class FilteringTests(unittest.TestCase):
    def test_fixed_and_adapter_trimming(self) -> None:
        record = FastqRecord(
            "@r1",
            "GGACGTAGATCGGAAGAGTT",
            "+",
            "I" * 20,
        )
        config = FilterConfig(
            min_length=1,
            trim_front=2,
            trim_tail=2,
            quality_trim=None,
            adapter="AGATCGGAAGAG",
        )
        trimmed = trim_record(record, config)
        self.assertEqual(trimmed.sequence, "ACGT")
        self.assertEqual(len(trimmed.sequence), len(trimmed.quality))

    def test_quality_trims_low_quality_tail(self) -> None:
        record = FastqRecord("@r1", "ACGTACGT", "+", "IIII!!!!")
        config = FilterConfig(
            min_length=1,
            min_mean_quality=0,
            quality_trim=20,
            quality_window=4,
        )
        self.assertEqual(trim_record(record, config).sequence, "ACGT")

    def test_rejection_reasons(self) -> None:
        short = FastqRecord("@short", "ACGT", "+", "IIII")
        self.assertEqual(rejection_reason(short, FilterConfig()), "too_short")
        n_read = FastqRecord("@n", "ACGTNN", "+", "IIIIII")
        config = FilterConfig(min_length=1, max_n_fraction=0.1, quality_trim=None)
        self.assertEqual(rejection_reason(n_read, config), "too_many_n")
        low = FastqRecord("@low", "ACGT", "+", "!!!!")
        config = FilterConfig(min_length=1, quality_trim=None)
        self.assertEqual(rejection_reason(low, config), "low_mean_quality")

    def test_filters_single_end(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "filtered.fastq.gz"
            stats = filter_fastq(
                DATA / "reads_R1.fastq",
                output,
                FilterConfig(min_length=10, max_n_fraction=1.0, quality_trim=None),
            )
            self.assertEqual(stats.input_reads, 3)
            self.assertEqual(stats.output_reads, 2)
            self.assertEqual(stats.discarded["low_mean_quality"], 1)
            self.assertEqual(len(list(read_fastq(output))), 2)

    def test_filters_pairs_together(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output1 = Path(directory) / "R1.fastq"
            output2 = Path(directory) / "R2.fastq"
            stats = filter_fastq(
                DATA / "reads_R1.fastq",
                output1,
                FilterConfig(min_length=10, max_n_fraction=1.0, quality_trim=None),
                read2=DATA / "reads_R2.fastq",
                output2=output2,
            )
            self.assertEqual(stats.input_reads, 6)
            self.assertEqual(stats.output_reads, 4)
            self.assertEqual(len(list(read_fastq(output1))), 2)
            self.assertEqual(len(list(read_fastq(output2))), 2)

    def test_requires_both_paired_outputs(self) -> None:
        with (
            tempfile.TemporaryDirectory() as directory,
            self.assertRaises(ConfigurationError),
        ):
            filter_fastq(
                DATA / "reads_R1.fastq",
                Path(directory) / "out.fastq",
                FilterConfig(),
                read2=DATA / "reads_R2.fastq",
            )

    def test_rejects_output_over_input(self) -> None:
        config = FilterConfig(min_length=1)
        with self.assertRaisesRegex(ConfigurationError, "must be different"):
            filter_fastq(DATA / "reads_R1.fastq", DATA / "reads_R1.fastq", config)

    def test_failed_pairing_does_not_leave_partial_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            read2 = directory_path / "bad_R2.fastq"
            read2.write_text("@wrong/2\nACGT\n+\nIIII\n", encoding="utf-8")
            output1 = directory_path / "output_R1.fastq"
            output2 = directory_path / "output_R2.fastq"
            with self.assertRaises(PairingError):
                filter_fastq(
                    DATA / "reads_R1.fastq",
                    output1,
                    FilterConfig(min_length=1),
                    read2=read2,
                    output2=output2,
                )
            self.assertFalse(output1.exists())
            self.assertFalse(output2.exists())


if __name__ == "__main__":
    unittest.main()
