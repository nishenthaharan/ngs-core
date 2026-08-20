"""Core primitives for streaming NGS file validation."""

__version__ = "0.1.0"

from .fastq import FastqRecord, read_fastq, read_paired_fastq

__all__ = ["FastqRecord", "read_fastq", "read_paired_fastq"]
