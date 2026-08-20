"""Core primitives for streaming NGS validation and quality control."""

__version__ = "0.2.0"

from .fastq import FastqRecord, read_fastq, read_paired_fastq
from .qc import QCResult, calculate_qc

__all__ = [
    "FastqRecord",
    "QCResult",
    "calculate_qc",
    "read_fastq",
    "read_paired_fastq",
]
