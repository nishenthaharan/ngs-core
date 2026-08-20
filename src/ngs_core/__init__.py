"""Core primitives for streaming NGS quality control and preprocessing."""

__version__ = "0.3.0"

from .fastq import FastqRecord, read_fastq, read_paired_fastq
from .filtering import FilterConfig, FilterStats, filter_fastq
from .qc import QCResult, calculate_qc

__all__ = [
    "FastqRecord",
    "FilterConfig",
    "FilterStats",
    "QCResult",
    "calculate_qc",
    "filter_fastq",
    "read_fastq",
    "read_paired_fastq",
]
