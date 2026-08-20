"""Memory-efficient quality-control metrics for FASTQ reads."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

from .fastq import FastqRecord, quality_scores, read_fastq

DEFAULT_ADAPTERS = {
    "illumina_universal": "AGATCGGAAGAG",
    "nextera_transposase": "CTGTCTCTTATA",
    "small_rna_3prime": "TGGAATTCTCGG",
}


def _safe_percentage(numerator: int | float, denominator: int | float) -> float:
    return round((100.0 * numerator / denominator), 4) if denominator else 0.0


@dataclass(slots=True)
class QCResult:
    """Serializable QC summary for one FASTQ input."""

    source: str
    reads: int
    bases: int
    gc_percent: float
    n_percent: float
    q20_percent: float
    q30_percent: float
    mean_quality: float
    mean_length: float
    median_length: float
    min_length: int
    max_length: int
    read_length_n50: int
    duplicate_estimate_percent: float
    adapter_hits: dict[str, int]
    per_cycle: list[dict[str, Any]]
    length_distribution: list[dict[str, int]]
    phred_offset: int

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class QCAccumulator:
    """Accumulate exact streaming metrics and bounded duplication estimates."""

    source: str
    phred_offset: int = 33
    duplicate_sample_size: int = 100_000
    reads: int = 0
    bases: int = 0
    gc_bases: int = 0
    n_bases: int = 0
    q20_bases: int = 0
    q30_bases: int = 0
    quality_sum: int = 0
    min_length: int | None = None
    max_length: int = 0
    length_counts: Counter[int] = field(default_factory=Counter)
    cycle_counts: list[int] = field(default_factory=list)
    cycle_quality_sum: list[int] = field(default_factory=list)
    cycle_q20: list[int] = field(default_factory=list)
    cycle_q30: list[int] = field(default_factory=list)
    cycle_bases: list[Counter[str]] = field(default_factory=list)
    adapter_hits: Counter[str] = field(default_factory=Counter)
    duplicate_sequences: Counter[str] = field(default_factory=Counter)

    def add(self, record: FastqRecord) -> None:
        """Add one record to the running statistics."""

        scores = quality_scores(record.quality, self.phred_offset)
        sequence = record.sequence.upper()
        length = len(sequence)

        self.reads += 1
        self.bases += length
        self.gc_bases += sequence.count("G") + sequence.count("C")
        self.n_bases += sequence.count("N")
        self.quality_sum += sum(scores)
        self.q20_bases += sum(score >= 20 for score in scores)
        self.q30_bases += sum(score >= 30 for score in scores)
        self.min_length = length if self.min_length is None else min(self.min_length, length)
        self.max_length = max(self.max_length, length)
        self.length_counts[length] += 1

        missing = length - len(self.cycle_counts)
        if missing > 0:
            self.cycle_counts.extend([0] * missing)
            self.cycle_quality_sum.extend([0] * missing)
            self.cycle_q20.extend([0] * missing)
            self.cycle_q30.extend([0] * missing)
            self.cycle_bases.extend(Counter() for _ in range(missing))

        for index, (base, score) in enumerate(zip(sequence, scores, strict=True)):
            self.cycle_counts[index] += 1
            self.cycle_quality_sum[index] += score
            self.cycle_q20[index] += score >= 20
            self.cycle_q30[index] += score >= 30
            normalized = base if base in {"A", "C", "G", "T", "N"} else "OTHER"
            self.cycle_bases[index][normalized] += 1

        for name, adapter in DEFAULT_ADAPTERS.items():
            if adapter in sequence:
                self.adapter_hits[name] += 1

        if self.reads <= self.duplicate_sample_size:
            self.duplicate_sequences[sequence[:50]] += 1

    def finish(self) -> QCResult:
        """Create an immutable summary from accumulated values."""

        duplicate_total = sum(self.duplicate_sequences.values())
        unique = len(self.duplicate_sequences)
        duplicate_estimate = duplicate_total - unique

        per_cycle: list[dict[str, Any]] = []
        for index, count in enumerate(self.cycle_counts):
            base_counts = self.cycle_bases[index]
            per_cycle.append(
                {
                    "cycle": index + 1,
                    "observations": count,
                    "mean_quality": round(self.cycle_quality_sum[index] / count, 4),
                    "q20_percent": _safe_percentage(self.cycle_q20[index], count),
                    "q30_percent": _safe_percentage(self.cycle_q30[index], count),
                    "base_percent": {
                        base: _safe_percentage(base_counts.get(base, 0), count)
                        for base in ("A", "C", "G", "T", "N", "OTHER")
                    },
                }
            )

        return QCResult(
            source=self.source,
            reads=self.reads,
            bases=self.bases,
            gc_percent=_safe_percentage(self.gc_bases, self.bases),
            n_percent=_safe_percentage(self.n_bases, self.bases),
            q20_percent=_safe_percentage(self.q20_bases, self.bases),
            q30_percent=_safe_percentage(self.q30_bases, self.bases),
            mean_quality=round(self.quality_sum / self.bases, 4) if self.bases else 0.0,
            mean_length=round(self.bases / self.reads, 4) if self.reads else 0.0,
            median_length=_weighted_median(self.length_counts),
            min_length=self.min_length or 0,
            max_length=self.max_length,
            read_length_n50=_read_length_n50(self.length_counts),
            duplicate_estimate_percent=_safe_percentage(
                duplicate_estimate, duplicate_total
            ),
            adapter_hits={
                name: self.adapter_hits.get(name, 0) for name in DEFAULT_ADAPTERS
            },
            per_cycle=per_cycle,
            length_distribution=[
                {"length": length, "reads": count}
                for length, count in sorted(self.length_counts.items())
            ],
            phred_offset=self.phred_offset,
        )


def _weighted_median(counts: Counter[int]) -> float:
    total = sum(counts.values())
    if total == 0:
        return 0.0
    lower_rank = (total - 1) // 2
    upper_rank = total // 2
    seen = 0
    lower_value = 0
    upper_value = 0
    for value, count in sorted(counts.items()):
        previous = seen
        seen += count
        if previous <= lower_rank < seen:
            lower_value = value
        if previous <= upper_rank < seen:
            upper_value = value
            break
    return (lower_value + upper_value) / 2


def _read_length_n50(counts: Counter[int]) -> int:
    total_bases = sum(length * count for length, count in counts.items())
    if total_bases == 0:
        return 0
    accumulated = 0
    for length, count in sorted(counts.items(), reverse=True):
        accumulated += length * count
        if accumulated * 2 >= total_bases:
            return length
    return 0


def calculate_qc(
    source: str | Path,
    *,
    phred_offset: int = 33,
    max_reads: int | None = None,
) -> QCResult:
    """Calculate a streaming QC result for one FASTQ file."""

    if max_reads is not None and max_reads < 1:
        raise ValueError("max_reads must be at least 1")
    accumulator = QCAccumulator(str(source), phred_offset=phred_offset)
    records: Iterable[FastqRecord] = read_fastq(source)
    for index, record in enumerate(records, start=1):
        accumulator.add(record)
        if max_reads is not None and index >= max_reads:
            break
    return accumulator.finish()
