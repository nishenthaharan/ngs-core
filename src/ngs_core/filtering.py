"""Read trimming and filtering for single- and paired-end FASTQ data."""

from __future__ import annotations

import os
import tempfile
from collections import Counter
from contextlib import ExitStack, contextmanager
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from .exceptions import ConfigurationError
from .fastq import FastqRecord, open_fastq_text, quality_scores, read_fastq, read_paired_fastq


@dataclass(frozen=True, slots=True)
class FilterConfig:
    """Parameters controlling deterministic trimming and filtering."""

    min_length: int = 50
    min_mean_quality: float = 20.0
    max_n_fraction: float = 0.05
    trim_front: int = 0
    trim_tail: int = 0
    quality_trim: int | None = 20
    quality_window: int = 4
    adapter: str | None = None
    phred_offset: int = 33

    def __post_init__(self) -> None:
        if self.min_length < 1:
            raise ConfigurationError("min_length must be at least 1")
        if not 0 <= self.min_mean_quality <= 93:
            raise ConfigurationError("min_mean_quality must be between 0 and 93")
        if not 0 <= self.max_n_fraction <= 1:
            raise ConfigurationError("max_n_fraction must be between 0 and 1")
        if self.trim_front < 0 or self.trim_tail < 0:
            raise ConfigurationError("fixed trimming values cannot be negative")
        if self.quality_trim is not None and not 0 <= self.quality_trim <= 93:
            raise ConfigurationError("quality_trim must be between 0 and 93")
        if self.quality_window < 1:
            raise ConfigurationError("quality_window must be at least 1")
        if self.phred_offset not in {33, 64}:
            raise ConfigurationError("phred_offset must be 33 or 64")
        if self.adapter is not None:
            adapter = self.adapter.upper()
            if len(adapter) < 6 or any(not char.isalpha() for char in adapter):
                raise ConfigurationError(
                    "adapter must contain at least six alphabetic bases"
                )


@dataclass(slots=True)
class FilterStats:
    """Summary of a filtering operation."""

    input_reads: int = 0
    output_reads: int = 0
    input_bases: int = 0
    output_bases: int = 0
    trimmed_bases: int = 0
    discarded: Counter[str] = field(default_factory=Counter)

    def as_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["discarded"] = dict(sorted(self.discarded.items()))
        result["retained_percent"] = (
            round(100 * self.output_reads / self.input_reads, 4)
            if self.input_reads
            else 0.0
        )
        return result


def trim_record(record: FastqRecord, config: FilterConfig) -> FastqRecord:
    """Apply fixed, adapter, and 3' sliding-window quality trimming."""

    start = min(config.trim_front, len(record.sequence))
    stop = len(record.sequence) - min(config.trim_tail, len(record.sequence) - start)
    sequence = record.sequence[start:stop]
    quality = record.quality[start:stop]

    if config.adapter:
        adapter_index = sequence.find(config.adapter.upper())
        if adapter_index >= 0:
            sequence = sequence[:adapter_index]
            quality = quality[:adapter_index]

    if config.quality_trim is not None and quality:
        scores = quality_scores(quality, config.phred_offset)
        end = len(scores)
        while end >= config.quality_window:
            window = scores[end - config.quality_window : end]
            if sum(window) / len(window) >= config.quality_trim:
                break
            end -= 1
        while end > 0 and scores[end - 1] < config.quality_trim:
            end -= 1
        sequence = sequence[:end]
        quality = quality[:end]

    return FastqRecord(record.name, sequence, record.separator, quality)


def rejection_reason(record: FastqRecord, config: FilterConfig) -> str | None:
    """Return the first deterministic reason a trimmed record should be removed."""

    length = len(record.sequence)
    if length < config.min_length:
        return "too_short"
    n_fraction = record.sequence.upper().count("N") / length
    if n_fraction > config.max_n_fraction:
        return "too_many_n"
    scores = quality_scores(record.quality, config.phred_offset)
    if sum(scores) / length < config.min_mean_quality:
        return "low_mean_quality"
    return None


def _record_input(stats: FilterStats, *records: FastqRecord) -> None:
    stats.input_reads += len(records)
    stats.input_bases += sum(len(record.sequence) for record in records)


def _record_output(
    stats: FilterStats,
    trimmed_records: tuple[FastqRecord, ...],
) -> None:
    stats.output_reads += len(trimmed_records)
    stats.output_bases += sum(len(record.sequence) for record in trimmed_records)


def _record_trimming(
    stats: FilterStats,
    original_records: tuple[FastqRecord, ...],
    trimmed_records: tuple[FastqRecord, ...],
) -> None:
    stats.trimmed_bases += sum(
        len(original.sequence) - len(trimmed.sequence)
        for original, trimmed in zip(original_records, trimmed_records, strict=True)
    )


def _same_file_path(first: str | Path, second: str | Path) -> bool:
    if str(first) == "-" or str(second) == "-":
        return str(first) == str(second)
    return Path(first).resolve(strict=False) == Path(second).resolve(strict=False)


@contextmanager
def _atomic_fastq_output(destination: str | Path):
    """Write a path atomically, leaving no incomplete output after an error."""

    if str(destination) == "-":
        with open_fastq_text(destination, "wt") as handle:
            yield handle
        return

    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = ".gz" if path.suffix.lower() == ".gz" else ".tmp"
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=suffix,
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    try:
        with open_fastq_text(temporary_path, "wt") as handle:
            yield handle
        # The temporary file lives beside the destination, so os.replace performs an
        # atomic same-filesystem handoff after the complete stream has been written.
        os.replace(temporary_path, path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def filter_fastq(
    read1: str | Path,
    output1: str | Path,
    config: FilterConfig,
    *,
    read2: str | Path | None = None,
    output2: str | Path | None = None,
    validate_names: bool = True,
) -> FilterStats:
    """Filter FASTQ input while keeping paired outputs synchronized.

    For paired data, a pair is retained only when both mates pass. Statistics
    count reads (not pairs), making single- and paired-end runs comparable.
    """

    if (read2 is None) != (output2 is None):
        raise ConfigurationError("read2 and output2 must be supplied together")
    inputs = [read1] + ([read2] if read2 is not None else [])
    outputs = [output1] + ([output2] if output2 is not None else [])
    # Check the complete input/output cross-product before opening a writer. This prevents
    # accidental truncation even when equivalent paths use different relative spellings.
    if any(
        _same_file_path(input_path, output_path)
        for input_path in inputs
        for output_path in outputs
    ):
        raise ConfigurationError("input and output paths must be different")
    if output2 is not None and _same_file_path(output1, output2):
        raise ConfigurationError("paired output paths must be different")
    if read2 is not None and _same_file_path(read1, read2):
        raise ConfigurationError("paired input paths must be different")
    if read2 is not None and (str(read1) == "-" or str(read2) == "-"):
        raise ConfigurationError("paired-end filtering does not support stdin")
    if output2 is not None and (str(output1) == "-" or str(output2) == "-"):
        raise ConfigurationError("paired-end filtering does not support stdout")

    stats = FilterStats()
    with ExitStack() as stack:
        handle1 = stack.enter_context(_atomic_fastq_output(output1))
        handle2 = (
            stack.enter_context(_atomic_fastq_output(output2))
            if output2 is not None
            else None
        )

        if read2 is None:
            for original in read_fastq(read1):
                _record_input(stats, original)
                trimmed = trim_record(original, config)
                _record_trimming(stats, (original,), (trimmed,))
                reason = rejection_reason(trimmed, config)
                if reason:
                    stats.discarded[reason] += 1
                    continue
                handle1.write(trimmed.to_fastq())
                _record_output(stats, (trimmed,))
            return stats

        assert handle2 is not None
        for original1, original2 in read_paired_fastq(
            read1, read2, validate_names=validate_names
        ):
            _record_input(stats, original1, original2)
            trimmed1 = trim_record(original1, config)
            trimmed2 = trim_record(original2, config)
            _record_trimming(
                stats,
                (original1, original2),
                (trimmed1, trimmed2),
            )
            reason1 = rejection_reason(trimmed1, config)
            reason2 = rejection_reason(trimmed2, config)
            # Never emit an orphaned mate: rejection of either read rejects the full pair.
            if reason1 or reason2:
                stats.discarded[f"pair_{reason1 or reason2}"] += 2
                continue
            handle1.write(trimmed1.to_fastq())
            handle2.write(trimmed2.to_fastq())
            _record_output(stats, (trimmed1, trimmed2))
    return stats
