"""Streaming FASTQ input/output and paired-read validation.

The parser deliberately handles one record at a time so multi-gigabyte inputs do
not need to fit in memory. Both plain-text and gzip-compressed files are detected
from their bytes, rather than trusting the filename alone.
"""

from __future__ import annotations

import gzip
import io
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from .exceptions import FastqFormatError, PairingError

GZIP_MAGIC = b"\x1f\x8b"


@dataclass(frozen=True, slots=True)
class FastqRecord:
    """A validated FASTQ record."""

    name: str
    sequence: str
    separator: str
    quality: str

    @property
    def identifier(self) -> str:
        """Return the whitespace-delimited read identifier without ``@``."""

        return self.name[1:].split(maxsplit=1)[0]

    @property
    def pair_key(self) -> str:
        """Return a pair-neutral identifier for common FASTQ naming schemes."""

        identifier = self.identifier
        if identifier.endswith(("/1", "/2")):
            return identifier[:-2]
        return identifier

    def to_fastq(self) -> str:
        """Serialize the record using canonical four-line FASTQ layout."""

        return f"{self.name}\n{self.sequence}\n{self.separator}\n{self.quality}\n"


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == GZIP_MAGIC


@contextmanager
def open_fastq_text(
    source: str | Path,
    mode: str = "rt",
    *,
    compresslevel: int = 6,
) -> Iterator[TextIO]:
    """Open a FASTQ path or ``-`` for stdin/stdout.

    Reading detects gzip by magic bytes. Writing uses gzip when the destination
    ends in ``.gz``. Text is decoded as UTF-8 with universal newline handling.
    """

    if mode not in {"rt", "wt"}:
        raise ValueError("mode must be 'rt' or 'wt'")

    if str(source) == "-":
        handle = sys.stdin if mode == "rt" else sys.stdout
        yield handle
        return

    path = Path(source)
    if mode == "rt":
        if not path.is_file():
            raise FastqFormatError(f"FASTQ input does not exist: {path}")
        if _is_gzip(path):
            with gzip.open(path, mode="rt", encoding="utf-8", newline=None) as handle:
                yield handle
        else:
            with path.open(mode="rt", encoding="utf-8", newline=None) as handle:
                yield handle
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".gz":
        with gzip.open(
            path,
            mode="wt",
            encoding="utf-8",
            newline="\n",
            compresslevel=compresslevel,
        ) as handle:
            yield handle
    else:
        with path.open(mode="wt", encoding="utf-8", newline="\n") as handle:
            yield handle


def _clean_line(line: str) -> str:
    return line.rstrip("\r\n")


def read_fastq(
    source: str | Path | TextIO,
    *,
    validate_bases: bool = True,
) -> Iterator[FastqRecord]:
    """Yield validated four-line FASTQ records from ``source``.

    ``source`` may be a path, ``-``, or an already-open text handle. Empty files
    are valid and yield no records. Wrapped FASTQ is rejected explicitly.
    """

    if isinstance(source, io.TextIOBase) or hasattr(source, "readline"):
        yield from _read_fastq_handle(source, "<stream>", validate_bases)
        return

    with open_fastq_text(source, "rt") as handle:
        yield from _read_fastq_handle(handle, str(source), validate_bases)


def _read_fastq_handle(
    handle: TextIO,
    source_name: str,
    validate_bases: bool,
) -> Iterator[FastqRecord]:
    record_number = 0
    while True:
        raw_name = handle.readline()
        if raw_name == "":
            return

        record_number += 1
        raw_sequence = handle.readline()
        raw_separator = handle.readline()
        raw_quality = handle.readline()
        if "" in (raw_sequence, raw_separator, raw_quality):
            raise FastqFormatError(f"{source_name}: truncated FASTQ record {record_number}")

        name = _clean_line(raw_name)
        sequence = _clean_line(raw_sequence).upper()
        separator = _clean_line(raw_separator)
        quality = _clean_line(raw_quality)

        if not name.startswith("@") or len(name) == 1:
            raise FastqFormatError(
                f"{source_name}: record {record_number} header must start with '@'"
            )
        if not separator.startswith("+"):
            raise FastqFormatError(
                f"{source_name}: record {record_number} separator must start with '+'"
            )
        if not sequence:
            raise FastqFormatError(f"{source_name}: record {record_number} has an empty sequence")
        if len(sequence) != len(quality):
            raise FastqFormatError(
                f"{source_name}: record {record_number} has sequence length "
                f"{len(sequence)} but quality length {len(quality)}"
            )
        if any(char.isspace() for char in sequence + quality):
            raise FastqFormatError(
                f"{source_name}: record {record_number} contains whitespace in sequence or quality"
            )
        if validate_bases and any(not (char.isalpha() or char in ".-") for char in sequence):
            raise FastqFormatError(
                f"{source_name}: record {record_number} contains invalid base characters"
            )

        yield FastqRecord(name, sequence, separator, quality)


def read_paired_fastq(
    read1: str | Path,
    read2: str | Path,
    *,
    validate_names: bool = True,
) -> Iterator[tuple[FastqRecord, FastqRecord]]:
    """Yield synchronized paired-end records without loading either file."""

    iterator1 = iter(read_fastq(read1))
    iterator2 = iter(read_fastq(read2))
    pair_number = 0

    while True:
        record1 = next(iterator1, None)
        record2 = next(iterator2, None)
        if record1 is None and record2 is None:
            return

        pair_number += 1
        if record1 is None or record2 is None:
            shorter = str(read1) if record1 is None else str(read2)
            raise PairingError(
                f"Paired inputs have different read counts; {shorter} ended "
                f"before pair {pair_number}"
            )
        if validate_names and record1.pair_key != record2.pair_key:
            raise PairingError(
                f"Pair {pair_number} is out of sync: {record1.identifier!r} != "
                f"{record2.identifier!r}"
            )
        yield record1, record2


def quality_scores(quality: str, phred_offset: int = 33) -> list[int]:
    """Decode an ASCII quality string and validate the resulting score range."""

    scores = [ord(char) - phred_offset for char in quality]
    if any(score < 0 or score > 93 for score in scores):
        raise FastqFormatError(f"Quality string is incompatible with Phred+{phred_offset} encoding")
    return scores
