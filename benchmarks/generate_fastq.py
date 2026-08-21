"""Generate deterministic synthetic FASTQ inputs for repeatable benchmarks."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import random
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO

DNA_ALPHABET = "ACGT"


@contextmanager
def _open_output(path: Path) -> Iterator[TextIO]:
    """Open plain text or byte-reproducible gzip output."""

    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() != ".gz":
        with path.open("wt", encoding="ascii", newline="\n") as handle:
            yield handle
        return

    # A fixed timestamp and empty embedded filename make gzip bytes reproducible across
    # output paths. This matters because manifests compare checksums between benchmark runs.
    with (
        path.open("wb") as raw_handle,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as gzip_handle,
        io.TextIOWrapper(gzip_handle, encoding="ascii", newline="\n") as text_handle,
    ):
        yield text_handle


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def generate_fastq(
    output: str | Path,
    *,
    records: int,
    read_length: int,
    seed: int,
    mate: int = 1,
    min_read_length: int | None = None,
) -> dict[str, int | str]:
    """Generate one deterministic FASTQ file and return its manifest entry."""

    minimum = read_length if min_read_length is None else min_read_length
    if records < 1:
        raise ValueError("records must be at least 1")
    if minimum < 1 or read_length < minimum:
        raise ValueError("read lengths must satisfy 1 <= minimum <= maximum")
    if mate not in {1, 2}:
        raise ValueError("mate must be 1 or 2")

    path = Path(output)
    rng = random.Random(seed + mate * 1_000_003)
    bases = 0
    with _open_output(path) as handle:
        for index in range(1, records + 1):
            length = rng.randint(minimum, read_length)
            sequence = "".join(rng.choice(DNA_ALPHABET) for _ in range(length))
            quality = "".join(chr(33 + rng.randint(20, 40)) for _ in range(length))
            handle.write(f"@synthetic:{index}/{mate}\n{sequence}\n+\n{quality}\n")
            bases += length

    return {
        "path": str(path),
        "mate": mate,
        "records": records,
        "bases": bases,
        "bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="R1 FASTQ or FASTQ.GZ destination")
    parser.add_argument("--output2", help="optional R2 FASTQ or FASTQ.GZ destination")
    parser.add_argument("--records", type=int, required=True)
    parser.add_argument("--read-length", type=int, required=True, help="maximum read length")
    parser.add_argument("--min-read-length", type=int, help="minimum variable read length")
    parser.add_argument("--seed", type=int, default=20260821)
    parser.add_argument("--manifest", help="write generation metadata as JSON")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.output2 and Path(args.output).resolve() == Path(args.output2).resolve():
        raise ValueError("paired outputs must be different paths")

    entries = [
        generate_fastq(
            args.output,
            records=args.records,
            read_length=args.read_length,
            min_read_length=args.min_read_length,
            seed=args.seed,
            mate=1,
        )
    ]
    if args.output2:
        entries.append(
            generate_fastq(
                args.output2,
                records=args.records,
                read_length=args.read_length,
                min_read_length=args.min_read_length,
                seed=args.seed,
                mate=2,
            )
        )

    manifest = {
        "schema_version": "1.0",
        "generator": "ngs-core",
        "seed": args.seed,
        "files": entries,
    }
    rendered = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    if args.manifest:
        manifest_path = Path(args.manifest)
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
