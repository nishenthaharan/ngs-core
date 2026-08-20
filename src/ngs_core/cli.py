"""Command-line interface for NGS Core."""

from __future__ import annotations

import argparse
import json
import sys
from typing import Sequence

from . import __version__
from .exceptions import NGSCoreError
from .fastq import read_fastq, read_paired_fastq


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ngs-core",
        description="Streaming FASTQ validation for reproducible NGS workflows.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser("validate", help="validate FASTQ structure and pairing")
    validate.add_argument("input", help="read 1 FASTQ/FASTQ.GZ path, or - for stdin")
    validate.add_argument("--read2", help="optional paired read 2 FASTQ/FASTQ.GZ")
    validate.add_argument(
        "--skip-pair-validation",
        action="store_true",
        help="validate files independently without matching read identifiers",
    )
    validate.set_defaults(handler=_run_validate)
    return parser


def _run_validate(args: argparse.Namespace) -> int:
    reads = 0
    bases = 0
    if args.read2:
        if args.input == "-" or args.read2 == "-":
            raise NGSCoreError("paired validation does not support stdin")
        for record1, record2 in read_paired_fastq(
            args.input,
            args.read2,
            validate_names=not args.skip_pair_validation,
        ):
            reads += 2
            bases += len(record1.sequence) + len(record2.sequence)
    else:
        for record in read_fastq(args.input):
            reads += 1
            bases += len(record.sequence)
    json.dump(
        {"valid": True, "files": 2 if args.read2 else 1, "reads": reads, "bases": bases},
        sys.stdout,
        indent=2,
    )
    sys.stdout.write("\n")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI and translate expected failures into concise diagnostics."""

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except (NGSCoreError, OSError, UnicodeError) as error:
        print(f"ngs-core: error: {error}", file=sys.stderr)
        return 2
