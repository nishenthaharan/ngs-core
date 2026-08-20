"""Command-line interface for NGS Core."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__
from .exceptions import ConfigurationError, NGSCoreError
from .fastq import read_fastq, read_paired_fastq
from .qc import QCAccumulator, calculate_qc
from .report import infer_report_format, render_report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ngs-core",
        description="Streaming FASTQ validation and quality control.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    qc = subparsers.add_parser("qc", help="calculate FASTQ quality-control metrics")
    qc.add_argument("input", help="read 1 FASTQ/FASTQ.GZ path, or - for stdin")
    qc.add_argument("--read2", help="optional paired read 2 FASTQ/FASTQ.GZ")
    qc.add_argument("-o", "--output", default="-", help="report path; default: stdout")
    qc.add_argument(
        "--format",
        choices=("auto", "json", "tsv", "html"),
        default="auto",
        help="report format inferred from output suffix by default",
    )
    qc.add_argument("--sample", help="sample label; defaults to the input filename")
    qc.add_argument("--phred-offset", type=int, choices=(33, 64), default=33)
    qc.add_argument("--max-reads", type=_positive_int, help="analyze only the first N reads")
    qc.add_argument(
        "--skip-pair-validation",
        action="store_true",
        help="do not compare paired read identifiers",
    )
    qc.set_defaults(handler=_run_qc)

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


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def _same_path(first: str, second: str) -> bool:
    if first == "-" or second == "-":
        return first == second
    return Path(first).resolve(strict=False) == Path(second).resolve(strict=False)


def _reject_report_collision(output: str, inputs: Sequence[str | None]) -> None:
    if output == "-":
        return
    if any(_same_path(output, input_path) for input_path in inputs if input_path is not None):
        raise ConfigurationError("report path must not overwrite a FASTQ input")


def _paired_qc(args: argparse.Namespace) -> dict[str, Any]:
    accumulator1 = QCAccumulator(args.input, phred_offset=args.phred_offset)
    accumulator2 = QCAccumulator(args.read2, phred_offset=args.phred_offset)
    for index, (record1, record2) in enumerate(
        read_paired_fastq(
            args.input,
            args.read2,
            validate_names=not args.skip_pair_validation,
        ),
        start=1,
    ):
        accumulator1.add(record1)
        accumulator2.add(record2)
        if args.max_reads is not None and index >= args.max_reads:
            break
    return {
        "read1": accumulator1.finish().as_dict(),
        "read2": accumulator2.finish().as_dict(),
    }


def _run_qc(args: argparse.Namespace) -> int:
    if args.read2 and args.input == "-":
        raise ConfigurationError("paired QC does not support stdin")
    _reject_report_collision(args.output, [args.input, args.read2])
    results = (
        _paired_qc(args)
        if args.read2
        else {
            "read1": calculate_qc(
                args.input,
                phred_offset=args.phred_offset,
                max_reads=args.max_reads,
            ).as_dict()
        }
    )
    sample = args.sample or ("stdin" if args.input == "-" else Path(args.input).name)
    payload = {
        "schema_version": "1.0",
        "tool": "ngs-core",
        "tool_version": __version__,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "sample": sample,
        "results": results,
    }
    report_format = infer_report_format(args.output, args.format)
    rendered = render_report(payload, report_format)
    if args.output == "-":
        sys.stdout.write(rendered)
    else:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
        print(f"QC report written to {output}", file=sys.stderr)
    return 0


def _run_validate(args: argparse.Namespace) -> int:
    reads = 0
    bases = 0
    if args.read2:
        if args.input == "-" or args.read2 == "-":
            raise ConfigurationError("paired validation does not support stdin")
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
