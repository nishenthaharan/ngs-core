"""Measure NGS Core commands and write auditable JSON observations."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import tempfile
import time
from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _linux_rss_bytes(process_id: int) -> int | None:
    """Read one Linux process RSS sample without adding a runtime dependency."""

    try:
        status = Path(f"/proc/{process_id}/status").read_text(encoding="utf-8")
    except OSError:
        return None
    for line in status.splitlines():
        if line.startswith("VmRSS:"):
            try:
                return int(line.split()[1]) * 1024
            except (IndexError, ValueError):
                return None
    return None


def _environment() -> dict[str, str | int]:
    return {
        "python": platform.python_version(),
        "operating_system": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor() or "unknown",
        "cpu_count": os.cpu_count() or 0,
    }


def _revision() -> str:
    if revision := os.environ.get("GITHUB_SHA"):
        return revision
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def measure_command(
    command: Sequence[str],
    *,
    operation: str,
    profile: str,
    repetition: int,
    reads: int,
    bases: int,
    poll_interval: float = 0.01,
) -> dict[str, Any]:
    """Execute one command and return timing, throughput, and resource metadata."""

    child_times_before = os.times()
    started = time.perf_counter()
    peak_rss: int | None = None
    with (
        tempfile.TemporaryFile(mode="w+", encoding="utf-8") as stdout_handle,
        tempfile.TemporaryFile(mode="w+", encoding="utf-8") as stderr_handle,
    ):
        process = subprocess.Popen(
            list(command),
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
        )
        while process.poll() is None:
            if (rss := _linux_rss_bytes(process.pid)) is not None:
                peak_rss = rss if peak_rss is None else max(peak_rss, rss)
            time.sleep(poll_interval)
        return_code = process.wait()
        wall_seconds = time.perf_counter() - started
        child_times_after = os.times()
        cpu_seconds = (
            child_times_after.children_user
            + child_times_after.children_system
            - child_times_before.children_user
            - child_times_before.children_system
        )
        stdout_handle.seek(0)
        stderr_handle.seek(0)
        stdout = stdout_handle.read()
        stderr = stderr_handle.read()

    return {
        "schema_version": "1.0",
        "revision": _revision(),
        "operation": operation,
        "profile": profile,
        "repetition": repetition,
        "command": list(command),
        "reads": reads,
        "bases": bases,
        "wall_seconds": round(wall_seconds, 6),
        "cpu_seconds": round(cpu_seconds, 6),
        "peak_rss_bytes": peak_rss,
        "reads_per_second": round(reads / wall_seconds, 4) if wall_seconds else 0.0,
        "bases_per_second": round(bases / wall_seconds, 4) if wall_seconds else 0.0,
        "return_code": return_code,
        "stdout_sha256": hashlib.sha256(stdout.encode()).hexdigest(),
        "stderr": stderr[-2000:] if return_code else "",
        "environment": _environment(),
    }


def _build_command(
    executable: str,
    operation: str,
    input_path: Path,
    read2: Path | None,
    output_directory: Path,
    repetition: int,
) -> list[str]:
    command = [executable, operation, str(input_path)]
    if read2 is not None:
        command.extend(["--read2", str(read2)])
    if operation == "qc":
        command.extend(["--output", str(output_directory / f"qc-{repetition}.json")])
    elif operation == "filter":
        command.extend(
            [
                "--output",
                str(output_directory / f"filtered-R1-{repetition}.fastq"),
                "--min-length",
                "1",
                "--min-mean-quality",
                "0",
                "--max-n-fraction",
                "1",
                "--quality-trim",
                "-1",
            ]
        )
        if read2 is not None:
            command.extend(
                ["--output2", str(output_directory / f"filtered-R2-{repetition}.fastq")]
            )
    return command


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--operation", choices=("validate", "qc", "filter"), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--read2", type=Path)
    parser.add_argument("--profile", default="custom")
    parser.add_argument("--records", type=int, required=True)
    parser.add_argument("--bases", type=int, required=True)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--ngs-core", default="ngs-core", help="installed CLI executable")
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.records < 1 or args.bases < 1 or args.repetitions < 1:
        raise ValueError("records, bases, and repetitions must be positive")

    observations = []
    with tempfile.TemporaryDirectory() as directory:
        output_directory = Path(directory)
        for repetition in range(1, args.repetitions + 1):
            command = _build_command(
                args.ngs_core,
                args.operation,
                args.input,
                args.read2,
                output_directory,
                repetition,
            )
            observations.append(
                measure_command(
                    command,
                    operation=args.operation,
                    profile=args.profile,
                    repetition=repetition,
                    reads=args.records * (2 if args.read2 else 1),
                    bases=args.bases,
                )
            )

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "observations": observations,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0 if all(row["return_code"] == 0 for row in observations) else 1


if __name__ == "__main__":
    raise SystemExit(main())
