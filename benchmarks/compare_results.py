"""Compare median benchmark observations from baseline and candidate revisions."""

from __future__ import annotations

import argparse
import json
import statistics
from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path
from typing import Any

BenchmarkKey = tuple[str, str]


def summarize(payload: dict[str, Any]) -> dict[BenchmarkKey, dict[str, float | int]]:
    """Summarize successful observations by profile and operation."""

    grouped: defaultdict[BenchmarkKey, list[dict[str, Any]]] = defaultdict(list)
    for observation in payload.get("observations", []):
        if observation.get("return_code") == 0:
            key = (str(observation["profile"]), str(observation["operation"]))
            grouped[key].append(observation)

    summaries = {}
    for key, rows in grouped.items():
        summaries[key] = {
            "observations": len(rows),
            "median_wall_seconds": statistics.median(row["wall_seconds"] for row in rows),
            "median_reads_per_second": statistics.median(
                row["reads_per_second"] for row in rows
            ),
            "median_bases_per_second": statistics.median(
                row["bases_per_second"] for row in rows
            ),
        }
    return summaries


def _percent_change(baseline: float, candidate: float) -> float:
    return round(100.0 * (candidate - baseline) / baseline, 4) if baseline else 0.0


def compare_payloads(
    baseline: dict[str, Any], candidate: dict[str, Any]
) -> list[dict[str, float | int | str]]:
    """Return comparable median changes for matching benchmark groups."""

    baseline_summary = summarize(baseline)
    candidate_summary = summarize(candidate)
    comparisons = []
    for profile, operation in sorted(baseline_summary.keys() & candidate_summary.keys()):
        before = baseline_summary[(profile, operation)]
        after = candidate_summary[(profile, operation)]
        comparisons.append(
            {
                "profile": profile,
                "operation": operation,
                "baseline_observations": before["observations"],
                "candidate_observations": after["observations"],
                "baseline_wall_seconds": before["median_wall_seconds"],
                "candidate_wall_seconds": after["median_wall_seconds"],
                "wall_regression_percent": _percent_change(
                    float(before["median_wall_seconds"]),
                    float(after["median_wall_seconds"]),
                ),
                "throughput_change_percent": _percent_change(
                    float(before["median_reads_per_second"]),
                    float(after["median_reads_per_second"]),
                ),
            }
        )
    return comparisons


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-wall-regression-percent", type=float, default=10.0)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    baseline = json.loads(args.baseline.read_text(encoding="utf-8"))
    candidate = json.loads(args.candidate.read_text(encoding="utf-8"))
    comparisons = compare_payloads(baseline, candidate)
    payload = {
        "schema_version": "1.0",
        "max_wall_regression_percent": args.max_wall_regression_percent,
        "comparisons": comparisons,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if not comparisons:
        return 2
    return int(
        any(
            row["wall_regression_percent"] > args.max_wall_regression_percent
            for row in comparisons
        )
    )


if __name__ == "__main__":
    raise SystemExit(main())
