"""JSON, TSV, and self-contained HTML report rendering."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

from .exceptions import ConfigurationError


def infer_report_format(output: str | Path, requested: str) -> str:
    """Resolve ``auto`` from a destination suffix."""

    if requested != "auto":
        return requested
    if str(output) == "-":
        return "json"
    suffix = Path(output).suffix.lower()
    return {".html": "html", ".htm": "html", ".tsv": "tsv"}.get(suffix, "json")


def render_report(payload: dict[str, Any], report_format: str) -> str:
    """Render a QC payload without requiring a plotting or template dependency."""

    if report_format == "json":
        return json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if report_format == "tsv":
        return _render_tsv(payload)
    if report_format == "html":
        return _render_html(payload)
    raise ConfigurationError(f"unsupported report format: {report_format}")


def _render_tsv(payload: dict[str, Any]) -> str:
    metrics = (
        "reads",
        "bases",
        "gc_percent",
        "n_percent",
        "q20_percent",
        "q30_percent",
        "mean_quality",
        "mean_length",
        "median_length",
        "min_length",
        "max_length",
        "read_length_n50",
        "duplicate_estimate_percent",
    )
    lines = ["mate\tsource\tmetric\tvalue"]
    for mate, result in payload["results"].items():
        for metric in metrics:
            lines.append(f"{mate}\t{result['source']}\t{metric}\t{result[metric]}")
        for name, count in result["adapter_hits"].items():
            lines.append(f"{mate}\t{result['source']}\tadapter_{name}_hits\t{count}")
    return "\n".join(lines) + "\n"


def _line_chart(
    values: list[float],
    *,
    width: int = 820,
    height: int = 260,
    maximum: float = 42.0,
) -> str:
    """Return an inline SVG so the HTML report remains a single portable file."""

    if not values:
        return '<p class="empty">No reads were available for this chart.</p>'
    left, right, top, bottom = 48, 16, 16, 38
    plot_width = width - left - right
    plot_height = height - top - bottom
    denominator = max(len(values) - 1, 1)
    points = " ".join(
        f"{left + (index / denominator) * plot_width:.2f},"
        f"{top + (1 - min(value, maximum) / maximum) * plot_height:.2f}"
        for index, value in enumerate(values)
    )
    guides = []
    for score in (0, 10, 20, 30, 40):
        y = top + (1 - score / maximum) * plot_height
        guides.append(
            f'<line x1="{left}" y1="{y:.2f}" x2="{width - right}" y2="{y:.2f}" '
            'class="grid"/><text x="8" y="'
            f'{y + 4:.2f}" class="axis">Q{score}</text>'
        )
    return (
        f'<svg viewBox="0 0 {width} {height}" role="img" '
        'aria-label="Mean quality by sequencing cycle">'
        + "".join(guides)
        + f'<polyline points="{points}" class="quality-line"/>'
        + f'<text x="{width / 2:.0f}" y="{height - 6}" class="axis">Sequencing cycle</text>'
        + "</svg>"
    )


def _metric_card(label: str, value: Any, suffix: str = "") -> str:
    return (
        '<div class="metric"><span>'
        + html.escape(label)
        + "</span><strong>"
        + html.escape(f"{value}{suffix}")
        + "</strong></div>"
    )


def _render_html(payload: dict[str, Any]) -> str:
    sections = []
    for mate, result in payload["results"].items():
        cycle_values = [row["mean_quality"] for row in result["per_cycle"]]
        cards = "".join(
            (
                _metric_card("Reads", f"{result['reads']:,}"),
                _metric_card("Bases", f"{result['bases']:,}"),
                _metric_card("GC", result["gc_percent"], "%"),
                _metric_card("Q30 bases", result["q30_percent"], "%"),
                _metric_card("Mean quality", result["mean_quality"]),
                _metric_card("Mean length", result["mean_length"], " bp"),
                _metric_card("Read N50", result["read_length_n50"], " bp"),
                _metric_card(
                    "Estimated duplicates",
                    result["duplicate_estimate_percent"],
                    "%",
                ),
            )
        )
        adapter_rows = "".join(
            f"<tr><td>{html.escape(name)}</td><td>{count:,}</td></tr>"
            for name, count in result["adapter_hits"].items()
        )
        sections.append(
            f'<section><div class="section-heading"><div><p class="eyebrow">{html.escape(mate)}</p>'
            f"<h2>{html.escape(result['source'])}</h2></div>"
            f'<span class="encoding">Phred+{result["phred_offset"]}</span></div>'
            f'<div class="metrics">{cards}</div>'
            '<div class="panel"><h3>Mean quality by cycle</h3>'
            f"{_line_chart(cycle_values)}</div>"
            '<div class="panel"><h3>Adapter observations</h3>'
            f"<table><thead><tr><th>Adapter signature</th><th>Reads</th></tr></thead>"
            f"<tbody>{adapter_rows}</tbody></table></div></section>"
        )

    sample = html.escape(str(payload.get("sample", "NGS sample")))
    generated = html.escape(str(payload.get("generated_at", "")))
    version = html.escape(str(payload.get("tool_version", "")))
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>NGS Core QC — {sample}</title>
<style>
:root {{ color-scheme: dark; --bg:#071018; --panel:#101d27; --line:#263947;
  --text:#edf7fb; --muted:#9bb0bc; --accent:#4de2c5; --blue:#68a7ff; }}
* {{ box-sizing:border-box }} body {{ margin:0; background:var(--bg); color:var(--text);
  font:15px/1.5 Inter,ui-sans-serif,system-ui,-apple-system,sans-serif }}
main {{ max-width:1050px; margin:auto; padding:56px 24px 80px }}
header {{ padding-bottom:32px; border-bottom:1px solid var(--line) }}
.eyebrow {{ color:var(--accent); letter-spacing:.16em; text-transform:uppercase;
  font-size:12px; font-weight:700; margin:0 0 6px }} h1 {{ font-size:42px; margin:0 }}
.subtitle,.meta,.empty {{ color:var(--muted) }}
.meta {{ margin-top:18px }} section {{ margin-top:48px }}
.section-heading {{ display:flex; align-items:end; justify-content:space-between; gap:20px }}
h2 {{ margin:0; font-size:24px; word-break:break-all }} h3 {{ margin:0 0 18px; font-size:17px }}
.encoding {{ color:var(--muted); border:1px solid var(--line); padding:5px 10px;
  border-radius:999px }}
.metrics {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr));
  gap:12px; margin:22px 0 }}
.metric,.panel {{ background:var(--panel); border:1px solid var(--line); border-radius:12px }}
.metric {{ padding:16px }} .metric span {{ display:block; color:var(--muted); font-size:12px }}
.metric strong {{ display:block; margin-top:5px; font-size:22px }}
.panel {{ padding:22px; margin-top:12px }}
svg {{ width:100%; height:auto }} .grid {{ stroke:var(--line); stroke-width:1 }}
.axis {{ fill:var(--muted); font-size:11px }} .quality-line {{ fill:none; stroke:var(--accent);
  stroke-width:3; stroke-linecap:round; stroke-linejoin:round }}
table {{ width:100%; border-collapse:collapse }} th,td {{ text-align:left; padding:10px 8px;
  border-bottom:1px solid var(--line) }} th {{ color:var(--muted); font-size:12px }}
footer {{ color:var(--muted); margin-top:48px; text-align:center }}
@media(max-width:600px) {{
  h1 {{ font-size:32px }} .section-heading {{ display:block }}
  .encoding {{ display:inline-block; margin-top:12px }}
}}
</style>
</head>
<body><main>
<header><p class="eyebrow">NGS Core quality control</p><h1>{sample}</h1>
<p class="subtitle">Streaming FASTQ metrics for reproducible sequencing assessment.</p>
<p class="meta">Generated {generated} · NGS Core {version}</p></header>
{"".join(sections)}
<footer>Generated locally by NGS Core. No sequence data leaves the analysis environment.</footer>
</main></body></html>
"""
