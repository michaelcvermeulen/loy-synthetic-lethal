"""Render the analysis report as Markdown and HTML without external converters."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

TEMPLATE_DIR = Path(__file__).parent / "templates"


@dataclass(frozen=True)
class ReportPaths:
    markdown: Path
    html: Path


def _format_cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if value != value:
            return "n/a"
        return f"{value:.6g}"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value)


def render_markdown_table(frame: pd.DataFrame) -> str:
    """Render a DataFrame as a GitHub-flavoured Markdown table."""
    if frame.empty:
        return "_No rows._"
    columns = [str(column) for column in frame.columns]
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join("---" for _ in columns) + "|"]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(_format_cell(row[column]) for column in frame.columns) + " |")
    return "\n".join(lines)


def build_markdown(data: Mapping[str, Any]) -> str:
    """Assemble the Markdown report body."""
    parts: list[str] = [f"# {data.get('title', 'A8 competition analysis')}", ""]

    if data.get("summary"):
        parts.extend([str(data["summary"]), ""])

    interval_note = data.get("interval_note")
    if interval_note:
        parts.extend(["## Interpretation of intervals", "", str(interval_note), ""])

    tables: Mapping[str, pd.DataFrame] = data.get("tables") or {}
    if tables:
        parts.extend(["## Tables", ""])
        for name, frame in tables.items():
            parts.extend([f"### {name}", "", render_markdown_table(frame), ""])

    figures: Sequence[Mapping[str, str]] = data.get("figures") or []
    if figures:
        parts.extend(["## Figures", ""])
        for figure in figures:
            caption = figure.get("caption", "")
            path = figure.get("path", "")
            parts.extend([f"![{caption}]({path})", "", f"*{caption}* — `{path}`", ""])

    limitations: Sequence[str] = data.get("limitations") or []
    if limitations:
        parts.extend(["## Limitations", ""])
        parts.extend(f"- {item}" for item in limitations)
        parts.append("")

    return "\n".join(parts)


def write_report(data: Mapping[str, Any], output: Path) -> ReportPaths:
    """Write ``report.md`` and ``report.html`` into ``output``."""
    destination = Path(output)
    destination.mkdir(parents=True, exist_ok=True)

    markdown = build_markdown(data)
    markdown_path = destination / "report.md"
    markdown_path.write_text(markdown, encoding="utf-8")

    environment = Environment(
        loader=FileSystemLoader(str(TEMPLATE_DIR)),
        autoescape=select_autoescape(["html"]),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    template = environment.get_template("report.html.j2")
    tables = {
        name: frame.to_html(index=False, border=0, na_rep="n/a", float_format=lambda v: f"{v:.6g}")
        for name, frame in (data.get("tables") or {}).items()
    }
    html = template.render(
        title=data.get("title", "A8 competition analysis"),
        summary=data.get("summary"),
        interval_note=data.get("interval_note"),
        tables=tables,
        figures=data.get("figures") or [],
        limitations=data.get("limitations") or [],
    )
    html_path = destination / "report.html"
    html_path.write_text(html, encoding="utf-8")

    return ReportPaths(markdown=markdown_path, html=html_path)
