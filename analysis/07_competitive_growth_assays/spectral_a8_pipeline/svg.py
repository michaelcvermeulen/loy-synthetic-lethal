from __future__ import annotations

from html import escape
from pathlib import Path

import numpy as np
import pandas as pd


def write_scatter_density_svg(
    events: pd.DataFrame,
    polygon: np.ndarray,
    output: Path,
    title: str,
    bins: int = 80,
) -> None:
    x = events["FSC-A"].to_numpy(dtype=float)
    y = events["SSC (Violet)-A"].to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    lx, ly = np.log10(x[valid]), np.log10(y[valid])
    xlo, xhi = np.quantile(lx, [0.001, 0.999])
    ylo, yhi = np.quantile(ly, [0.001, 0.999])
    density, xedges, yedges = np.histogram2d(lx, ly, bins=bins, range=[[xlo, xhi], [ylo, yhi]])
    transformed = np.log1p(density)
    maximum = float(transformed.max()) or 1.0

    width, height = 760, 620
    left, top, plot_width, plot_height = 75, 55, 620, 500

    def px(value: float) -> float:
        return left + (value - xlo) / (xhi - xlo) * plot_width

    def py(value: float) -> float:
        return top + plot_height - (value - ylo) / (yhi - ylo) * plot_height

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="18">{escape(title)}</text>',
    ]
    cell_width = plot_width / bins + 0.2
    cell_height = plot_height / bins + 0.2
    for ix, iy in np.argwhere(density > 0):
        intensity = transformed[ix, iy] / maximum
        blue = int(245 - 180 * intensity)
        red = int(245 - 220 * intensity)
        color = f"rgb({red},{blue},245)"
        parts.append(
            f'<rect x="{px(xedges[ix]):.2f}" y="{py(yedges[iy + 1]):.2f}" '
            f'width="{cell_width:.2f}" height="{cell_height:.2f}" fill="{color}"/>'
        )
    vertices = np.asarray(polygon, dtype=float)
    polygon_points = " ".join(
        f"{px(np.log10(vx)):.2f},{py(np.log10(vy)):.2f}" for vx, vy in vertices if vx > 0 and vy > 0
    )
    parts.extend(
        [
            f'<polygon points="{polygon_points}" fill="none" stroke="#d62728" stroke-width="3"/>',
            f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#222"/>',
            f'<text x="{left + plot_width / 2}" y="600" text-anchor="middle" font-family="sans-serif">log10 FSC-A</text>',
            f'<text x="18" y="{top + plot_height / 2}" transform="rotate(-90 18 {top + plot_height / 2})" text-anchor="middle" font-family="sans-serif">log10 SSC (Violet)-A</text>',
            "</svg>",
        ]
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(parts), encoding="utf-8")


_FONT = 'font-family="sans-serif"'


def _frame(width: int, height: int, title: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="26" text-anchor="middle" {_FONT} font-size="17">{escape(title)}</text>',
    ]


def _write(parts: list[str], output: Path) -> None:
    parts.append("</svg>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(parts), encoding="utf-8")


def write_reporter_density_svg(
    values: np.ndarray,
    threshold: float,
    output: Path,
    title: str,
    bins: int = 90,
    secondary_threshold: float | None = None,
) -> None:
    """Histogram of transformed reporter values with the accepted threshold."""
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    counts, edges = np.histogram(data, bins=bins)
    lo, hi = float(edges[0]), float(edges[-1])
    peak = float(counts.max()) or 1.0

    width, height = 760, 420
    left, top, plot_width, plot_height = 70, 50, 650, 300

    def px(value: float) -> float:
        return left + (value - lo) / (hi - lo) * plot_width

    def py(count: float) -> float:
        return top + plot_height - (count / peak) * plot_height

    parts = _frame(width, height, title)
    step = plot_width / len(counts)
    for index, count in enumerate(counts):
        if count <= 0:
            continue
        x = left + index * step
        y = py(float(count))
        parts.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{step + 0.4:.2f}" '
            f'height="{top + plot_height - y:.2f}" fill="#6b8fd4"/>'
        )
    for value, color in ((threshold, "#d62728"), (secondary_threshold, "#2ca02c")):
        if value is None or not (lo <= value <= hi):
            continue
        parts.append(
            f'<line x1="{px(value):.2f}" y1="{top}" x2="{px(value):.2f}" '
            f'y2="{top + plot_height}" stroke="{color}" stroke-width="2.5" stroke-dasharray="7,4"/>'
        )
    parts.append(
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#222"/>'
    )
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        value = lo + fraction * (hi - lo)
        parts.append(
            f'<text x="{px(value):.2f}" y="{top + plot_height + 18}" text-anchor="middle" '
            f'{_FONT} font-size="11">{value:.1f}</text>'
        )
    parts.append(
        f'<text x="{left + plot_width / 2}" y="{top + plot_height + 42}" text-anchor="middle" '
        f'{_FONT} font-size="12">arcsinh transformed value</text>'
    )
    parts.append(
        f'<text x="{left}" y="{height - 18}" {_FONT} font-size="11">'
        f'dashed red: accepted threshold at {threshold:.3f}</text>'
    )
    _write(parts, output)


def write_time_qc_svg(bins: pd.DataFrame, output: Path, title: str) -> None:
    """Binned event rate over acquisition order, with flagged bins highlighted."""
    frame = bins.reset_index(drop=True)
    rates = frame["event_rate"].to_numpy(dtype=float)
    finite = rates[np.isfinite(rates)]
    top_rate = float(finite.max()) if finite.size else 1.0

    width, height = 820, 380
    left, top, plot_width, plot_height = 70, 50, 700, 260
    step = plot_width / max(len(frame), 1)

    parts = _frame(width, height, title)
    segments = list(dict.fromkeys(frame["segment"].tolist()))
    for index, row in frame.iterrows():
        rate = float(row["event_rate"])
        if not np.isfinite(rate):
            rate = top_rate
        bar = (rate / top_rate) * plot_height if top_rate > 0 else 0.0
        color = "#d62728" if bool(row["bad"]) else "#4c78a8"
        parts.append(
            f'<rect x="{left + index * step:.2f}" y="{top + plot_height - bar:.2f}" '
            f'width="{max(step - 1.0, 1.0):.2f}" height="{bar:.2f}" fill="{color}"/>'
        )
    boundary = 0
    for segment in segments[:-1]:
        boundary += int((frame["segment"] == segment).sum())
        x = left + boundary * step
        parts.append(
            f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_height}" '
            f'stroke="#333" stroke-width="1.5" stroke-dasharray="4,3"/>'
        )
    parts.append(
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#222"/>'
    )
    flagged = int(frame["bad"].sum())
    parts.append(
        f'<text x="{left}" y="{top + plot_height + 24}" {_FONT} font-size="12">'
        f'{len(frame)} bins, {flagged} flagged (red); dashed lines separate acquisition segments</text>'
    )
    parts.append(
        f'<text x="18" y="{top + plot_height / 2}" transform="rotate(-90 18 {top + plot_height / 2})" '
        f'text-anchor="middle" {_FONT} font-size="12">events per second</text>'
    )
    _write(parts, output)


def write_saturation_composition_svg(saturation: pd.DataFrame, output: Path, title: str) -> None:
    """Grouped bars comparing saturation categories per recording."""
    frame = saturation.reset_index(drop=True)
    categories = [column for column in frame.columns if column != "sample"]
    colors = {
        "any_saturated": "#8c8c8c",
        "spectral_saturated": "#d62728",
        "imgb3_saturated": "#ff9d5c",
    }

    width, height = 760, 400
    left, top, plot_width, plot_height = 70, 50, 640, 260
    group_width = plot_width / max(len(frame), 1)
    bar_width = group_width / (len(categories) + 1)
    top_value = float(np.nanmax(frame[categories].to_numpy(dtype=float))) or 1.0

    parts = _frame(width, height, title)
    for row_index, row in frame.iterrows():
        for category_index, category in enumerate(categories):
            value = float(row[category])
            bar = (value / top_value) * plot_height
            x = left + row_index * group_width + category_index * bar_width + bar_width / 2
            fill = colors.get(category, "#4c78a8")
            parts.append(
                f'<rect x="{x:.2f}" y="{top + plot_height - bar:.2f}" width="{bar_width * 0.9:.2f}" '
                f'height="{bar:.2f}" fill="{fill}"/>'
            )
            parts.append(
                f'<text x="{x + bar_width * 0.45:.2f}" y="{top + plot_height - bar - 5:.2f}" '
                f'text-anchor="middle" {_FONT} font-size="10">{value:.2f}</text>'
            )
        parts.append(
            f'<text x="{left + row_index * group_width + group_width / 2:.2f}" '
            f'y="{top + plot_height + 18}" text-anchor="middle" {_FONT} font-size="12">'
            f'{escape(str(row["sample"]))}</text>'
        )
    parts.append(
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#222"/>'
    )
    for index, category in enumerate(categories):
        fill = colors.get(category, "#4c78a8")
        parts.append(
            f'<rect x="{left + index * 190}" y="{height - 34}" width="12" height="12" fill="{fill}"/>'
        )
        parts.append(
            f'<text x="{left + index * 190 + 18}" y="{height - 24}" {_FONT} font-size="11">'
            f'{escape(category)}</text>'
        )
    parts.append(
        f'<text x="18" y="{top + plot_height / 2}" transform="rotate(-90 18 {top + plot_height / 2})" '
        f'text-anchor="middle" {_FONT} font-size="12">percent of events</text>'
    )
    _write(parts, output)


def write_sensitivity_delta_svg(summary: pd.DataFrame, output: Path, title: str) -> None:
    """Horizontal bars of the maximum absolute change contributed by each family."""
    frame = summary.reset_index(drop=True)
    values = frame["max_abs_delta_loy_fraction"].to_numpy(dtype=float)
    top_value = float(np.nanmax(values)) if values.size else 1.0
    top_value = top_value if top_value > 0 else 1.0

    width = 760
    row_height = 30
    left, top, plot_width = 150, 50, 480
    height = top + len(frame) * row_height + 70

    parts = _frame(width, height, title)
    for index, row in frame.iterrows():
        value = float(row["max_abs_delta_loy_fraction"])
        bar = (value / top_value) * plot_width
        y = top + index * row_height
        parts.append(
            f'<rect x="{left}" y="{y:.2f}" width="{bar:.2f}" height="{row_height - 9}" fill="#4c78a8"/>'
        )
        parts.append(
            f'<text x="{left - 8}" y="{y + row_height - 14}" text-anchor="end" {_FONT} font-size="12">'
            f'{escape(str(row["family"]))}</text>'
        )
        log2_delta = float(row["max_abs_delta_log2_ratio"])
        parts.append(
            f'<text x="{left + bar + 6:.2f}" y="{y + row_height - 14}" {_FONT} font-size="11">'
            f'{value:.4f} (log2 {log2_delta:.3f})</text>'
        )
    parts.append(
        f'<text x="{left}" y="{height - 28}" {_FONT} font-size="12">'
        f'maximum absolute change in LOY fraction versus the baseline analysis</text>'
    )
    _write(parts, output)


def write_pulse_relationship_svg(
    area: np.ndarray,
    height_values: np.ndarray,
    fit,
    mad_limit: float,
    output: Path,
    title: str,
    sample_limit: int = 6000,
) -> None:
    """Log-log pulse area against height with the robust fit and residual band."""
    a = np.asarray(area, dtype=float)
    h = np.asarray(height_values, dtype=float)
    valid = np.isfinite(a) & np.isfinite(h) & (a > 0) & (h > 0)
    x = np.log10(a[valid])
    y = np.log10(h[valid])
    if x.size > sample_limit:
        step = max(1, x.size // sample_limit)
        x, y = x[::step], y[::step]

    xlo, xhi = float(np.quantile(x, 0.001)), float(np.quantile(x, 0.999))
    ylo, yhi = float(np.quantile(y, 0.001)), float(np.quantile(y, 0.999))

    width, svg_height = 720, 560
    left, top, plot_width, plot_height = 75, 50, 590, 440

    def px(value: float) -> float:
        return left + (value - xlo) / (xhi - xlo) * plot_width

    def py(value: float) -> float:
        return top + plot_height - (value - ylo) / (yhi - ylo) * plot_height

    parts = _frame(width, svg_height, title)
    for xv, yv in zip(x, y):
        if xlo <= xv <= xhi and ylo <= yv <= yhi:
            parts.append(
                f'<circle cx="{px(xv):.1f}" cy="{py(yv):.1f}" r="0.7" fill="#3b6ea5" fill-opacity="0.35"/>'
            )
    tolerance = mad_limit * fit.residual_scale
    for offset, dash in (
        (0.0, ""),
        (fit.residual_center + tolerance, ' stroke-dasharray="6,4"'),
        (fit.residual_center - tolerance, ' stroke-dasharray="6,4"'),
    ):
        y1 = fit.intercept + fit.slope * xlo + offset
        y2 = fit.intercept + fit.slope * xhi + offset
        parts.append(
            f'<line x1="{px(xlo):.2f}" y1="{py(y1):.2f}" x2="{px(xhi):.2f}" y2="{py(y2):.2f}" '
            f'stroke="#d62728" stroke-width="2"{dash}/>'
        )
    parts.append(
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#222"/>'
    )
    parts.append(
        f'<text x="{left}" y="{svg_height - 16}" {_FONT} font-size="11">'
        f'slope {fit.slope:.4f}, residual scale {fit.residual_scale:.4f}, band {mad_limit:g} MAD</text>'
    )
    _write(parts, output)
