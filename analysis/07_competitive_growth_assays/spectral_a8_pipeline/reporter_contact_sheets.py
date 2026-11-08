"""GFP-versus-mCherry contact sheets for every Run 1 and Run 2 sample."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import math
from pathlib import Path
import textwrap

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from .competition_report import _file_map, _primary_features, fit_lda, load_cef_gates
from .fcs_io import load_fcs
from .longitudinal_report import build_new_run_design


@dataclass(frozen=True)
class ReporterPlotSample:
    path: Path
    sheet: str
    timepoint: str
    role: str
    cef_path: Path


def format_panel_title(stem: str) -> str:
    """Wrap a sample identifier without allowing it to cross panel boundaries."""
    readable = stem.replace("_CTRL_", " / CTRL_")
    return "\n".join(textwrap.wrap(readable, width=34)[:2])


def inventory_reporter_samples(repo_root: Path) -> list[ReporterPlotSample]:
    """Assign every longitudinal-study FCS file to exactly one contact sheet."""
    repo = Path(repo_root).resolve()
    d0_root = repo / "spectral" / "20260828_113753_Export" / "Mike_QC1"
    d0_cef = d0_root / "20260828_113753_Mike_QC1.cef"
    d4_root = repo / "spectral" / "20260901_125230_Mike_Run1" / "Mike_Run1"
    d4_cef = d4_root / "20260901_125230_Mike_Run1.cef"
    new_root = repo / "spectral" / "Mike_Run 2"
    new_cef = new_root / "20260904_171431_Mike_Run 2.cef"
    samples = [
        ReporterPlotSample(path, "run1_d0", "D0", "d0", d0_cef)
        for path in sorted(d0_root.rglob("*.fcs"))
    ]
    for path in sorted((d4_root / "1").glob("*.fcs")):
        prefix = path.name[0].upper()
        if prefix == "A":
            sheet, role = "run1_d4_controls", "control"
        elif prefix in {"B", "C"}:
            sheet, role = "run1_d4_ratio_standards", "ratio_standard"
        else:
            sheet, role = "run1_d4_experimental", "experimental"
        samples.append(ReporterPlotSample(path, sheet, "D4", role, d4_cef))
    for parsed in build_new_run_design((new_root / "run2").glob("*.fcs")):
        sheet = "run1_d7" if parsed.experiment_run == "run1" else "run2_d0"
        samples.append(
            ReporterPlotSample(parsed.path, sheet, parsed.timepoint, "experimental", new_cef)
        )
    return samples


def generate_contact_sheets(repo_root: Path, output: Path) -> list[Path]:
    """Render six shared-axis contact sheets and a complete plot manifest."""
    repo = Path(repo_root).resolve()
    destination = Path(output)
    if destination.exists():
        raise FileExistsError(f"output already exists: {destination}")
    destination.mkdir(parents=True)

    samples = inventory_reporter_samples(repo)
    d4_root = repo / "spectral" / "20260901_125230_Mike_Run1" / "Mike_Run1"
    d4_files = _file_map(d4_root / "1")
    d4_gates = load_cef_gates(d4_root / "20260901_125230_Mike_Run1.cef")
    control_roles = {
        "A03-CTRL_GFP": 1,
        "A05-LOY_GFP": 1,
        "A07-CTRL_cherry": 0,
        "A09-LOY_cherry": 0,
    }
    training_features = []
    training_labels = []
    for name, label in control_roles.items():
        values, count, _saturated = _primary_features(d4_files[name], d4_gates)
        training_features.append(values)
        training_labels.append(np.full(count, label, dtype=int))
    classifier = fit_lda(np.vstack(training_features), np.concatenate(training_labels))

    gates_by_cef = {}
    plotted = []
    manifest = []
    axis_values = []
    for sample in samples:
        loaded = load_fcs(sample.path)
        if not {"eGFP-A", "mCherry-A"}.issubset(loaded.events.columns):
            plotted.append((sample, None, None))
            manifest.append(
                {
                    "sheet": sample.sheet,
                    "timepoint": sample.timepoint,
                    "role": sample.role,
                    "source_path": str(sample.path.relative_to(repo)),
                    "status": "reporter_channels_unavailable",
                    "primary_singlets": "",
                    "gfp_fraction": "",
                    "mcherry_fraction": "",
                }
            )
            continue
        if sample.cef_path not in gates_by_cef:
            gates_by_cef[sample.cef_path] = load_cef_gates(sample.cef_path)
        values, count, _saturated = _primary_features(
            sample.path, gates_by_cef[sample.cef_path]
        )
        predictions = classifier.predict(values)
        unstained = "unstained" in sample.path.name.lower()
        if unstained:
            gfp_fraction = mcherry_fraction = 0.0
            classes = np.full(count, -1, dtype=int)
        else:
            gfp_fraction = float(predictions.mean())
            mcherry_fraction = 1.0 - gfp_fraction
            classes = predictions
        plotted.append((sample, values, classes))
        axis_values.append(values)
        manifest.append(
            {
                "sheet": sample.sheet,
                "timepoint": sample.timepoint,
                "role": sample.role,
                "source_path": str(sample.path.relative_to(repo)),
                "status": "plotted",
                "primary_singlets": count,
                "gfp_fraction": gfp_fraction,
                "mcherry_fraction": mcherry_fraction,
            }
        )

    combined = np.vstack(axis_values)
    x_low, y_low = np.quantile(combined, 0.002, axis=0)
    x_high, y_high = np.quantile(combined, 0.998, axis=0)
    x_pad = 0.05 * (x_high - x_low)
    y_pad = 0.05 * (y_high - y_low)
    limits = (x_low - x_pad, x_high + x_pad, y_low - y_pad, y_high + y_pad)
    manifest_by_path = {row["source_path"]: row for row in manifest}

    sheet_specs = [
        ("run1_d0", "01_run1_d0.png", "Run 1 — D0", 3),
        ("run1_d4_controls", "02_run1_d4_controls.png", "Run 1 — D4 controls", 3),
        (
            "run1_d4_ratio_standards",
            "03_run1_d4_ratio_standards.png",
            "Run 1 — D4 ratio standards",
            5,
        ),
        (
            "run1_d4_experimental",
            "04_run1_d4_experimental.png",
            "Run 1 — D4 experimental wells",
            4,
        ),
        ("run1_d7", "05_run1_d7.png", "Run 1 — D7", 4),
        ("run2_d0", "06_run2_d0.png", "Run 2 — D0", 4),
    ]
    generated = []
    for sheet, filename, title, columns in sheet_specs:
        members = [row for row in plotted if row[0].sheet == sheet]
        rows = math.ceil(len(members) / columns)
        panel_width, panel_height = 640, 430
        header_height, footer_height = 90, 58
        image = Image.new(
            "RGBA",
            (panel_width * columns, header_height + panel_height * rows + footer_height),
            "white",
        )
        draw = ImageDraw.Draw(image, "RGBA")
        title_font = ImageFont.truetype("arialbd.ttf", 25)
        panel_title_font = ImageFont.truetype("arialbd.ttf", 14)
        label_font = ImageFont.truetype("arial.ttf", 13)
        small_font = ImageFont.truetype("arial.ttf", 11)
        draw.text(
            (image.width / 2, 26),
            f"{title}: primary singlets, shared reporter-space axes",
            fill="#111827",
            font=title_font,
            anchor="ma",
        )
        for member_index, (sample, values, classes) in enumerate(members):
            grid_row, grid_column = divmod(member_index, columns)
            panel_left = grid_column * panel_width
            panel_top = header_height + grid_row * panel_height
            plot_left, plot_top = panel_left + 72, panel_top + 62
            plot_right, plot_bottom = panel_left + panel_width - 22, panel_top + panel_height - 62
            plot_width, plot_height = plot_right - plot_left, plot_bottom - plot_top
            draw.rectangle(
                (plot_left, plot_top, plot_right, plot_bottom),
                fill="white",
                outline="#9ca3af",
                width=1,
            )
            for tick_index in range(5):
                fraction = tick_index / 4
                tick_x = plot_left + fraction * plot_width
                tick_y = plot_bottom - fraction * plot_height
                draw.line((tick_x, plot_top, tick_x, plot_bottom), fill="#e5e7eb", width=1)
                draw.line((plot_left, tick_y, plot_right, tick_y), fill="#e5e7eb", width=1)
                x_value = limits[0] + fraction * (limits[1] - limits[0])
                y_value = limits[2] + fraction * (limits[3] - limits[2])
                draw.text(
                    (tick_x, plot_bottom + 5),
                    f"{x_value:.1f}",
                    fill="#4b5563",
                    font=small_font,
                    anchor="ma",
                )
                draw.text(
                    (plot_left - 7, tick_y),
                    f"{y_value:.1f}",
                    fill="#4b5563",
                    font=small_font,
                    anchor="rm",
                )
            if values is None:
                draw.multiline_text(
                    ((plot_left + plot_right) / 2, (plot_top + plot_bottom) / 2 - 18),
                    "Reporter-space channels\nunavailable",
                    fill="#7f1d1d",
                    font=panel_title_font,
                    anchor="mm",
                    align="center",
                    spacing=5,
                )
                draw.text(
                    ((plot_left + plot_right) / 2, (plot_top + plot_bottom) / 2 + 32),
                    "No eGFP-A + mCherry-A pair",
                    fill="#6b7280",
                    font=small_font,
                    anchor="ma",
                )
            else:
                rng = np.random.default_rng(20260904)
                selected = (
                    rng.choice(len(values), 3000, replace=False)
                    if len(values) > 3000
                    else np.arange(len(values))
                )
                shown = values[selected]
                shown_classes = classes[selected]
                plot_layer = Image.new("RGBA", (plot_width + 1, plot_height + 1), (0, 0, 0, 0))
                plot_draw = ImageDraw.Draw(plot_layer, "RGBA")
                for class_value, colour in ((0, "#be185d"), (1, "#15803d"), (-1, "#6b7280")):
                    keep = shown_classes == class_value
                    if np.any(keep):
                        x_pixels = (
                            (shown[keep, 0] - limits[0])
                            / (limits[1] - limits[0])
                            * plot_width
                        )
                        y_pixels = (
                            plot_height
                            - (shown[keep, 1] - limits[2])
                            / (limits[3] - limits[2])
                            * plot_height
                        )
                        points = [
                            (int(x), int(y))
                            for x, y in zip(x_pixels, y_pixels)
                            if 0 <= x <= plot_width and 0 <= y <= plot_height
                        ]
                        plot_draw.point(points, fill=colour + "70")
                if abs(classifier.weights[1]) > 1e-12:
                    boundary_x = np.array([limits[0], limits[1]])
                    boundary_y = -(
                        classifier.weights[0] * boundary_x + classifier.intercept
                    ) / classifier.weights[1]
                    boundary_pixels = [
                        (
                            int((x - limits[0]) / (limits[1] - limits[0]) * plot_width),
                            int(
                                plot_height
                                - (y - limits[2]) / (limits[3] - limits[2]) * plot_height
                            ),
                        )
                        for x, y in zip(boundary_x, boundary_y)
                    ]
                    plot_draw.line(
                        boundary_pixels,
                        fill="#111827",
                        width=2,
                    )
                image.alpha_composite(plot_layer, (plot_left, plot_top))
                row = manifest_by_path[str(sample.path.relative_to(repo))]
                if "unstained" in sample.path.name.lower():
                    annotation = f"n={row['primary_singlets']:,} · unstained background"
                else:
                    annotation = (
                        f"n={row['primary_singlets']:,} · GFP {100 * row['gfp_fraction']:.1f}% · "
                        f"mCherry {100 * row['mcherry_fraction']:.1f}%"
                    )
                annotation_box = draw.textbbox(
                    (plot_left + 5, plot_bottom - 5),
                    annotation,
                    font=small_font,
                    anchor="ld",
                )
                draw.rectangle(annotation_box, fill=(255, 255, 255, 205))
                draw.text(
                    (plot_left + 5, plot_bottom - 5),
                    annotation,
                    fill="#111827",
                    font=small_font,
                    anchor="ld",
                )
            wrapped_title = format_panel_title(sample.path.stem)
            draw.multiline_text(
                (panel_left + panel_width / 2, panel_top + 8),
                wrapped_title,
                fill="#111827",
                font=panel_title_font,
                anchor="ma",
                align="center",
                spacing=1,
            )
            draw.text(
                ((plot_left + plot_right) / 2, panel_top + panel_height - 27),
                "arcsinh(eGFP-A / 1000)",
                fill="#374151",
                font=label_font,
                anchor="ma",
            )
            vertical_label = Image.new("RGBA", (plot_height, 26), (0, 0, 0, 0))
            vertical_draw = ImageDraw.Draw(vertical_label)
            vertical_draw.text(
                (plot_height / 2, 13),
                "arcsinh(mCherry-A / 1000)",
                fill="#374151",
                font=label_font,
                anchor="mm",
            )
            vertical_label = vertical_label.rotate(90, expand=True)
            image.alpha_composite(vertical_label, (panel_left + 5, plot_top))
        draw.text(
            (image.width / 2, image.height - 27),
            "green = GFP-class  ·  magenta = mCherry-class  ·  black = shared LDA boundary",
            fill="#4b5563",
            font=label_font,
            anchor="mm",
        )
        path = destination / filename
        image.convert("RGB").save(path, format="PNG", optimize=True)
        generated.append(path)

    with (destination / "plot_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)
    return generated


def main(argv: list[str] | None = None) -> int:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Plot GFP versus mCherry for all Run 1/2 samples")
    parser.add_argument(
        "--output",
        type=Path,
        default=repo / "analysis" / "spectral_competition_longitudinal_20260904" / "reporter_contact_sheets",
    )
    args = parser.parse_args(argv)
    generated = generate_contact_sheets(repo, args.output.resolve())
    for path in generated:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
