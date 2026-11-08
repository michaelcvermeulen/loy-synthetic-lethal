"""Longitudinal analysis for Run 1 (D0/D4/D7) and Run 2 (D0)."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Iterable, Mapping

import numpy as np

from .competition_report import (
    _file_map,
    _primary_features,
    analyze_repository as analyze_prior_repository,
    delta_log_odds,
    fit_lda,
    load_cef_gates,
    summarise_effects,
)
from .fcs_io import load_fcs


@dataclass(frozen=True)
class LongitudinalSample:
    path: Path
    experiment_run: str
    timepoint: str
    orientation: str
    biological_replicate: int
    technical_replicate: int


_NEW_RUN_PATTERN = re.compile(
    r"^(?P<well>[ABDE])\d{2}-"
    r"(?P<mixture>LOY_GFP_CTRL_CHERRY|LOY_CHERRY_CTRL_GFP)_"
    r"rep(?P<bio>[1-3])_(?P<run>run[12])(?:_+\d+|_)?$",
    re.IGNORECASE,
)


def parse_new_run_sample(path: Path) -> LongitudinalSample:
    """Parse one Mike_Run 2 FCS filename into experimental metadata."""
    source = Path(path)
    match = _NEW_RUN_PATTERN.match(source.stem)
    if match is None:
        raise ValueError(f"unrecognised Mike_Run 2 sample name: {source.name}")
    well = match.group("well").upper()
    experiment_run = match.group("run").lower()
    expected_letters = {"run1": {"A", "B"}, "run2": {"D", "E"}}
    if well not in expected_letters[experiment_run]:
        raise ValueError(
            f"well {well} is inconsistent with {experiment_run}: {source.name}"
        )
    mixture = match.group("mixture").upper()
    orientation = "LOY_GFP" if mixture == "LOY_GFP_CTRL_CHERRY" else "LOY_mCherry"
    return LongitudinalSample(
        path=source,
        experiment_run=experiment_run,
        timepoint="D7" if experiment_run == "run1" else "D0",
        orientation=orientation,
        biological_replicate=int(match.group("bio")),
        technical_replicate=1 if well in {"A", "D"} else 2,
    )


def build_new_run_design(paths: Iterable[Path]) -> list[LongitudinalSample]:
    """Build and validate the two-technical-replicate sample design."""
    samples = sorted(
        (parse_new_run_sample(path) for path in paths), key=lambda sample: sample.path.name
    )
    groups: dict[tuple[str, str, int], list[int]] = {}
    for sample in samples:
        key = (
            sample.experiment_run,
            sample.orientation,
            sample.biological_replicate,
        )
        groups.setdefault(key, []).append(sample.technical_replicate)
    for key, replicates in groups.items():
        if sorted(replicates) != [1, 2]:
            raise ValueError(
                "expected technical replicates 1 and 2 for "
                f"{key[0]} {key[1]} biological replicate {key[2]}; "
                f"found {sorted(replicates)}"
            )
    return samples


def pool_new_run_measurements(
    samples: Iterable[LongitudinalSample],
    measurements: Mapping[Path, Mapping[str, float | int]],
) -> list[dict]:
    """Pool each biological replicate's two technical wells by event counts."""
    groups: dict[tuple[str, str, str, int], list[LongitudinalSample]] = {}
    for sample in samples:
        key = (
            sample.experiment_run,
            sample.timepoint,
            sample.orientation,
            sample.biological_replicate,
        )
        groups.setdefault(key, []).append(sample)
    rows: list[dict] = []
    for key, wells in sorted(groups.items()):
        well_rows = []
        for sample in sorted(wells, key=lambda item: item.technical_replicate):
            measurement = measurements[sample.path]
            count = int(measurement["primary_singlets"])
            loy_count = int(measurement["loy_count"])
            well_rows.append(
                {
                    "sample": sample.path.stem,
                    "technical_replicate": sample.technical_replicate,
                    "loy_count": loy_count,
                    "primary_singlets": count,
                    "loy_fraction": loy_count / count,
                    "saturated_fraction": float(measurement["saturated_fraction"]),
                }
            )
        total = sum(row["primary_singlets"] for row in well_rows)
        positive = sum(row["loy_count"] for row in well_rows)
        rows.append(
            {
                "experiment_run": key[0],
                "timepoint": key[1],
                "orientation": key[2],
                "biological_replicate": key[3],
                "loy_count": positive,
                "primary_singlets": total,
                "loy_fraction": positive / total,
                "technical_difference_pp": 100.0
                * abs(well_rows[0]["loy_fraction"] - well_rows[1]["loy_fraction"]),
                "saturated_fraction": sum(
                    row["saturated_fraction"] * row["primary_singlets"]
                    for row in well_rows
                )
                / total,
                "technical_wells": well_rows,
            }
        )
    return rows


def _summarise_interval(deltas: Mapping[str, list[float]], days: int) -> dict:
    if all(len(values) >= 2 for values in deltas.values()):
        summary = summarise_effects(deltas)
    else:
        common = float(np.mean([np.mean(values) for values in deltas.values()]))
        summary = {
            "common_delta": common,
            "interval_odds_ratio": float(np.exp(common)),
        }
    summary["days"] = days
    summary["per_day_log_ratio"] = summary["common_delta"] / days
    return summary


def build_longitudinal_result(old_result: Mapping, new_rows: Iterable[Mapping]) -> dict:
    """Combine the established Run 1 D0/D4 result with new D7 and Run 2 D0 rows."""
    pooled = [dict(row) for row in new_rows]
    lookup = {
        (row["experiment_run"], row["orientation"], row["biological_replicate"]): row
        for row in pooled
    }
    trajectories = []
    d0_d7: dict[str, list[float]] = {"LOY_GFP": [], "LOY_mCherry": []}
    d4_d7: dict[str, list[float]] = {"LOY_GFP": [], "LOY_mCherry": []}
    for orientation, orientation_data in old_result["orientations"].items():
        d0 = float(orientation_data["d0_loy"])
        for plate in orientation_data["plates"]:
            replicate = int(plate["replicate"])
            d4 = float(plate["d4_loy"])
            d7_row = lookup[("run1", orientation, replicate)]
            d7 = float(d7_row["loy_fraction"])
            delta_07 = delta_log_odds(d0, d7)
            delta_47 = delta_log_odds(d4, d7)
            d0_d7[orientation].append(delta_07)
            d4_d7[orientation].append(delta_47)
            trajectories.append(
                {
                    "orientation": orientation,
                    "biological_replicate": replicate,
                    "timepoints": {"D0": d0, "D4": d4, "D7": d7},
                    "delta_log_odds_d0_d4": delta_log_odds(d0, d4),
                    "delta_log_odds_d0_d7": delta_07,
                    "delta_log_odds_d4_d7": delta_47,
                    "d7_measurement": d7_row,
                }
            )
    return {
        "title": "Spectral competition assay: Run 1 D0/D4/D7 and Run 2 D0",
        "source_runs": {
            **dict(old_result["source_runs"]),
            "new_acquisition": str(Path("spectral") / "Mike_Run 2"),
        },
        "run1": trajectories,
        "run2_d0": sorted(
            (row for row in pooled if row["experiment_run"] == "run2"),
            key=lambda row: (row["orientation"], row["biological_replicate"]),
        ),
        "effects": {
            "d0_to_d4": dict(old_result["summary"]),
            "d0_to_d7": _summarise_interval(d0_d7, 7),
            "d4_to_d7": _summarise_interval(d4_d7, 3),
        },
    }


def verify_unmixing_compatibility(reference: Mapping, candidate: Mapping) -> None:
    """Require the same named spectral matrix and numerical spillover payload."""
    keys = ("bdspectral matrix name", "spill")
    if any(reference.get(key) != candidate.get(key) for key in keys):
        raise ValueError("spectral unmixing differs between acquisitions")


def analyze_new_acquisition(repo_root: Path) -> dict:
    """Gate, classify, pool, and inventory every Mike_Run 2 sample."""
    repo = Path(repo_root).resolve()
    prior_root = repo / "spectral" / "20260901_125230_Mike_Run1" / "Mike_Run1"
    prior_files = _file_map(prior_root / "1")
    prior_gates = load_cef_gates(prior_root / "20260901_125230_Mike_Run1.cef")
    control_roles = {
        "A03-CTRL_GFP": 1,
        "A05-LOY_GFP": 1,
        "A07-CTRL_cherry": 0,
        "A09-LOY_cherry": 0,
    }
    features = []
    labels = []
    for name, label in control_roles.items():
        values, count, _saturated = _primary_features(prior_files[name], prior_gates)
        features.append(values)
        labels.append(np.full(count, label, dtype=int))
    classifier = fit_lda(np.vstack(features), np.concatenate(labels))

    new_root = repo / "spectral" / "Mike_Run 2"
    fcs_root = new_root / "run2"
    samples = build_new_run_design(fcs_root.glob("*.fcs"))
    gates = load_cef_gates(new_root / "20260904_171431_Mike_Run 2.cef")
    reference_metadata = load_fcs(prior_files["A03-CTRL_GFP"]).metadata
    measurements: dict[Path, dict] = {}
    manifest = []
    for sample in samples:
        loaded = load_fcs(sample.path)
        verify_unmixing_compatibility(reference_metadata, loaded.metadata)
        sample_features, count, saturated = _primary_features(sample.path, gates)
        predictions = classifier.predict(sample_features)
        loy_count = (
            int(predictions.sum())
            if sample.orientation == "LOY_GFP"
            else int((1 - predictions).sum())
        )
        measurements[sample.path] = {
            "loy_count": loy_count,
            "primary_singlets": count,
            "saturated_fraction": saturated,
        }
        digest = sha256()
        with sample.path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        manifest.append(
            {
                "path": str(sample.path.relative_to(repo)),
                "size_bytes": sample.path.stat().st_size,
                "sha256": digest.hexdigest(),
            }
        )
    return {
        "sample_count": len(samples),
        "source_manifest": manifest,
        "unmixing_compatible": True,
        "unmixing_matrix": reference_metadata.get("bdspectral matrix name"),
        "pooled_rows": pool_new_run_measurements(samples, measurements),
    }


def _write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _effect_text(effect: Mapping) -> str:
    text = (
        f"{effect['common_delta']:.3f} log-odds; interval odds ratio "
        f"{effect['interval_odds_ratio']:.3f}"
    )
    if "ci_low" in effect:
        text += (
            f" (95% plate-variance interval {effect['odds_ratio_ci_low']:.3f}–"
            f"{effect['odds_ratio_ci_high']:.3f})"
        )
    return text


def _build_markdown(result: Mapping) -> str:
    effect_07 = result["effects"]["d0_to_d7"]
    effect_47 = result["effects"]["d4_to_d7"]
    new = result["new_acquisition"]
    maximum_difference = max(
        (row["technical_difference_pp"] for row in new["pooled_rows"]), default=0.0
    )
    lines = [
        f"# {result['title']}",
        "",
        "## Executive summary",
        "",
        f"For Run 1, the common D0-to-D7 change was {_effect_text(effect_07)}. "
        f"The paired D4-to-D7 change was {_effect_text(effect_47)}.",
        "",
        "Run 2 currently contributes D0 measurements only, so those values are reported as "
        "baseline biological replicates and are not interpreted as a longitudinal effect.",
        "",
        "## Sample mapping",
        "",
        "- Run 1: original August 28 acquisition = D0; September 1 acquisition = D4; "
        "A/B wells in `Mike_Run 2` = D7.",
        "- Run 2: D/E wells in `Mike_Run 2` = D0.",
        "- `rep1`–`rep3` are biological replicates. Each A/B or D/E pair contains two "
        "technical acquisition replicates pooled by classified event counts.",
        "",
        "## Run 1 longitudinal results",
        "",
        "| Dye orientation | Biological replicate | D0 LOY | D4 LOY | D7 LOY | D0→D7 Δ log-odds | D4→D7 Δ log-odds |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in result["run1"]:
        values = row["timepoints"]
        lines.append(
            f"| {row['orientation']} | {row['biological_replicate']} | "
            f"{100 * values['D0']:.2f}% | {100 * values['D4']:.2f}% | "
            f"{100 * values['D7']:.2f}% | {row['delta_log_odds_d0_d7']:.3f} | "
            f"{row['delta_log_odds_d4_d7']:.3f} |"
        )
    lines.extend(
        [
            "",
            "![Run 1 D0/D4/D7 trajectories and Run 2 baselines](figures/longitudinal_loy_fraction.svg)",
            "",
            "## Run 2 D0 baseline",
            "",
            "| Dye orientation | Biological replicate | LOY fraction | Primary singlets | Technical difference |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in result["run2_d0"]:
        lines.append(
            f"| {row['orientation']} | {row['biological_replicate']} | "
            f"{100 * row['loy_fraction']:.2f}% | {row['primary_singlets']} | "
            f"{row['technical_difference_pp']:.2f} pp |"
        )
    lines.extend(
        [
            "",
            "## Quality control and methods",
            "",
            f"- Loaded {new['sample_count']} FCS files from `spectral/Mike_Run 2`; the "
            "source manifest records size and SHA-256 for every file.",
            f"- Spectral unmixing matched the prior acquisition exactly: "
            f"`{new['unmixing_matrix']}` with an identical spillover payload.",
            f"- Maximum disagreement between new technical wells was {maximum_difference:.2f} percentage points.",
            "- Each acquisition used its own CEF Scatter → SSC Singlets → FSC Singlets gate hierarchy.",
            "- The existing equal-prior GFP-versus-mCherry LDA classifier was retained because the "
            "instrument, reporter parameters, and spectral unmixing payload were identical.",
            "- Technical wells were pooled by summing classified LOY events and primary singlets; "
            "biological replicates were not pooled.",
            "",
            "## Interpretation note",
            "",
            "Run 1 D0 remains a shared starting mixture within each dye orientation, whereas Run 2 "
            "has independently measured D0 biological replicates. Run 2 should not be compared with "
            "an endpoint until its matching later timepoints are acquired.",
            "",
        ]
    )
    return "\n".join(lines)


def _write_longitudinal_svg(result: Mapping, output: Path) -> None:
    width, height = 920, 460
    left, right, top, bottom = 80, 880, 55, 400
    x_positions = {"D0": 150, "D4": 430, "D7": 710}
    y = lambda value: bottom - float(value) * (bottom - top)
    colours = {"LOY_GFP": "#2b6cb0", "LOY_mCherry": "#c53030"}
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img">',
        "<title>LOY fraction through Run 1 and Run 2 baseline</title>",
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="460" y="28" text-anchor="middle" font-family="sans-serif" font-size="18">Run 1 trajectories; Run 2 D0 baselines shown as diamonds</text>',
    ]
    for tick in range(0, 101, 20):
        ypos = y(tick / 100)
        parts.append(f'<line x1="{left}" y1="{ypos:.1f}" x2="{right}" y2="{ypos:.1f}" stroke="#e2e8f0"/>')
        parts.append(f'<text x="68" y="{ypos + 4:.1f}" text-anchor="end" font-family="sans-serif" font-size="11">{tick}%</text>')
    for label, xpos in x_positions.items():
        parts.append(f'<text x="{xpos}" y="425" text-anchor="middle" font-family="sans-serif" font-size="13">{label}</text>')
    for row in result["run1"]:
        points = " ".join(
            f"{x_positions[label]},{y(row['timepoints'][label]):.1f}" for label in ("D0", "D4", "D7")
        )
        colour = colours[row["orientation"]]
        parts.append(f'<polyline points="{points}" fill="none" stroke="{colour}" stroke-width="2" opacity="0.7"/>')
        for label in ("D0", "D4", "D7"):
            parts.append(f'<circle cx="{x_positions[label]}" cy="{y(row["timepoints"][label]):.1f}" r="4" fill="{colour}"/>')
    for row in result["run2_d0"]:
        xpos = x_positions["D0"] + 45 + 10 * (row["biological_replicate"] - 2)
        ypos = y(row["loy_fraction"])
        colour = colours[row["orientation"]]
        parts.append(f'<rect x="{xpos - 4}" y="{ypos - 4:.1f}" width="8" height="8" transform="rotate(45 {xpos} {ypos:.1f})" fill="{colour}" stroke="white"/>')
    parts.extend(
        [
            '<text x="790" y="75" font-family="sans-serif" font-size="12" fill="#2b6cb0">LOY-GFP</text>',
            '<text x="790" y="94" font-family="sans-serif" font-size="12" fill="#c53030">LOY-mCherry</text>',
            "</svg>",
        ]
    )
    output.write_text("\n".join(parts), encoding="utf-8")


def write_longitudinal_outputs(result: Mapping, output: Path) -> None:
    """Write the longitudinal report and its machine-readable companion artifacts."""
    destination = Path(output)
    if destination.exists():
        raise FileExistsError(f"output already exists: {destination}")
    figures = destination / "figures"
    figures.mkdir(parents=True)
    (destination / "results.json").write_text(
        json.dumps(result, indent=2, sort_keys=True), encoding="utf-8"
    )
    (destination / "report.md").write_text(_build_markdown(result), encoding="utf-8")
    trajectory_rows = []
    for row in result["run1"]:
        trajectory_rows.append(
            {
                "orientation": row["orientation"],
                "biological_replicate": row["biological_replicate"],
                **row["timepoints"],
                "delta_log_odds_d0_d4": row["delta_log_odds_d0_d4"],
                "delta_log_odds_d0_d7": row["delta_log_odds_d0_d7"],
                "delta_log_odds_d4_d7": row["delta_log_odds_d4_d7"],
            }
        )
    _write_csv(
        destination / "run1_trajectories.csv",
        trajectory_rows,
        [
            "orientation",
            "biological_replicate",
            "D0",
            "D4",
            "D7",
            "delta_log_odds_d0_d4",
            "delta_log_odds_d0_d7",
            "delta_log_odds_d4_d7",
        ],
    )
    _write_csv(
        destination / "pooled_new_measurements.csv",
        list(result["new_acquisition"]["pooled_rows"]),
        [
            "experiment_run",
            "timepoint",
            "orientation",
            "biological_replicate",
            "loy_count",
            "primary_singlets",
            "loy_fraction",
            "technical_difference_pp",
            "saturated_fraction",
        ],
    )
    _write_csv(
        destination / "source_manifest.csv",
        list(result["new_acquisition"]["source_manifest"]),
        ["path", "size_bytes", "sha256"],
    )
    _write_longitudinal_svg(result, figures / "longitudinal_loy_fraction.svg")


def analyze_repository(repo_root: Path) -> dict:
    """Run the complete prior-plus-new longitudinal analysis."""
    old_result = analyze_prior_repository(repo_root)
    new_acquisition = analyze_new_acquisition(repo_root)
    result = build_longitudinal_result(old_result, new_acquisition["pooled_rows"])
    result["analysis_date"] = "2026-09-04"
    result["classifier"] = old_result["classifier"]
    result["new_acquisition"] = new_acquisition
    return result


def main(argv: list[str] | None = None) -> int:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Generate the D0/D4/D7 longitudinal report")
    parser.add_argument(
        "--output",
        type=Path,
        default=repo / "analysis" / "spectral_competition_longitudinal_20260904",
    )
    args = parser.parse_args(argv)
    result = analyze_repository(repo)
    write_longitudinal_outputs(result, args.output.resolve())
    print(args.output.resolve() / "report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
