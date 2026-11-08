r"""Reproducible D0-to-D4 report for the 2026-09-01 A8 competition assay.

The raw FCS and CEF files are opened read-only.  The report deliberately uses
the run-specific FACSChorus Scatter -> SSC Singlets -> FSC Singlets polygons,
then a two-dimensional linear discriminant learned from the four current-day
single-colour controls.

Run from the repository root::

    .\codex\flow_env\python.exe -m spectral.a8_pipeline.competition_report
"""

from __future__ import annotations

import argparse
import json
import math
import re
import zipfile
from dataclasses import dataclass
from html import escape
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
from scipy.stats import t

from .fcs_io import load_fcs
from .gating import points_in_polygon


@dataclass(frozen=True)
class GateDefinition:
    name: str
    parameters: tuple[str, str]
    vertices: np.ndarray
    parent_population: str
    gate_id: str


@dataclass(frozen=True)
class LDAModel:
    weights: np.ndarray
    intercept: float

    def predict(self, features: np.ndarray) -> np.ndarray:
        values = np.asarray(features, dtype=float)
        if values.ndim != 2 or values.shape[1] != len(self.weights):
            raise ValueError("features have the wrong shape for this LDA model")
        return (values @ self.weights + self.intercept > 0).astype(int)


def background_gate_thresholds(unstained_features: np.ndarray, percentile: float = 99.0) -> np.ndarray:
    """Return per-channel transformed-intensity thresholds from an unstained cloud."""
    values = np.asarray(unstained_features, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2 or len(values) == 0:
        raise ValueError("unstained features must be a non-empty N x 2 array")
    if not 0.0 < float(percentile) < 100.0:
        raise ValueError("percentile must lie between zero and one hundred")
    return np.quantile(values, float(percentile) / 100.0, axis=0)


def apply_background_gate(features: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
    """Classify events: 0 background, 1 GFP+, 2 mCherry+, 3 double-positive."""
    values = np.asarray(features, dtype=float)
    limits = np.asarray(thresholds, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2 or limits.shape != (2,):
        raise ValueError("features must be N x 2 and thresholds must have two values")
    above_gfp = values[:, 0] > limits[0]
    above_mcherry = values[:, 1] > limits[1]
    return above_gfp.astype(int) + 2 * above_mcherry.astype(int)


def load_cef_gates(path: Path) -> dict[str, GateDefinition]:
    """Read polygon gates from the experiment JSON inside a CEF archive."""
    with zipfile.ZipFile(Path(path)) as archive:
        payload = json.loads(archive.read("experiment.json"))
    definitions: dict[str, GateDefinition] = {}
    gates = payload["experiment"]["panels"][0]["analysis"]["gates"]
    for gate in gates:
        if gate.get("gateKind") != "Polygon":
            continue
        parameters = tuple(item["name"] for item in gate.get("parameters", []))
        if len(parameters) != 2:
            raise ValueError(f"polygon gate {gate.get('name')!r} does not have two parameters")
        vertices = np.asarray(
            [[item["x"], item["y"]] for item in gate.get("vertices", [])], dtype=float
        )
        if vertices.shape[0] < 3:
            raise ValueError(f"polygon gate {gate.get('name')!r} has fewer than three vertices")
        definition = GateDefinition(
            name=str(gate["name"]),
            parameters=(str(parameters[0]), str(parameters[1])),
            vertices=vertices,
            parent_population=str(gate.get("parentPopulationId", "")),
            gate_id=str(gate.get("gateId", "")),
        )
        definitions[definition.name] = definition
    return definitions


def fit_lda(features: np.ndarray, labels: np.ndarray) -> LDAModel:
    """Fit equal-prior two-class LDA with a pooled covariance estimate."""
    values = np.asarray(features, dtype=float)
    classes = np.asarray(labels, dtype=int)
    if values.ndim != 2 or len(values) != len(classes):
        raise ValueError("features and labels must have matching rows")
    if set(np.unique(classes)) != {0, 1}:
        raise ValueError("labels must contain both classes 0 and 1")
    groups = [values[classes == value] for value in (0, 1)]
    if any(len(group) < 2 for group in groups):
        raise ValueError("each LDA class needs at least two observations")
    means = [group.mean(axis=0) for group in groups]
    scatter = sum((group - mean).T @ (group - mean) for group, mean in zip(groups, means))
    covariance = scatter / (len(values) - 2)
    weights = np.linalg.pinv(covariance) @ (means[1] - means[0])
    intercept = -0.5 * float((means[1] + means[0]) @ weights)
    return LDAModel(weights=np.asarray(weights, dtype=float), intercept=intercept)


def _sample_features(features: np.ndarray, maximum: int) -> np.ndarray:
    values = np.asarray(features, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("reporter features must have two columns")
    if maximum <= 0:
        raise ValueError("maximum sample size must be positive")
    if len(values) <= maximum:
        return values
    indices = np.linspace(0, len(values) - 1, maximum, dtype=int)
    return values[indices]


def build_classifier_plot_payload(
    classifier: LDAModel,
    control_features: Mapping[str, np.ndarray],
    d0_features: Mapping[str, np.ndarray],
    d4_features: Mapping[str, np.ndarray],
    maximum_points_per_panel: int = 800,
) -> dict:
    """Build compact, deterministic points for classifier diagnostic plots."""
    control_specs = [
        ("A03-CTRL_GFP", "CTRL-GFP", 1),
        ("A05-LOY_GFP", "LOY-GFP", 1),
        ("A07-CTRL_cherry", "CTRL-mCherry", 0),
        ("A09-LOY_cherry", "LOY-mCherry", 0),
    ]
    per_control = max(1, maximum_points_per_panel // len(control_specs))
    control_series = []
    for sample, label, reporter_class in control_specs:
        points = _sample_features(control_features[sample], per_control)
        control_series.append(
            {
                "label": label,
                "reporter_class": reporter_class,
                "points": np.round(points, 4).tolist(),
            }
        )

    labels = {
        "LOY_GFP": "LOY-GFP / CTRL-mCherry",
        "LOY_mCherry": "LOY-mCherry / CTRL-GFP",
    }
    panels = [{"title": "Training controls", "series": control_series}]
    for orientation in ("LOY_GFP", "LOY_mCherry"):
        for day, source in (("D0", d0_features), ("D4", d4_features)):
            sampled = _sample_features(source[orientation], maximum_points_per_panel)
            predictions = classifier.predict(sampled)
            panels.append(
                {
                    "title": f"{day} | {labels[orientation]}",
                    "series": [
                        {
                            "label": "GFP-like",
                            "reporter_class": 1,
                            "points": np.round(sampled[predictions == 1], 4).tolist(),
                        },
                        {
                            "label": "mCherry-like",
                            "reporter_class": 0,
                            "points": np.round(sampled[predictions == 0], 4).tolist(),
                        },
                    ],
                }
            )
    return {
        "weights": np.asarray(classifier.weights, dtype=float).tolist(),
        "intercept": float(classifier.intercept),
        "maximum_points_per_panel": int(maximum_points_per_panel),
        "panels": panels,
    }


def build_all_sample_plot_payload(
    classifier: LDAModel,
    samples: Mapping[str, Mapping],
    maximum_points_per_sample: int = 300,
) -> dict:
    """Build per-sample reporter plots while retaining full-event fractions."""
    group_names = ("controls", "ratio_ladder", "competition_endpoints")
    groups: dict[str, list[dict]] = {name: [] for name in group_names}
    displayed: list[np.ndarray] = []
    for sample_name in sorted(samples):
        sample = samples[sample_name]
        category = str(sample["category"])
        if category not in groups:
            raise ValueError(f"unknown all-sample plot category: {category}")
        features = np.asarray(sample["features"], dtype=float)
        is_unstained = bool(sample.get("unstained", False)) or "unstained" in sample_name.lower()
        total = int(len(features))
        plotted = _sample_features(features, maximum_points_per_sample)
        if is_unstained:
            groups[category].append(
                {
                    "sample": sample_name,
                    "primary_singlets": total,
                    "gfp_count": 0,
                    "gfp_fraction": 0.0,
                    "mcherry_fraction": 0.0,
                    "unstained": True,
                    "series": [
                        {
                            "label": "Unstained background",
                            "reporter_class": -1,
                            "points": np.round(plotted, 4).tolist(),
                        }
                    ],
                }
            )
            displayed.append(plotted)
            continue
        full_predictions = classifier.predict(features)
        plot_predictions = classifier.predict(plotted)
        gfp_count = int(full_predictions.sum())
        groups[category].append(
            {
                "sample": sample_name,
                "primary_singlets": total,
                "gfp_count": gfp_count,
                "gfp_fraction": gfp_count / total,
                "mcherry_fraction": 1.0 - gfp_count / total,
                "unstained": False,
                "series": [
                    {
                        "label": "GFP-like",
                        "reporter_class": 1,
                        "points": np.round(plotted[plot_predictions == 1], 4).tolist(),
                    },
                    {
                        "label": "mCherry-like",
                        "reporter_class": 0,
                        "points": np.round(plotted[plot_predictions == 0], 4).tolist(),
                    },
                ],
            }
        )
        displayed.append(plotted)
    if not displayed:
        x_limits, y_limits = (0.0, 1.0), (0.0, 1.0)
    else:
        values = np.vstack(displayed)
        xmin, ymin = values.min(axis=0)
        xmax, ymax = values.max(axis=0)
        xpad = max(0.25, 0.04 * (xmax - xmin))
        ypad = max(0.25, 0.04 * (ymax - ymin))
        x_limits = (float(xmin - xpad), float(xmax + xpad))
        y_limits = (float(ymin - ypad), float(ymax + ypad))
    return {
        "weights": np.asarray(classifier.weights, dtype=float).tolist(),
        "intercept": float(classifier.intercept),
        "maximum_points_per_sample": int(maximum_points_per_sample),
        "x_limits": list(x_limits),
        "y_limits": list(y_limits),
        "groups": groups,
    }


def expected_gfp_fraction(loy_fraction: float, orientation: str) -> float:
    """Convert a LOY:CTRL mixture fraction to its expected GFP fraction."""
    if orientation == "LOY_GFP":
        return float(loy_fraction)
    if orientation == "CTRL_GFP":
        return float(1.0 - loy_fraction)
    raise ValueError(f"unknown dye orientation: {orientation}")


def pool_technical_counts(counts: Iterable[tuple[int, int]]) -> tuple[int, int, float]:
    """Pool technical aliquots by summing event counts."""
    rows = list(counts)
    positive = int(sum(row[0] for row in rows))
    total = int(sum(row[1] for row in rows))
    if total <= 0:
        raise ValueError("technical replicate total must be positive")
    if positive < 0 or positive > total:
        raise ValueError("positive count must lie between zero and total")
    return positive, total, positive / total


def delta_log_odds(d0_fraction: float, d4_fraction: float) -> float:
    """Return logit(D4) - logit(D0)."""
    if not (0 < d0_fraction < 1 and 0 < d4_fraction < 1):
        raise ValueError("fractions must lie strictly between zero and one")
    return math.log(d4_fraction / (1.0 - d4_fraction)) - math.log(
        d0_fraction / (1.0 - d0_fraction)
    )


def summarise_effects(deltas_by_orientation: Mapping[str, Iterable[float]]) -> dict[str, float]:
    """Estimate the balanced common effect with dye orientation as a fixed factor."""
    groups = [np.asarray(list(values), dtype=float) for values in deltas_by_orientation.values()]
    if len(groups) != 2 or any(len(group) < 2 for group in groups):
        raise ValueError("two dye orientations with at least two plates each are required")
    common = float(np.mean([group.mean() for group in groups]))
    residual_df = int(sum(len(group) - 1 for group in groups))
    sse = float(sum(np.square(group - group.mean()).sum() for group in groups))
    mse = sse / residual_df
    variance = mse * 0.25 * sum(1.0 / len(group) for group in groups)
    half_width = float(t.ppf(0.975, residual_df) * math.sqrt(variance))
    return {
        "common_delta": common,
        "ci_low": common - half_width,
        "ci_high": common + half_width,
        "interval_odds_ratio": math.exp(common),
        "odds_ratio_ci_low": math.exp(common - half_width),
        "odds_ratio_ci_high": math.exp(common + half_width),
        "per_day_log_ratio": common / 4.0,
        "residual_df": float(residual_df),
    }


def _primary_features(path: Path, gates: Mapping[str, GateDefinition]) -> tuple[np.ndarray, int, float]:
    required = {
        "eGFP-A",
        "mCherry-A",
        "Saturated",
        "LightLoss (Violet)-A",
        "SSC (Imaging)-A",
        "SSC (Imaging)-H",
        "SSC (Imaging)-W",
        "LightLoss (Violet)-H",
        "LightLoss (Violet)-W",
    }
    data = load_fcs(path, required=sorted(required))
    events = data.events
    keep = np.ones(len(events), dtype=bool)
    for name in ("Scatter", "SSC Singlets", "FSC Singlets"):
        gate = gates[name]
        inside = points_in_polygon(
            events.loc[:, list(gate.parameters)].to_numpy(dtype=float), gate.vertices
        )
        keep &= inside
    primary = events.loc[keep]
    values = primary.loc[:, ["eGFP-A", "mCherry-A"]].to_numpy(dtype=float)
    transformed = np.arcsinh(values / 1000.0)
    saturated = float((primary["Saturated"].to_numpy(dtype=float) > 0).mean())
    return transformed, len(primary), saturated


def _file_map(root: Path) -> dict[str, Path]:
    return {path.stem: path for path in Path(root).rglob("*.fcs")}


def analyze_repository(repo_root: Path) -> dict:
    """Run the D0/D4 analysis and return JSON-serialisable results."""
    repo = Path(repo_root).resolve()
    d0_root = repo / "spectral" / "20260828_113753_Export" / "Mike_QC1"
    d4_root = repo / "spectral" / "20260901_125230_Mike_Run1" / "Mike_Run1"
    d0_cef = d0_root / "20260828_113753_Mike_QC1.cef"
    d4_cef = d4_root / "20260901_125230_Mike_Run1.cef"
    d0_gates, d4_gates = load_cef_gates(d0_cef), load_cef_gates(d4_cef)
    d0_files, d4_files = _file_map(d0_root), _file_map(d4_root / "1")

    control_roles = {
        "A03-CTRL_GFP": 1,
        "A05-LOY_GFP": 1,
        "A07-CTRL_cherry": 0,
        "A09-LOY_cherry": 0,
    }
    training_features: list[np.ndarray] = []
    training_labels: list[np.ndarray] = []
    control_cache: dict[str, tuple[np.ndarray, int, float]] = {}
    for name, label in control_roles.items():
        sample = _primary_features(d4_files[name], d4_gates)
        control_cache[name] = sample
        training_features.append(sample[0])
        training_labels.append(np.full(sample[1], label, dtype=int))
    classifier = fit_lda(np.vstack(training_features), np.concatenate(training_labels))
    d4_sample_cache = dict(control_cache)

    def d4_sample(name: str) -> tuple[np.ndarray, int, float]:
        if name not in d4_sample_cache:
            d4_sample_cache[name] = _primary_features(d4_files[name], d4_gates)
        return d4_sample_cache[name]

    control_results = []
    for name, expected in control_roles.items():
        features, count, saturated = control_cache[name]
        observed = float(classifier.predict(features).mean())
        control_results.append(
            {
                "sample": name,
                "expected_gfp_class": expected,
                "observed_gfp_fraction": observed,
                "classification_error": abs(observed - expected),
                "primary_singlets": count,
                "saturated_fraction": saturated,
            }
        )

    ratio_rows = []
    ratio_pattern = re.compile(
        r"^[BC]\d+-(10|25|50|75|90)-(10|25|50|75|90)-(CTRL|LOY)_GFP[12](?:_001)?$"
    )
    for name, path in sorted(d4_files.items()):
        match = ratio_pattern.match(name)
        if not match:
            continue
        nominal_loy = float(match.group(1)) / 100.0
        orientation = f"{match.group(3)}_GFP"
        features, count, saturated = d4_sample(name)
        gfp_count = int(classifier.predict(features).sum())
        observed_gfp = gfp_count / count
        observed_loy = observed_gfp if orientation == "LOY_GFP" else 1.0 - observed_gfp
        ratio_rows.append(
            {
                "sample": name,
                "orientation": orientation,
                "nominal_loy": nominal_loy,
                "expected_gfp": expected_gfp_fraction(nominal_loy, orientation),
                "observed_gfp": observed_gfp,
                "observed_loy": observed_loy,
                "primary_singlets": count,
                "saturated_fraction": saturated,
            }
        )

    d0_specs = {
        "LOY_GFP": ("F11-LOY_GFP___CTRL_mCherry", True),
        "LOY_mCherry": ("G11-CTRL_GFP__LOY_mCherry", False),
    }
    orientations: dict[str, dict] = {}
    d0_plot_features: dict[str, np.ndarray] = {}
    for orientation, (name, loy_is_gfp) in d0_specs.items():
        features, count, saturated = _primary_features(d0_files[name], d0_gates)
        d0_plot_features[orientation] = features
        predictions = classifier.predict(features)
        loy_count = int(predictions.sum()) if loy_is_gfp else int((1 - predictions).sum())
        orientations[orientation] = {
            "d0_sample": name,
            "d0_loy_count": loy_count,
            "d0_singlets": count,
            "d0_loy": loy_count / count,
            "d0_saturated_fraction": saturated,
            "plates": [],
        }

    endpoint_rows = []
    d4_plot_parts: dict[str, list[np.ndarray]] = {"LOY_GFP": [], "LOY_mCherry": []}
    for name, path in sorted(d4_files.items()):
        if "_exp_rep" not in name:
            continue
        orientation = "LOY_GFP" if name.startswith("E") else "LOY_mCherry"
        loy_is_gfp = orientation == "LOY_GFP"
        replicate = int(re.search(r"rep(\d)", name).group(1))
        features, count, saturated = d4_sample(name)
        d4_plot_parts[orientation].append(features)
        predictions = classifier.predict(features)
        loy_count = int(predictions.sum()) if loy_is_gfp else int((1 - predictions).sum())
        endpoint_rows.append(
            {
                "sample": name,
                "orientation": orientation,
                "replicate": replicate,
                "loy_count": loy_count,
                "primary_singlets": count,
                "loy_fraction": loy_count / count,
                "saturated_fraction": saturated,
            }
        )

    maximum_technical_difference = 0.0
    deltas_by_orientation: dict[str, list[float]] = {name: [] for name in orientations}
    for orientation, data in orientations.items():
        for replicate in (1, 2, 3):
            wells = [
                row
                for row in endpoint_rows
                if row["orientation"] == orientation and row["replicate"] == replicate
            ]
            if len(wells) != 2:
                raise ValueError(f"expected two technical wells for {orientation} replicate {replicate}")
            positive, total, fraction = pool_technical_counts(
                (int(row["loy_count"]), int(row["primary_singlets"])) for row in wells
            )
            technical_difference = abs(wells[0]["loy_fraction"] - wells[1]["loy_fraction"])
            maximum_technical_difference = max(maximum_technical_difference, technical_difference)
            delta = delta_log_odds(float(data["d0_loy"]), fraction)
            deltas_by_orientation[orientation].append(delta)
            data["plates"].append(
                {
                    "replicate": replicate,
                    "loy_count": positive,
                    "primary_singlets": total,
                    "d4_loy": fraction,
                    "technical_difference_pp": 100.0 * technical_difference,
                    "delta_log_odds": delta,
                    "interval_odds_ratio": math.exp(delta),
                    "technical_wells": wells,
                }
            )
        data["mean_delta_log_odds"] = float(np.mean(deltas_by_orientation[orientation]))
        data["sd_delta_log_odds"] = float(np.std(deltas_by_orientation[orientation], ddof=1))
        data["mean_interval_odds_ratio"] = math.exp(data["mean_delta_log_odds"])

    summary = summarise_effects(deltas_by_orientation)
    d4_plot_features = {
        orientation: np.vstack(parts) for orientation, parts in d4_plot_parts.items()
    }
    classifier_plot = build_classifier_plot_payload(
        classifier,
        {name: sample[0] for name, sample in control_cache.items()},
        d0_plot_features,
        d4_plot_features,
    )
    all_sample_inputs: dict[str, dict] = {}
    for name in sorted(d4_files):
        features, _count, _saturated = d4_sample(name)
        if name.startswith("A"):
            category = "controls"
        elif name.startswith(("B", "C")):
            category = "ratio_ladder"
        elif name.startswith(("D", "E")):
            category = "competition_endpoints"
        else:
            raise ValueError(f"unrecognised September sample group: {name}")
        all_sample_inputs[name] = {"category": category, "features": features, "unstained": "unstained" in name.lower()}
    all_sample_plots = build_all_sample_plot_payload(classifier, all_sample_inputs)
    return {
        "title": "Spectral competition assay: D0 to D4",
        "analysis_date": "2026-09-01",
        "source_runs": {
            "d0": str(d0_root.relative_to(repo)),
            "d4": str(d4_root.relative_to(repo)),
        },
        "classifier": {
            "name": "two-dimensional equal-prior LDA on arcsinh(eGFP-A/1000, mCherry-A/1000)",
            "controls": control_results,
            "maximum_control_error": max(row["classification_error"] for row in control_results),
        },
        "classifier_plot": classifier_plot,
        "all_sample_plots": all_sample_plots,
        "ratio_ladder": ratio_rows,
        "orientations": orientations,
        "summary": summary,
        "technical_max_difference_pp": 100.0 * maximum_technical_difference,
        "primary_gate_order": ["Scatter", "SSC Singlets", "FSC Singlets"],
        "unmixing_matrix": "mcherry_GFP_matrix",
    }


def _svg_header(width: int, height: int, title: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img">',
        f"<title>{escape(title)}</title>",
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2:.0f}" y="30" text-anchor="middle" font-family="sans-serif" font-size="18">{escape(title)}</text>',
    ]


def _write_svg(parts: list[str], output: Path) -> None:
    parts.append("</svg>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(parts), encoding="utf-8")


def write_classifier_workflow_svg(output: Path) -> None:
    width, height = 1040, 290
    parts = _svg_header(width, height, "How reporter identity was assigned")
    parts.extend(
        [
            '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#555"/></marker></defs>',
        ]
    )
    stages = [
        ("1", "Primary singlets", ["Scatter", "SSC Singlets", "FSC Singlets"]),
        ("2", "Reporter transform", ["eGFP-A and mCherry-A", "arcsinh(value / 1000)"]),
        ("3", "Learn boundary", ["Four D4 single-colour controls", "Equal-prior 2D LDA"]),
        (
            "4",
            "Assign every event",
            ["Same boundary at D0 and D4", "score > 0: GFP-like", "score ≤ 0: mCherry-like"],
        ),
    ]
    stages[-1] = (stages[-1][0], stages[-1][1], list(stages[-1][2]) + ["Unstained controls: 0%/0% background"])
    box_width, box_height, gap, x0, y = 220, 156, 36, 26, 64
    for index, (number, heading, details) in enumerate(stages):
        x = x0 + index * (box_width + gap)
        parts.append(
            f'<g data-stage="{number}"><rect x="{x}" y="{y}" width="{box_width}" height="{box_height}" rx="10" fill="#f7f7f7" stroke="#777"/>'
            f'<circle cx="{x + 23}" cy="{y + 25}" r="14" fill="#333"/>'
            f'<text x="{x + 23}" y="{y + 30}" text-anchor="middle" font-family="sans-serif" font-size="13" fill="white">{number}</text>'
            f'<text x="{x + 45}" y="{y + 30}" font-family="sans-serif" font-size="14" font-weight="bold">{escape(heading)}</text>'
        )
        for line_index, line in enumerate(details):
            parts.append(
                f'<text x="{x + 16}" y="{y + 66 + line_index * 24}" font-family="sans-serif" font-size="12">{escape(line)}</text>'
            )
        parts.append("</g>")
        if index < len(stages) - 1:
            x1 = x + box_width + 5
            x2 = x + box_width + gap - 6
            parts.append(
                f'<line data-connector="{number}" x1="{x1}" y1="{y + box_height / 2}" x2="{x2}" y2="{y + box_height / 2}" stroke="#555" stroke-width="2" marker-end="url(#arrow)"/>'
            )
    _write_svg(parts, output)


def _boundary_segment(
    weights: Iterable[float],
    intercept: float,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
) -> tuple[tuple[float, float], tuple[float, float]] | None:
    wx, wy = (float(value) for value in weights)
    xmin, xmax = x_limits
    ymin, ymax = y_limits
    candidates: list[tuple[float, float]] = []
    if abs(wy) > 1e-12:
        for x in (xmin, xmax):
            y = -(wx * x + float(intercept)) / wy
            if ymin <= y <= ymax:
                candidates.append((x, y))
    if abs(wx) > 1e-12:
        for y in (ymin, ymax):
            x = -(wy * y + float(intercept)) / wx
            if xmin <= x <= xmax:
                candidates.append((x, y))
    unique: list[tuple[float, float]] = []
    for point in candidates:
        if not any(np.allclose(point, other) for other in unique):
            unique.append(point)
    return (unique[0], unique[1]) if len(unique) >= 2 else None


def write_reporter_space_svg(payload: Mapping, output: Path) -> None:
    panels = list(payload.get("panels", []))
    all_points = [
        point
        for panel in panels
        for series in panel.get("series", [])
        for point in series.get("points", [])
    ]
    if all_points:
        values = np.asarray(all_points, dtype=float)
        xmin, ymin = values.min(axis=0)
        xmax, ymax = values.max(axis=0)
        xpad = max(0.25, 0.04 * (xmax - xmin))
        ypad = max(0.25, 0.04 * (ymax - ymin))
        x_limits = (float(xmin - xpad), float(xmax + xpad))
        y_limits = (float(ymin - ypad), float(ymax + ypad))
    else:
        x_limits, y_limits = (0.0, 1.0), (0.0, 1.0)

    width, height = 1040, 1040
    plot_width, plot_height = 390, 210
    positions = [(325, 110), (85, 430), (565, 430), (85, 750), (565, 750)]
    parts = _svg_header(width, height, "Reporter-space classification with one shared boundary")
    parts.extend(
        [
            '<g aria-label="legend">',
            '<circle cx="245" cy="50" r="5" fill="#247a3b"/><text x="257" y="54" font-family="sans-serif" font-size="12">GFP-like</text>',
            '<circle cx="350" cy="50" r="5" fill="#b33a63"/><text x="362" y="54" font-family="sans-serif" font-size="12">mCherry-like</text>',
            '<line x1="475" y1="50" x2="510" y2="50" stroke="#222" stroke-width="2" stroke-dasharray="6,4"/><text x="519" y="54" font-family="sans-serif" font-size="12">LDA boundary</text>',
            '<circle cx="650" cy="50" r="5" fill="#777"/><text x="662" y="54" font-family="sans-serif" font-size="12">CTRL control</text>',
            '<polygon points="792,44 786,55 798,55" fill="#777"/><text x="806" y="54" font-family="sans-serif" font-size="12">LOY control</text>',
            '</g>',
        ]
    )

    def px(value: float, left: float) -> float:
        return left + (value - x_limits[0]) / (x_limits[1] - x_limits[0]) * plot_width

    def py(value: float, top: float) -> float:
        return top + plot_height - (value - y_limits[0]) / (y_limits[1] - y_limits[0]) * plot_height

    boundary = _boundary_segment(
        payload.get("weights", [1.0, -1.0]),
        float(payload.get("intercept", 0.0)),
        x_limits,
        y_limits,
    )
    for panel_index, panel in enumerate(panels):
        left, top = positions[panel_index]
        parts.append(f'<g data-panel="{panel_index + 1}">')
        parts.append(
            f'<text x="{left + plot_width / 2}" y="{top - 14}" text-anchor="middle" font-family="sans-serif" font-size="14">{escape(str(panel["title"]))}</text>'
        )
        parts.append(
            f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="#fbfbfb" stroke="#777"/>'
        )
        for tick in np.linspace(x_limits[0], x_limits[1], 4):
            x = px(float(tick), left)
            parts.append(
                f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_height}" stroke="#e5e5e5"/>'
            )
            parts.append(
                f'<text data-axis-tick="x" x="{x:.2f}" y="{top + plot_height + 17}" text-anchor="middle" font-family="sans-serif" font-size="11">{tick:.1f}</text>'
            )
        for tick in np.linspace(y_limits[0], y_limits[1], 4):
            y = py(float(tick), top)
            parts.append(
                f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" stroke="#e5e5e5"/>'
            )
            parts.append(
                f'<text data-axis-tick="y" x="{left - 8}" y="{y + 4:.2f}" text-anchor="end" font-family="sans-serif" font-size="11">{tick:.1f}</text>'
            )
        for series in panel.get("series", []):
            reporter_class = int(series.get("reporter_class", 0))
            color = "#247a3b" if reporter_class == 1 else "#b33a63"
            label = str(series.get("label", ""))
            is_loy = label.startswith("LOY-")
            for x_value, y_value in series.get("points", []):
                x, y = px(float(x_value), left), py(float(y_value), top)
                common = f'data-reporter-class="{reporter_class}" data-series="{escape(label)}" fill="{color}" fill-opacity="0.28"'
                if is_loy:
                    points = f"{x:.2f},{y-2.8:.2f} {x-2.8:.2f},{y+2.3:.2f} {x+2.8:.2f},{y+2.3:.2f}"
                    parts.append(f'<polygon {common} points="{points}"/>')
                else:
                    parts.append(f'<circle {common} cx="{x:.2f}" cy="{y:.2f}" r="2.1"/>')
        if boundary is not None:
            (x1, y1), (x2, y2) = boundary
            parts.append(
                f'<line data-boundary="shared" x1="{px(x1, left):.2f}" y1="{py(y1, top):.2f}" x2="{px(x2, left):.2f}" y2="{py(y2, top):.2f}" stroke="#222" stroke-width="2" stroke-dasharray="6,4"/>'
            )
        parts.extend(
            [
                f'<text x="{left + plot_width / 2}" y="{top + plot_height + 38}" text-anchor="middle" font-family="sans-serif" font-size="12">arcsinh(eGFP-A / 1000)</text>',
                f'<text x="{left - 48}" y="{top + plot_height / 2}" transform="rotate(-90 {left - 48} {top + plot_height / 2})" text-anchor="middle" font-family="sans-serif" font-size="12">arcsinh(mCherry-A / 1000)</text>',
                "</g>",
            ]
        )
    _write_svg(parts, output)


def write_all_sample_sheet_svg(payload: Mapping, group_name: str, output: Path) -> None:
    specifications = {
        "controls": ("September 1 controls", 3),
        "ratio_ladder": ("September 1 known-ratio samples", 4),
        "competition_endpoints": ("September 1 competition endpoint wells", 3),
    }
    if group_name not in specifications:
        raise ValueError(f"unknown all-sample sheet: {group_name}")
    title, columns = specifications[group_name]
    panels = list(payload.get("groups", {}).get(group_name, []))
    rows = max(1, math.ceil(len(panels) / columns))
    cell_width, cell_height = 340, 250
    grid_left, grid_top = 20, 78
    width = grid_left * 2 + columns * cell_width
    height = grid_top + rows * cell_height + 24
    plot_width, plot_height = 250, 155
    x_limits = tuple(float(value) for value in payload.get("x_limits", [0.0, 1.0]))
    y_limits = tuple(float(value) for value in payload.get("y_limits", [0.0, 1.0]))
    boundary = _boundary_segment(
        payload.get("weights", [1.0, -1.0]),
        float(payload.get("intercept", 0.0)),
        x_limits,
        y_limits,
    )
    parts = _svg_header(width, height, title)
    parts.extend(
        [
            f'<desc>Each panel uses identical transformed reporter axes and the same LDA boundary. Percentages use all primary singlets; displayed points are a deterministic subsample.</desc>',
            '<g aria-label="legend">',
            f'<circle cx="{width / 2 - 155:.1f}" cy="52" r="5" fill="#247a3b"/><text x="{width / 2 - 143:.1f}" y="56" font-family="sans-serif" font-size="12">GFP-like</text>',
            f'<circle cx="{width / 2 - 55:.1f}" cy="52" r="5" fill="#b33a63"/><text x="{width / 2 - 43:.1f}" y="56" font-family="sans-serif" font-size="12">mCherry-like</text>',
            f'<line x1="{width / 2 + 72:.1f}" y1="52" x2="{width / 2 + 107:.1f}" y2="52" stroke="#222" stroke-width="2" stroke-dasharray="6,4"/><text x="{width / 2 + 116:.1f}" y="56" font-family="sans-serif" font-size="12">shared LDA boundary</text>',
            f'<circle cx="{width - 195:.1f}" cy="52" r="5" fill="#777777"/><text x="{width - 183:.1f}" y="56" font-family="sans-serif" font-size="12">unstained background</text>',
            "</g>",
        ]
    )

    def px(value: float, left: float) -> float:
        return left + (value - x_limits[0]) / (x_limits[1] - x_limits[0]) * plot_width

    def py(value: float, top: float) -> float:
        return top + plot_height - (value - y_limits[0]) / (y_limits[1] - y_limits[0]) * plot_height

    for index, panel in enumerate(panels):
        column, row = index % columns, index // columns
        cell_x = grid_left + column * cell_width
        cell_y = grid_top + row * cell_height
        left, top = cell_x + 62, cell_y + 42
        sample_name = str(panel["sample"])
        parts.append(
            f'<g data-sample="{escape(sample_name)}" data-x-min="{x_limits[0]:g}" data-x-max="{x_limits[1]:g}" data-y-min="{y_limits[0]:g}" data-y-max="{y_limits[1]:g}">'
        )
        parts.append(
            f'<text x="{cell_x + cell_width / 2:.1f}" y="{cell_y + 13}" text-anchor="middle" font-family="sans-serif" font-size="11">{escape(sample_name)}</text>'
        )
        parts.append(
            f'<text x="{cell_x + cell_width / 2:.1f}" y="{cell_y + 29}" text-anchor="middle" font-family="sans-serif" font-size="11">{("n=" + format(int(panel["primary_singlets"]), ",") + " | GFP+ 0.0% | mCherry+ 0.0% | background") if panel.get("unstained") else ("n=" + format(int(panel["primary_singlets"]), ",") + " | GFP-like " + format(100.0 * float(panel["gfp_fraction"]), ".1f") + "% | mCherry-like " + format(100.0 * float(panel["mcherry_fraction"]), ".1f") + "%")}</text>'
        )
        parts.append(
            f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="#fbfbfb" stroke="#777"/>'
        )
        for tick in np.linspace(x_limits[0], x_limits[1], 3):
            x = px(float(tick), left)
            parts.append(
                f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_height}" stroke="#e5e5e5"/>'
            )
            parts.append(
                f'<text data-axis-tick="x" x="{x:.2f}" y="{top + plot_height + 15}" text-anchor="middle" font-family="sans-serif" font-size="11">{tick:.1f}</text>'
            )
        for tick in np.linspace(y_limits[0], y_limits[1], 3):
            y = py(float(tick), top)
            parts.append(
                f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" stroke="#e5e5e5"/>'
            )
            parts.append(
                f'<text data-axis-tick="y" x="{left - 7}" y="{y + 3:.2f}" text-anchor="end" font-family="sans-serif" font-size="11">{tick:.1f}</text>'
            )
        for series in panel.get("series", []):
            reporter_class = int(series["reporter_class"])
            color = "#777777" if reporter_class == -1 else ("#247a3b" if reporter_class == 1 else "#b33a63")
            for x_value, y_value in series.get("points", []):
                parts.append(
                    f'<circle data-reporter-class="{reporter_class}" cx="{px(float(x_value), left):.2f}" cy="{py(float(y_value), top):.2f}" r="1.8" fill="{color}" fill-opacity="0.32"/>'
                )
        if boundary is not None:
            (x1, y1), (x2, y2) = boundary
            parts.append(
                f'<line data-boundary="shared" x1="{px(x1, left):.2f}" y1="{py(y1, top):.2f}" x2="{px(x2, left):.2f}" y2="{py(y2, top):.2f}" stroke="#222" stroke-width="1.6" stroke-dasharray="5,4"/>'
            )
        parts.extend(
            [
                f'<text x="{left + plot_width / 2}" y="{top + plot_height + 34}" text-anchor="middle" font-family="sans-serif" font-size="11">arcsinh(eGFP-A / 1000)</text>',
                f'<text x="{left - 43}" y="{top + plot_height / 2}" transform="rotate(-90 {left - 43} {top + plot_height / 2})" text-anchor="middle" font-family="sans-serif" font-size="11">arcsinh(mCherry-A / 1000)</text>',
                "</g>",
            ]
        )
    if not panels:
        parts.append(
            f'<text x="{width / 2:.1f}" y="{grid_top + 80}" text-anchor="middle" font-family="sans-serif" font-size="13">No samples</text>'
        )
    _write_svg(parts, output)


def write_ratio_ladder_svg(rows: Iterable[Mapping], output: Path) -> None:
    data = list(rows)
    width, height = 760, 620
    left, top, span = 85, 60, 480
    parts = _svg_header(width, height, "Known-ratio recovery by dye orientation")

    def px(value: float) -> float:
        return left + value * span

    def py(value: float) -> float:
        return top + span - value * span

    parts.extend(
        [
            f'<rect x="{left}" y="{top}" width="{span}" height="{span}" fill="none" stroke="#333"/>',
            f'<line x1="{px(0):.1f}" y1="{py(0):.1f}" x2="{px(1):.1f}" y2="{py(1):.1f}" stroke="#777" stroke-width="2" stroke-dasharray="7,5"/>',
            f'<text x="{px(.77):.1f}" y="{py(.83):.1f}" font-family="sans-serif" font-size="11" fill="#555">expected</text>',
        ]
    )
    for tick in np.linspace(0, 1, 6):
        parts.append(
            f'<text x="{px(tick):.1f}" y="{top + span + 20}" text-anchor="middle" font-family="sans-serif" font-size="11">{100*tick:.0f}</text>'
        )
        parts.append(
            f'<text x="{left - 12}" y="{py(tick) + 4:.1f}" text-anchor="end" font-family="sans-serif" font-size="11">{100*tick:.0f}</text>'
        )
    styles = {
        "LOY_GFP": ("#247a3b", "circle", "LOY-GFP / CTRL-mCherry"),
        "CTRL_GFP": ("#b33a3a", "triangle", "LOY-mCherry / CTRL-GFP"),
    }
    for row in data:
        color, shape, _label = styles[str(row["orientation"])]
        x, y = px(float(row["nominal_loy"])), py(float(row["observed_loy"]))
        if shape == "circle":
            parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="{color}" fill-opacity="0.78"/>')
        else:
            points = f"{x:.1f},{y-6:.1f} {x-6:.1f},{y+5:.1f} {x+6:.1f},{y+5:.1f}"
            parts.append(f'<polygon points="{points}" fill="{color}" fill-opacity="0.78"/>')
    parts.extend(
        [
            f'<text x="{left + span / 2}" y="{height - 38}" text-anchor="middle" font-family="sans-serif" font-size="13">Nominal LOY (%)</text>',
            f'<text x="24" y="{top + span / 2}" transform="rotate(-90 24 {top + span / 2})" text-anchor="middle" font-family="sans-serif" font-size="13">Observed LOY (%)</text>',
        ]
    )
    # Boxed legend inside the axes: the ladder hugs the 1:1 diagonal, so the
    # upper-left corner is always free and nothing spills past the canvas.
    box_x, box_y, box_w, box_h = 97, 72, 252, 64
    marker_x, label_x, legend_y = 115, 130, 96
    parts.append('<g data-legend="ratio-ladder">')
    parts.append(
        f'<rect x="{box_x}" y="{box_y}" width="{box_w}" height="{box_h}" rx="4" '
        f'fill="#ffffff" fill-opacity="0.92" stroke="#555" stroke-width="1"/>'
    )
    for index, (_key, (color, shape, label)) in enumerate(styles.items()):
        y = legend_y + 26 * index
        if shape == "circle":
            parts.append(f'<circle cx="{marker_x}" cy="{y}" r="5" fill="{color}"/>')
        else:
            parts.append(
                f'<polygon points="{marker_x},{y-6} {marker_x-6},{y+5} {marker_x+6},{y+5}" fill="{color}"/>'
            )
        parts.append(
            f'<text x="{label_x}" y="{y+4}" font-family="sans-serif" font-size="12">{escape(label)}</text>'
        )
    parts.append("</g>")
    _write_svg(parts, output)


def write_competition_svg(result: Mapping, output: Path) -> None:
    width, height = 820, 620
    left, top, plot_width, plot_height = 75, 60, 680, 470
    parts = _svg_header(width, height, "LOY fraction declines from D0 to D4")

    def py(value: float) -> float:
        return top + plot_height - value * plot_height

    x_positions = {"LOY_GFP": (155, 300), "LOY_mCherry": (505, 650)}
    styles = {"LOY_GFP": "#247a3b", "LOY_mCherry": "#b33a3a"}
    for tick in np.linspace(0, 1, 6):
        y = py(float(tick))
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left+plot_width}" y2="{y:.1f}" stroke="#e2e2e2"/>')
        parts.append(f'<text x="{left-10}" y="{y+4:.1f}" text-anchor="end" font-family="sans-serif" font-size="11">{100*tick:.0f}</text>')
    for orientation, data in result["orientations"].items():
        x0, x4 = x_positions[orientation]
        color = styles[orientation]
        d0 = float(data["d0_loy"])
        for plate in data["plates"]:
            d4 = float(plate["d4_loy"])
            parts.append(f'<line x1="{x0}" y1="{py(d0):.1f}" x2="{x4}" y2="{py(d4):.1f}" stroke="{color}" stroke-width="2" stroke-opacity="0.65"/>')
            parts.append(f'<circle cx="{x4}" cy="{py(d4):.1f}" r="6" fill="{color}"/>')
            for index, well in enumerate(plate.get("technical_wells", [])):
                offset = -8 if index == 0 else 8
                parts.append(f'<circle cx="{x4+offset}" cy="{py(float(well["loy_fraction"])):.1f}" r="2.5" fill="white" stroke="{color}"/>')
        parts.append(f'<circle cx="{x0}" cy="{py(d0):.1f}" r="7" fill="white" stroke="{color}" stroke-width="3"/>')
    for orientation, (x0, x4) in x_positions.items():
        parts.append(f'<text x="{x0}" y="{top+plot_height+22}" text-anchor="middle" font-family="sans-serif" font-size="12">Day 0</text>')
        parts.append(f'<text x="{x4}" y="{top+plot_height+22}" text-anchor="middle" font-family="sans-serif" font-size="12">Day 4</text>')
    parts.extend(
        [
            f'<text x="{(155+300)/2}" y="{height-30}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#247a3b">LOY-GFP / CTRL-mCherry</text>',
            f'<text x="{(505+650)/2}" y="{height-30}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#b33a3a">LOY-mCherry / CTRL-GFP</text>',
            f'<text x="22" y="{top+plot_height/2}" transform="rotate(-90 22 {top+plot_height/2})" text-anchor="middle" font-family="sans-serif" font-size="13">LOY fraction (%)</text>',
        ]
    )
    _write_svg(parts, output)


def _percent(value: float) -> str:
    return f"{100.0 * value:.2f}%"


def build_report_markdown(result: Mapping) -> str:
    summary = result["summary"]
    orientations = result.get("orientations", {})
    lines = [
        "# Spectral competition assay: D0 to D4",
        "",
        "**Run dates:** 28 August 2026 (D0) and 1 September 2026 (D4)  ",
        "**Instrument:** BD FACSDiscover A8  ",
        "**Comparison:** CTRL versus LOY HT1080 cells with reciprocal GFP/mCherry labelling",
        "",
        "## Executive summary",
        "",
        f"LOY representation decreased in all six independently cultured plates. The dye-orientation-adjusted four-day change was **{summary['common_delta']:.3f} log-odds** (plate-variance 95% interval {summary['ci_low']:.3f} to {summary['ci_high']:.3f}), corresponding to a LOY:CTRL interval odds ratio of **{summary['interval_odds_ratio']:.3f}** ({summary.get('odds_ratio_ci_low', math.nan):.3f} to {summary.get('odds_ratio_ci_high', math.nan):.3f}). Thus, the LOY-to-CTRL ratio approximately halved over four days under this co-culture condition.",
        "",
        "This is evidence for a net competitive disadvantage of the tested LOY cells under these conditions. It does not distinguish proliferation, death, attachment, or harvest recovery, and it does not imply that the absolute LOY cell number halved.",
        "",
        "## Experimental design",
        "",
        "Each reciprocal dye orientation had one measured D0 mixture. Each D0 mixture was split into three independently cultured plates on 28 August. Cultures remained continuously in the original plates and were harvested near 80% confluence on 1 September. Each plate culture was split into two technical assay wells only at harvest.",
        "",
        "The statistical unit is therefore the plate culture (`n=3` per dye orientation). Technical wells are pooled by event counts and are not treated as independent replicates.",
        "",
        "## Methods",
        "",
        "1. Raw FCS files were read without modification. The run-specific FACSChorus polygon hierarchy—Scatter, SSC Singlets, then FSC Singlets—was extracted from each CEF and applied in raw units.",
        "2. Unmixed `eGFP-A` and `mCherry-A` values were transformed as `arcsinh(value / 1000)`.",
        "3. An equal-prior two-dimensional linear discriminant was trained from the four September single-colour controls: CTRL-GFP, LOY-GFP, CTRL-mCherry, and LOY-mCherry.",
        "4. The same boundary was applied to D0 and D4 because both exports use the identical `mcherry_GFP_matrix` spectral-unmixing hash and optical configuration.",
        "Unstained controls were treated as background controls rather than forced into the two-class reporter classifier; they are reported as 0% GFP-positive / 0% mCherry-positive.",
        "5. Technical aliquots were pooled by summing classified LOY events and total primary singlets.",
        "6. Plate effects were calculated as `logit(D4 LOY fraction) - logit(D0 LOY fraction)`. The common estimate balances the two dye orientations and uses within-orientation plate variance.",
        "",
        "The primary analysis retains instrument-flagged `Saturated` events. Excluding all such events strongly distorted known-ratio recovery; prior mask decoding showed that the flag is dominated by the ImgB3 imaging detector rather than spectral detectors.",
        "",
        "## Plate-level results",
        "",
        "| Dye orientation | Plate | D0 LOY | D4 LOY | Technical difference | Δ log-odds | Four-day odds ratio |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for orientation, data in orientations.items():
        label = "LOY-GFP / CTRL-mCherry" if orientation == "LOY_GFP" else "LOY-mCherry / CTRL-GFP"
        for plate in data["plates"]:
            lines.append(
                f"| {label} | {plate['replicate']} | {_percent(data['d0_loy'])} | {_percent(plate['d4_loy'])} | {plate['technical_difference_pp']:.2f} pp | {plate['delta_log_odds']:.3f} | {plate['interval_odds_ratio']:.3f} |"
            )
    lines.extend(["", "## Dye-orientation summaries", ""])
    for orientation, data in orientations.items():
        label = "LOY-GFP / CTRL-mCherry" if orientation == "LOY_GFP" else "LOY-mCherry / CTRL-GFP"
        mean_d4 = float(np.mean([plate["d4_loy"] for plate in data["plates"]]))
        lines.append(
            f"- **{label}:** D0 {_percent(data['d0_loy'])}; mean D4 {_percent(mean_d4)}; mean Δ log-odds {data['mean_delta_log_odds']:.3f}; interval odds ratio {data['mean_interval_odds_ratio']:.3f}."
        )
    lines.extend(
        [
            "",
            "## Quality control",
            "",
            f"- Maximum technical-well disagreement was **{result.get('technical_max_difference_pp', math.nan):.2f} percentage points**.",
            f"- Maximum pure-control classification error was **{100.0 * result.get('classifier', {}).get('maximum_control_error', math.nan):.2f}%**.",
            "- Three alternative discriminant models were checked during exploratory analysis and produced a common four-day odds ratio of approximately 0.50.",
            "- Both technical aliquots from LOY-GFP plate 3 gave nearly identical LOY proportions despite its lower singlet recovery, supporting the stability of its composition estimate.",
            "- All six plate-level Δ log-odds estimates were negative.",
            "",
            "### Known-ratio ladder",
            "",
            "| Orientation | Nominal LOY | Mean observed LOY |",
            "|---|---:|---:|",
        ]
    )
    ladder = result.get("ratio_ladder", [])
    for orientation in ("LOY_GFP", "CTRL_GFP"):
        label = "LOY-GFP / CTRL-mCherry" if orientation == "LOY_GFP" else "LOY-mCherry / CTRL-GFP"
        for nominal in (0.10, 0.25, 0.50, 0.75, 0.90):
            values = [row["observed_loy"] for row in ladder if row["orientation"] == orientation and abs(row["nominal_loy"] - nominal) < 1e-9]
            if values:
                lines.append(f"| {label} | {_percent(nominal)} | {_percent(float(np.mean(values)))} |")
    lines.extend(
        [
            "",
            "The LOY-GFP ladder tracks the intended mixtures closely. The reciprocal ladder is monotonic but under-recovers LOY at parts of the range, especially the nominal 90% point. Agreement between its duplicate wells suggests a mixture-preparation or cell-counting difference rather than random classification noise. The ladder was therefore used as QC, not as a mathematical correction.",
            "",
            "## Figures",
            "",
            "![Reporter classification workflow](figures/classification_workflow.svg)",
            "",
            "*Reporter classification workflow.* Primary singlets were transformed, used to learn one linear boundary from the four single-colour controls, and then assigned to one of two mutually exclusive reporter classes.",
            "",
            "![Reporter-space classification](figures/reporter_classification.svg)",
            "",
            "*Reporter-space classification.* Pure controls and both dye orientations at D0 and D4 are shown on shared transformed axes. The dashed decision boundary is identical in every panel.",
            "",
            "### All September 1 samples",
            "",
            "Every panel below uses identical transformed axes and the same classifier boundary. Displayed points are proportion-preserving subsamples; the annotated reporter percentages use all primary singlets.",
            "",
            "- [Controls — 6 samples](figures/all_samples/controls.svg)",
            "- [Known-ratio ladder — 20 samples](figures/all_samples/ratio_ladder.svg)",
            "- [Competition endpoint wells — 12 samples](figures/all_samples/competition_endpoints.svg)",
            "",
            "![Known-ratio recovery](figures/ratio_ladder.svg)",
            "",
            "*Known-ratio recovery.* Duplicate wells are shown separately; the dashed line is perfect recovery.",
            "",
            "![D0 to D4 LOY fraction](figures/d0_to_d4.svg)",
            "",
            "*D0 to D4 LOY fraction.* Each coloured line is one independent plate culture; small open points are technical assay wells.",
            "",
            "## Limitations",
            "",
            "- Each dye orientation has one shared D0 mixture. The three plates quantify post-split culture variability but not variability among independently prepared D0 mixtures.",
            "- The interval reflects plate-to-plate variation conditional on the shared D0 values and classification model; it is not a full uncertainty interval for all preparation and gating choices.",
            "- D0 LOY-GFP/CTRL-mCherry was acquired at a different imaging flow-rate setting from D4.",
            "- There was no viability dye, so the assay measures net competitive representation rather than a specific proliferation or survival mechanism.",
            "- Cultures were near 80% confluence at harvest; the result applies to this density and four-day continuous-culture condition.",
            "- Results concern the tested cell preparations and do not by themselves generalize to all LOY states or clones.",
            "",
            "## Recommendation",
            "",
            "Repeat with three independently prepared and measured D0 mixtures per dye orientation. Match acquisition flow mode across days, equalize viable cell concentrations before preparing ratio standards, preserve the reciprocal dye swap, and retain technical harvest aliquots as QC rather than independent replicates.",
            "",
            "## Reproduction",
            "",
            "From the repository root:",
            "",
            "```powershell",
            ".\\codex\\flow_env\\python.exe -m spectral.a8_pipeline.competition_report",
            "```",
            "",
            "The command reads the raw FCS and CEF inputs and regenerates this report, `results.json`, and the SVG figures. It refuses to overwrite an existing output directory unless `--force` is supplied.",
            "",
        ]
    )
    return "\n".join(lines)


def clear_generated_output(output: Path) -> None:
    """Remove generated files while retaining sync-managed directories."""
    destination = Path(output)
    for path in destination.rglob("*"):
        if path.is_file():
            path.unlink()


def write_outputs(result: Mapping, output: Path) -> None:
    destination = Path(output)
    if destination.exists():
        if not destination.is_dir() or any(path.is_file() for path in destination.rglob("*")):
            raise FileExistsError(f"output is not an empty directory: {destination}")
    else:
        destination.mkdir(parents=True)
    figures = destination / "figures"
    write_classifier_workflow_svg(figures / "classification_workflow.svg")
    write_reporter_space_svg(result["classifier_plot"], figures / "reporter_classification.svg")
    all_sample_figures = figures / "all_samples"
    write_all_sample_sheet_svg(
        result["all_sample_plots"], "controls", all_sample_figures / "controls.svg"
    )
    write_all_sample_sheet_svg(
        result["all_sample_plots"], "ratio_ladder", all_sample_figures / "ratio_ladder.svg"
    )
    write_all_sample_sheet_svg(
        result["all_sample_plots"],
        "competition_endpoints",
        all_sample_figures / "competition_endpoints.svg",
    )
    write_ratio_ladder_svg(result["ratio_ladder"], figures / "ratio_ladder.svg")
    write_competition_svg(result, figures / "d0_to_d4.svg")
    (destination / "results.json").write_text(
        json.dumps(result, indent=2, sort_keys=True), encoding="utf-8"
    )
    (destination / "report.md").write_text(build_report_markdown(result), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Generate the 2026-09-01 competition report")
    parser.add_argument(
        "--output",
        type=Path,
        default=repo / "analysis" / "spectral_competition_20260901",
    )
    parser.add_argument("--force", action="store_true", help="replace generated files")
    args = parser.parse_args(argv)
    output = args.output.resolve()
    if output.exists():
        if not args.force:
            parser.error(f"output already exists: {output}; pass --force to replace it")
        expected = {"report.md", "results.json", "figures"}
        present = {path.name for path in output.iterdir()}
        if not present.issubset(expected):
            parser.error(f"refusing to replace output with unexpected files: {sorted(present - expected)}")
        clear_generated_output(output)
    result = analyze_repository(repo)
    write_outputs(result, output)
    print(output / "report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
