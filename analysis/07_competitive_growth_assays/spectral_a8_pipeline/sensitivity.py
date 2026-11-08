"""Predefined sensitivity families for the A8 competition analysis.

The matrix is fixed in advance and crosses only interpretable alternatives. It is
deliberately not an exhaustive product of every option, and no scenario is
selected or rejected by reference to the LOY:CTRL result it produces.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from .gating import points_in_polygon, pulse_singlet_mask, scale_polygon
from .reporters import classify_reporters
from .statistics import competition_metrics

SATURATION_RULES = ("spectral", "none", "relevant", "any")


@dataclass(frozen=True)
class Scenario:
    """One predefined analysis alternative."""

    family: str
    label: str
    time_qc: bool = True
    segment: str = "all"
    saturation_rule: str = "spectral"
    polygon_scale: float = 1.0
    singlet_mad: float = 3.5
    cofactor: float = 10000.0
    gfp_shift: float = 0.0
    mcherry_shift: float = 0.0
    is_baseline: bool = False

    def __post_init__(self) -> None:
        if self.saturation_rule not in SATURATION_RULES:
            raise ValueError(f"unknown saturation rule: {self.saturation_rule}")
        if self.polygon_scale <= 0:
            raise ValueError("polygon scale must be positive")
        if self.singlet_mad <= 0:
            raise ValueError("singlet MAD limit must be positive")
        if self.cofactor <= 0:
            raise ValueError("cofactor must be positive")


@dataclass(frozen=True)
class SampleInputs:
    """Everything a scenario needs from one recording, computed once."""

    label: str
    loy_reporter: str
    events: pd.DataFrame
    segment_labels: np.ndarray
    time_good: np.ndarray
    saturation: pd.DataFrame


def default_sensitivity_config() -> dict:
    """Return the frozen sensitivity matrix described in the design document."""
    return {
        "time_qc": [True, False],
        "segments": ["Run 1", "Run 2"],
        "polygon_scales": [0.95, 1.05],
        "singlet_mads": [3.0, 4.0],
        "saturation_rules": ["none", "relevant", "any"],
        "threshold_shifts": [-1.0, 0.0, 1.0],
        "cofactors": [5000.0, 20000.0],
    }


def build_sensitivity_scenarios(config: Mapping[str, Sequence]) -> list[Scenario]:
    """Expand the configuration into the fixed list of scenarios.

    Exactly one scenario is marked ``is_baseline``. Every other family lists only
    the alternatives that differ from it, except the threshold family, which
    includes its unshifted centre as an internal consistency check.
    """
    scenarios: list[Scenario] = [
        Scenario(family="time", label="baseline: time QC on", is_baseline=True)
    ]

    for enabled in config["time_qc"]:
        if not enabled:
            scenarios.append(Scenario(family="time", label="time QC off", time_qc=False))

    for segment in config["segments"]:
        scenarios.append(Scenario(family="segment", label=f"segment {segment}", segment=segment))

    for scale in config["polygon_scales"]:
        scenarios.append(
            Scenario(family="cell", label=f"cell polygon x{scale:g}", polygon_scale=float(scale))
        )

    for mad in config["singlet_mads"]:
        scenarios.append(
            Scenario(family="singlet", label=f"singlet {mad:g} MAD", singlet_mad=float(mad))
        )

    for rule in config["saturation_rules"]:
        scenarios.append(
            Scenario(family="saturation", label=f"saturation rule {rule}", saturation_rule=rule)
        )

    for gfp in config["threshold_shifts"]:
        for mcherry in config["threshold_shifts"]:
            scenarios.append(
                Scenario(
                    family="threshold",
                    label=f"threshold shift GFP {gfp:+g} mCherry {mcherry:+g}",
                    gfp_shift=float(gfp),
                    mcherry_shift=float(mcherry),
                )
            )

    for cofactor in config["cofactors"]:
        scenarios.append(
            Scenario(
                family="cofactor", label=f"cofactor {cofactor:g}", cofactor=float(cofactor)
            )
        )

    return scenarios


def saturation_exclusion(saturation: pd.DataFrame, rule: str) -> np.ndarray:
    """Return the boolean mask of events a saturation rule excludes."""
    if rule == "none":
        return np.zeros(len(saturation), dtype=bool)
    if rule == "spectral":
        return saturation["spectral_saturated"].to_numpy(dtype=bool)
    if rule == "relevant":
        return (
            saturation["spectral_saturated"].to_numpy(dtype=bool)
            | saturation["imgb1_saturated"].to_numpy(dtype=bool)
            | saturation["imgb2_saturated"].to_numpy(dtype=bool)
        )
    if rule == "any":
        return saturation["any_saturated"].to_numpy(dtype=bool)
    raise ValueError(f"unknown saturation rule: {rule}")


def evaluate_scenario(
    sample: SampleInputs,
    scenario: Scenario,
    thresholds: tuple[float, float],
    polygon: np.ndarray,
    gate_x: str = "FSC-A",
    gate_y: str = "SSC (Violet)-A",
) -> dict:
    """Apply one scenario to one recording and return its counts and metrics."""
    events = sample.events
    total = len(events)

    keep = sample.time_good.copy() if scenario.time_qc else np.ones(total, dtype=bool)
    after_time = int(keep.sum())

    keep &= ~saturation_exclusion(sample.saturation, scenario.saturation_rule)
    after_saturation = int(keep.sum())

    vertices = (
        polygon if scenario.polygon_scale == 1.0 else scale_polygon(polygon, scenario.polygon_scale)
    )
    inside = points_in_polygon(events.loc[:, [gate_x, gate_y]].to_numpy(dtype=float), vertices)
    keep &= inside
    after_polygon = int(keep.sum())

    if scenario.segment != "all":
        keep &= sample.segment_labels == scenario.segment
    after_segment = int(keep.sum())

    result = {
        "sample": sample.label,
        "family": scenario.family,
        "label": scenario.label,
        "is_baseline": scenario.is_baseline,
        "time_qc": scenario.time_qc,
        "segment": scenario.segment,
        "saturation_rule": scenario.saturation_rule,
        "polygon_scale": scenario.polygon_scale,
        "singlet_mad": scenario.singlet_mad,
        "cofactor": scenario.cofactor,
        "gfp_shift": scenario.gfp_shift,
        "mcherry_shift": scenario.mcherry_shift,
        "gfp_threshold": thresholds[0],
        "mcherry_threshold": thresholds[1],
        "events_total": total,
        "after_time_qc": after_time,
        "after_saturation": after_saturation,
        "after_polygon": after_polygon,
        "after_segment": after_segment,
    }

    if after_segment < 3:
        result.update(
            {
                "singlets": 0,
                "gfp_only": 0,
                "mcherry_only": 0,
                "double_positive": 0,
                "double_negative": 0,
                "assignable": 0,
                "loy_fraction": float("nan"),
                "ratio": float("nan"),
                "log2_ratio": float("nan"),
                "wilson_low": float("nan"),
                "wilson_high": float("nan"),
            }
        )
        return result

    subset = events.loc[keep]
    singlet, _ = pulse_singlet_mask(subset, mad_limit=scenario.singlet_mad)
    primary = subset.loc[singlet]

    classes = classify_reporters(
        primary["eGFP-A"].to_numpy(dtype=float),
        primary["mCherry-A"].to_numpy(dtype=float),
        thresholds[0],
        thresholds[1],
    )
    metrics = competition_metrics(classes, loy_reporter=sample.loy_reporter)
    result.update(
        {
            "singlets": int(singlet.sum()),
            "gfp_only": metrics.gfp_only,
            "mcherry_only": metrics.mcherry_only,
            "double_positive": metrics.double_positive,
            "double_negative": metrics.double_negative,
            "assignable": metrics.assignable,
            "loy_fraction": metrics.loy_fraction,
            "ratio": metrics.ratio,
            "log2_ratio": metrics.log2_ratio,
            "wilson_low": metrics.wilson_low,
            "wilson_high": metrics.wilson_high,
        }
    )
    return result


def run_sensitivity(
    samples: Iterable[SampleInputs],
    scenarios: Sequence[Scenario],
    threshold_resolver: Callable[[Scenario], tuple[float, float]],
    polygon: np.ndarray,
) -> pd.DataFrame:
    """Run every scenario against every recording.

    ``threshold_resolver`` maps a scenario to the raw-unit GFP and mCherry
    thresholds it should use, so cofactor refits and directional shifts stay out
    of the gating logic.
    """
    rows = []
    for sample in samples:
        for scenario in scenarios:
            if scenario.segment != "all" and scenario.segment not in set(sample.segment_labels):
                continue
            rows.append(
                evaluate_scenario(sample, scenario, threshold_resolver(scenario), polygon)
            )
    return pd.DataFrame(rows)


def summarise_families(
    results: pd.DataFrame,
    baseline_fraction: float,
    baseline_log2: float,
) -> pd.DataFrame:
    """Report the maximum absolute change per sensitivity family."""
    frame = results.copy()
    frame["delta_loy_fraction"] = frame["loy_fraction"] - baseline_fraction
    frame["delta_log2_ratio"] = frame["log2_ratio"] - baseline_log2
    grouped = frame.groupby("family", sort=False)

    def _peak(values: pd.Series) -> float:
        magnitudes = np.abs(values.to_numpy(dtype=float))
        if not np.any(np.isfinite(magnitudes)):
            return float("nan")
        return float(np.nanmax(magnitudes))

    summary = grouped.agg(
        scenarios=("delta_loy_fraction", "size"),
        max_abs_delta_loy_fraction=("delta_loy_fraction", _peak),
        max_abs_delta_log2_ratio=("delta_log2_ratio", _peak),
    ).reset_index()
    return summary


def shifted_threshold(
    transformed_threshold: float,
    mean_low: float,
    mean_high: float,
    shift: float,
    shift_fraction: float,
    cofactor: float,
) -> float:
    """Move a threshold by a fraction of the distance between component means."""
    offset = shift * shift_fraction * (mean_high - mean_low)
    return float(np.sinh(transformed_threshold + offset) * cofactor)


def relabel_baseline(scenarios: Sequence[Scenario], family: str) -> list[Scenario]:
    """Return the scenarios with the baseline reassigned to another family."""
    out = []
    for scenario in scenarios:
        out.append(replace(scenario, family=family) if scenario.is_baseline else scenario)
    return out
