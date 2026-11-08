"""Command-line orchestration for the A8 spectral competition analysis.

Everything is generated into a sibling staging directory and only renamed into
place after source hashes have been re-verified. Sources are opened read-only and
an existing output directory is never deleted or overwritten.

Usage::

    python -m spectral.a8_pipeline.run_analysis \
        --input "A:\\flow_jo_AUG28" \
        --output "A:\\flow_jo_AUG28\\programmatic_analysis" \
        --config ".\\spectral\\a8_pipeline\\default_config.json" [--dry-run]
"""

from __future__ import annotations

import argparse
import json
import shutil
import time
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .fcs_io import load_fcs
from .gating import points_in_polygon, pulse_singlet_mask
from .images import (
    DEFAULT_CHANNEL_MAP,
    build_image_manifest,
    render_montage,
    select_population_images,
)
from .manifest import build_manifest, verify_manifest
from .report import write_report
from .reporters import classify_reporters, equal_posterior_threshold, fit_two_gaussian
from .saturation import decode_saturation
from .sensitivity import (
    SampleInputs,
    build_sensitivity_scenarios,
    default_sensitivity_config,
    run_sensitivity,
    summarise_families,
)
from .statistics import competition_metrics
from .svg import (
    write_pulse_relationship_svg,
    write_reporter_density_svg,
    write_saturation_composition_svg,
    write_scatter_density_svg,
    write_sensitivity_delta_svg,
    write_time_qc_svg,
)
from .time_qc import detect_segments, time_qc

STAGING_NAME = ".programmatic_analysis_staging"

SAMPLE_ROLES = {
    "F11": ("F11-LOY_GFP___CTRL_mCherry.fcs", "F11-LOY_GFP___CTRL_mCherry_images", "GFP"),
    "G11": ("G11-CTRL_GFP__LOY_mCherry.fcs", "G11-CTRL_GFP__LOY_mCherry_images", "mCherry"),
}

CONTROL_FILES = (
    "Plate_002-H07-mCherry_reporter_mCherry_001.fcs",
    "Plate_002-H08-GFP_reporter_eGFP_003.fcs",
    "Plate_002-H11-Unstained_001.fcs",
)

MONTAGE_POPULATIONS = (
    "Scatter-kept",
    "Scatter-excluded",
    "Pulse singlets",
    "Pulse non-singlets",
    "GFP-only",
    "mCherry-only",
    "Double-positive",
    "Double-negative",
    "GFP boundary",
    "mCherry boundary",
    "Spectral saturation",
)


class OutputLocationError(RuntimeError):
    """Raised when the requested output directory is unusable."""


class PipelineValidationError(RuntimeError):
    """Raised when an input or model validity check fails."""


@dataclass
class PipelineResult:
    raw_counts: dict[str, int] = field(default_factory=dict)
    segment_counts: dict[str, list[int]] = field(default_factory=dict)
    matrix_names: dict[str, str | None] = field(default_factory=dict)
    saturation_agreement: dict[str, bool] = field(default_factory=dict)
    image_coverage: dict[str, tuple[int, int]] = field(default_factory=dict)
    thresholds: dict[str, float] = field(default_factory=dict)
    gate_counts: pd.DataFrame | None = None
    competition: pd.DataFrame | None = None
    sensitivity: pd.DataFrame | None = None
    family_summary: pd.DataFrame | None = None
    output: Path | None = None
    staging: Path | None = None
    dry_run: bool = False
    log: list[str] = field(default_factory=list)


def _retry_filesystem(action, description: str, attempts: int = 8, delay: float = 0.75):
    """Retry a filesystem operation that a sync client or scanner may briefly lock.

    The repository lives inside a Dropbox folder, so a directory that was just
    written can be held open by the sync client for a moment. Windows raises
    ``PermissionError`` (WinError 32) rather than blocking, so a short retry is
    the appropriate response.
    """
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            return action()
        except (PermissionError, OSError) as error:
            last = error
            if attempt == attempts - 1:
                break
            time.sleep(delay * (attempt + 1))
    raise RuntimeError(
        f"{description} failed after {attempts} attempts; a sync client or virus "
        f"scanner may be holding it open: {last}"
    ) from last


def resolve_output(input_root: Path, output: Path) -> tuple[Path, Path]:
    """Validate the output location and return it with its staging sibling."""
    source = Path(input_root).resolve()
    destination = Path(output).resolve()

    if destination == source:
        raise OutputLocationError("output must not be the input root itself")
    if destination in source.parents:
        raise OutputLocationError("output must not be a parent of the input root")

    staging = destination.parent / STAGING_NAME
    return destination, staging


def _log(result: PipelineResult, message: str) -> None:
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    line = f"{stamp}  {message}"
    result.log.append(line)
    print(line, flush=True)


def _load_config(config_path: Path) -> dict[str, Any]:
    return json.loads(Path(config_path).read_text(encoding="utf-8"))


def _prepare_sample(
    root: Path, label: str, config: dict, result: PipelineResult
) -> tuple[SampleInputs, pd.DataFrame, np.ndarray, np.ndarray, dict]:
    filename, _image_dir, loy_reporter = SAMPLE_ROLES[label]
    reporters = config["reporter_parameters"]
    required = (reporters["gfp"], reporters["mcherry"], reporters["autofluorescence"])
    data = load_fcs(root / filename, required=required)
    events = data.events

    result.raw_counts[label] = len(events)
    result.matrix_names[label] = data.matrix_name
    if data.matrix_name != config["matrix_name"]:
        raise PipelineValidationError(
            f"{label} matrix name is {data.matrix_name!r}, expected {config['matrix_name']!r}"
        )

    segments = detect_segments(events["Time"].to_numpy(dtype=float))
    result.segment_counts[label] = [segment.size for segment in segments]

    qc = time_qc(
        events,
        segments,
        min_events=config["time_qc"]["minimum_events_per_bin"],
        z_limit=config["time_qc"]["robust_z_limit"],
    )
    saturation = decode_saturation(events)

    decoded_or = saturation["any_saturated"].to_numpy(dtype=bool)
    recorded = events["Saturated"].to_numpy(dtype=float) > 0
    result.saturation_agreement[label] = bool(np.array_equal(decoded_or, recorded))
    if not result.saturation_agreement[label]:
        raise PipelineValidationError(f"{label} decoded saturation mask disagrees with Saturated")

    segment_labels = np.empty(len(events), dtype=object)
    for segment in segments:
        segment_labels[segment.start : segment.stop] = segment.label

    sample = SampleInputs(
        label=label,
        loy_reporter=loy_reporter,
        events=events,
        segment_labels=segment_labels,
        time_good=qc.good_events,
        saturation=saturation,
    )
    return sample, qc.bins, qc.good_events, segment_labels, {"segments": segments}


def run_pipeline(
    input_root: Path,
    output: Path,
    config_path: Path,
    dry_run: bool = False,
) -> PipelineResult:
    """Run the whole analysis. Returns a result even for a dry run."""
    root = Path(input_root).resolve()
    if not root.is_dir():
        raise PipelineValidationError(f"input root does not exist: {root}")

    config = _load_config(config_path)
    destination, staging = resolve_output(root, output)

    result = PipelineResult(dry_run=dry_run, output=destination, staging=staging)
    _log(result, f"input root {root}")
    _log(result, f"output {destination} (dry run: {dry_run})")

    if not dry_run and destination.exists():
        raise OutputLocationError(
            f"output already exists and will not be overwritten: {destination}"
        )
    if staging.exists():
        raise OutputLocationError(f"staging directory already exists: {staging}")

    for control in CONTROL_FILES:
        if not (root / control).is_file():
            raise PipelineValidationError(f"expected control file is missing: {control}")

    _log(result, "hashing sources")
    manifest_before = build_manifest(root)
    _log(result, f"hashed {len(manifest_before)} source files")

    reporters = config["reporter_parameters"]
    gfp_name, mcherry_name = reporters["gfp"], reporters["mcherry"]
    cofactor = float(config["reporters"]["cofactor"])
    polygon = np.asarray(config["cell_gate"]["polygon_raw"], dtype=float)
    gate_x, gate_y = config["cell_gate"]["x"], config["cell_gate"]["y"]
    mad_limit = float(config["singlets"]["primary_mad_limit"])

    samples: list[SampleInputs] = []
    qc_bins: dict[str, pd.DataFrame] = {}
    for label in SAMPLE_ROLES:
        sample, bins, _good, _labels, _extra = _prepare_sample(root, label, config, result)
        samples.append(sample)
        qc_bins[label] = bins
        _log(
            result,
            f"{label}: {result.raw_counts[label]} events, segments {result.segment_counts[label]}, "
            f"saturation decode agrees",
        )

    # ---- primary gating ----
    primary: dict[str, pd.DataFrame] = {}
    pulse_fits: dict[str, dict] = {}
    masks: dict[str, dict[str, np.ndarray]] = {}
    gate_rows: list[dict] = []

    for sample in samples:
        events = sample.events
        inside = points_in_polygon(
            events.loc[:, [gate_x, gate_y]].to_numpy(dtype=float), polygon
        )
        spectral = sample.saturation["spectral_saturated"].to_numpy(dtype=bool)
        pre = sample.time_good & ~spectral & inside
        subset = events.loc[pre]
        singlet_local, fits = pulse_singlet_mask(subset, mad_limit=mad_limit)
        singlet = np.zeros(len(events), dtype=bool)
        singlet[np.flatnonzero(pre)] = singlet_local

        primary[sample.label] = subset.loc[singlet_local]
        pulse_fits[sample.label] = fits
        masks[sample.label] = {
            "inside": inside,
            "spectral": spectral,
            "pre": pre,
            "singlet": singlet,
        }
        gate_rows.append(
            {
                "sample": sample.label,
                "raw_events": len(events),
                "after_time_qc": int(sample.time_good.sum()),
                "after_spectral_saturation": int((sample.time_good & ~spectral).sum()),
                "after_cell_polygon": int(pre.sum()),
                "primary_singlets": int(singlet.sum()),
                "scatter_excluded": int((~inside).sum()),
                "spectral_saturated": int(spectral.sum()),
                "time_qc_excluded": int((~sample.time_good).sum()),
            }
        )
        _log(
            result,
            f"{sample.label}: {int(pre.sum())} pre-singlet, {int(singlet.sum())} primary singlets",
        )

    result.gate_counts = pd.DataFrame(gate_rows)

    # ---- pooled reporter thresholds ----
    models, transformed, thresholds_raw = {}, {}, {}
    for column in (gfp_name, mcherry_name):
        pooled = np.concatenate(
            [primary[s.label][column].to_numpy(dtype=float) for s in samples]
        )
        model = fit_two_gaussian(
            np.arcsinh(pooled / cofactor),
            seed=int(config["reporters"]["seed"]),
            n_init=int(config["reporters"]["initializations"]),
        )
        threshold = equal_posterior_threshold(
            model,
            minimum_weight=float(config["reporters"]["minimum_component_weight"]),
            minimum_separation=float(config["reporters"]["minimum_separation"]),
        )
        models[column] = model
        transformed[column] = threshold
        thresholds_raw[column] = float(np.sinh(threshold) * cofactor)
        _log(
            result,
            f"{column}: threshold {threshold:.6f} transformed = {thresholds_raw[column]:.2f} raw, "
            f"separation {model.separation:.4f}, weights "
            f"{model.weights[0]:.4f}/{model.weights[1]:.4f}",
        )
    result.thresholds = dict(thresholds_raw)

    # ---- classification and competition ----
    classes_by_sample: dict[str, np.ndarray] = {}
    competition_rows: list[dict] = []
    for sample in samples:
        frame = primary[sample.label]
        classes = classify_reporters(
            frame[gfp_name].to_numpy(dtype=float),
            frame[mcherry_name].to_numpy(dtype=float),
            thresholds_raw[gfp_name],
            thresholds_raw[mcherry_name],
        )
        classes_by_sample[sample.label] = classes
        positions = frame.index.to_numpy()
        scopes = ["all"] + list(dict.fromkeys(sample.segment_labels[positions].tolist()))
        for scope in scopes:
            selector = (
                np.ones(len(frame), dtype=bool)
                if scope == "all"
                else sample.segment_labels[positions] == scope
            )
            metrics = competition_metrics(classes[selector], loy_reporter=sample.loy_reporter)
            competition_rows.append(
                {
                    "sample": sample.label,
                    "scope": scope,
                    "loy_reporter": sample.loy_reporter,
                    "events": int(selector.sum()),
                    "gfp_only": metrics.gfp_only,
                    "mcherry_only": metrics.mcherry_only,
                    "double_positive": metrics.double_positive,
                    "double_negative": metrics.double_negative,
                    "loy": metrics.loy,
                    "ctrl": metrics.ctrl,
                    "assignable": metrics.assignable,
                    "loy_fraction": metrics.loy_fraction,
                    "ratio": metrics.ratio,
                    "log2_ratio": metrics.log2_ratio,
                    "wilson_low": metrics.wilson_low,
                    "wilson_high": metrics.wilson_high,
                }
            )
    result.competition = pd.DataFrame(competition_rows)

    dye_swap = result.competition[result.competition["scope"] == "all"]["log2_ratio"]
    dye_swap_mean = float(dye_swap.mean()) if len(dye_swap) == 2 else float("nan")
    _log(result, f"dye-swap mean log2 ratio {dye_swap_mean:.4f}")

    # ---- sensitivity ----
    shift_fraction = float(config["reporters"]["threshold_shift_fraction"])
    cofactor_cache: dict[float, tuple[float, float]] = {}

    def thresholds_for_cofactor(value: float) -> tuple[float, float]:
        if value not in cofactor_cache:
            fitted = []
            for column in (gfp_name, mcherry_name):
                pooled = np.concatenate(
                    [primary[s.label][column].to_numpy(dtype=float) for s in samples]
                )
                model = fit_two_gaussian(
                    np.arcsinh(pooled / value),
                    seed=int(config["reporters"]["seed"]),
                    n_init=int(config["reporters"]["initializations"]),
                )
                fitted.append(float(np.sinh(equal_posterior_threshold(model)) * value))
            cofactor_cache[value] = (fitted[0], fitted[1])
        return cofactor_cache[value]

    def resolve_thresholds(scenario) -> tuple[float, float]:
        if scenario.cofactor != cofactor:
            return thresholds_for_cofactor(scenario.cofactor)
        values = []
        for column, shift in (
            (gfp_name, scenario.gfp_shift),
            (mcherry_name, scenario.mcherry_shift),
        ):
            model = models[column]
            offset = shift * shift_fraction * (model.means[1] - model.means[0])
            values.append(float(np.sinh(transformed[column] + offset) * cofactor))
        return (values[0], values[1])

    scenarios = build_sensitivity_scenarios(default_sensitivity_config())
    _log(result, f"running {len(scenarios)} sensitivity scenarios")
    result.sensitivity = run_sensitivity(samples, scenarios, resolve_thresholds, polygon)

    baseline = result.sensitivity[result.sensitivity["is_baseline"]].set_index("sample")
    summaries = []
    for label in SAMPLE_ROLES:
        subset = result.sensitivity[result.sensitivity["sample"] == label]
        summary = summarise_families(
            subset,
            float(baseline.loc[label, "loy_fraction"]),
            float(baseline.loc[label, "log2_ratio"]),
        )
        summary.insert(0, "sample", label)
        summaries.append(summary)
    result.family_summary = pd.concat(summaries, ignore_index=True)

    # ---- images ----
    image_manifests: dict[str, pd.DataFrame] = {}
    for label in SAMPLE_ROLES:
        _fcs, image_dir, _role = SAMPLE_ROLES[label]
        manifest = build_image_manifest(root / image_dir, event_count=result.raw_counts[label])
        image_manifests[label] = manifest
        coverage = int(manifest["available"].sum())
        result.image_coverage[label] = (coverage, result.raw_counts[label])
        _log(result, f"{label}: image coverage {coverage}/{result.raw_counts[label]}")

    # ---- artifacts ----
    _log(result, f"staging artifacts in {staging}")
    staging.mkdir(parents=True, exist_ok=False)
    try:
        _write_artifacts(
            staging=staging,
            root=root,
            config=config,
            config_path=Path(config_path),
            result=result,
            samples=samples,
            primary=primary,
            masks=masks,
            classes_by_sample=classes_by_sample,
            qc_bins=qc_bins,
            pulse_fits=pulse_fits,
            image_manifests=image_manifests,
            models=models,
            transformed=transformed,
            thresholds_raw=thresholds_raw,
            manifest_before=manifest_before,
            dye_swap_mean=dye_swap_mean,
            gfp_name=gfp_name,
            mcherry_name=mcherry_name,
            cofactor=cofactor,
            polygon=polygon,
            mad_limit=mad_limit,
        )

        _log(result, "re-verifying source hashes")
        verify_manifest(manifest_before)
        _log(result, "source hashes unchanged")

        (staging / "logs").mkdir(exist_ok=True)
        (staging / "logs" / "run.log").write_text("\n".join(result.log), encoding="utf-8")

        if dry_run:
            _retry_filesystem(
                lambda: shutil.rmtree(staging), f"removing staging directory {staging}"
            )
            _log(result, "dry run complete; staging discarded and no output written")
        else:
            _retry_filesystem(
                lambda: staging.rename(destination),
                f"renaming {staging} to {destination}",
            )
            _log(result, f"analysis written to {destination}")
    except Exception:
        if staging.exists():
            shutil.rmtree(staging, ignore_errors=True)
        raise

    return result


def _write_artifacts(
    *,
    staging: Path,
    root: Path,
    config: dict,
    config_path: Path,
    result: PipelineResult,
    samples: list[SampleInputs],
    primary: dict[str, pd.DataFrame],
    masks: dict[str, dict[str, np.ndarray]],
    classes_by_sample: dict[str, np.ndarray],
    qc_bins: dict[str, pd.DataFrame],
    pulse_fits: dict[str, dict],
    image_manifests: dict[str, pd.DataFrame],
    models: dict,
    transformed: dict[str, float],
    thresholds_raw: dict[str, float],
    manifest_before: list,
    dye_swap_mean: float,
    gfp_name: str,
    mcherry_name: str,
    cofactor: float,
    polygon: np.ndarray,
    mad_limit: float,
) -> None:
    figures = staging / "figures"
    for name in ("acquisition", "scatter", "singlets", "reporters", "sensitivity"):
        (figures / name).mkdir(parents=True, exist_ok=True)
    montage_dir = staging / "image_montages"
    montage_dir.mkdir(parents=True, exist_ok=True)

    # frozen configuration, including the accepted thresholds
    frozen = dict(config)
    frozen["frozen_thresholds_raw"] = thresholds_raw
    frozen["frozen_thresholds_transformed"] = transformed
    frozen["reporter_models"] = {
        column: {
            "means": [float(v) for v in model.means],
            "sigmas": [float(v) for v in model.sigmas],
            "weights": [float(v) for v in model.weights],
            "separation": float(model.separation),
            "log_likelihood": float(model.log_likelihood),
        }
        for column, model in models.items()
    }
    frozen["source_config_file"] = str(config_path)
    (staging / "analysis_config.json").write_text(
        json.dumps(frozen, indent=2), encoding="utf-8"
    )

    pd.DataFrame(
        [
            {
                "path": str(record.path.relative_to(root)),
                "size": record.size,
                "sha256": record.sha256,
            }
            for record in manifest_before
        ]
    ).to_csv(staging / "input_manifest.csv", index=False)

    result.gate_counts.to_csv(staging / "gate_counts.csv", index=False)
    result.competition.to_csv(staging / "competition_results.csv", index=False)
    result.sensitivity.to_csv(staging / "sensitivity_results.csv", index=False)
    result.family_summary.to_csv(staging / "sensitivity_family_summary.csv", index=False)

    # ---- event labels ----
    label_frames = []
    for sample in samples:
        events = sample.events
        mask = masks[sample.label]
        positions = primary[sample.label].index.to_numpy()
        reporter_class = np.full(len(events), "not-classified", dtype=object)
        reporter_class[positions] = classes_by_sample[sample.label]
        manifest = image_manifests[sample.label]
        label_frames.append(
            pd.DataFrame(
                {
                    "sample": sample.label,
                    "source_file": SAMPLE_ROLES[sample.label][0],
                    "event_number": np.arange(len(events), dtype=int) + 1,
                    "segment": sample.segment_labels,
                    "time_qc_pass": sample.time_good,
                    "spectral_saturated": sample.saturation["spectral_saturated"].to_numpy(),
                    "any_saturated": sample.saturation["any_saturated"].to_numpy(),
                    "imgb3_saturated": sample.saturation["imgb3_saturated"].to_numpy(),
                    "in_cell_gate": mask["inside"],
                    "pulse_singlet": mask["singlet"],
                    "reporter_class": reporter_class,
                    "image_available": manifest["available"].to_numpy(),
                    gfp_name: events[gfp_name].to_numpy(dtype=float),
                    mcherry_name: events[mcherry_name].to_numpy(dtype=float),
                }
            )
        )
    labels = pd.concat(label_frames, ignore_index=True)
    if labels.duplicated(subset=["sample", "event_number"]).any():
        raise PipelineValidationError("event labels are not unique by sample and event number")
    labels.to_csv(staging / "event_labels.csv.gz", index=False, compression="gzip")

    image_rows = []
    for label, manifest in image_manifests.items():
        frame = manifest.copy()
        frame.insert(0, "sample", label)
        frame["path"] = [None if p is None else Path(p).name for p in frame["path"]]
        image_rows.append(frame)
    pd.concat(image_rows, ignore_index=True).to_csv(staging / "image_manifest.csv", index=False)

    # ---- figures ----
    saturation_frame = pd.DataFrame(
        [
            {
                "sample": sample.label,
                "any_saturated": 100.0 * float(sample.saturation["any_saturated"].mean()),
                "spectral_saturated": 100.0 * float(sample.saturation["spectral_saturated"].mean()),
                "imgb3_saturated": 100.0 * float(sample.saturation["imgb3_saturated"].mean()),
            }
            for sample in samples
        ]
    )
    write_saturation_composition_svg(
        saturation_frame, figures / "acquisition" / "saturation_composition.svg",
        "Saturation composition (% of events)",
    )

    for sample in samples:
        write_time_qc_svg(
            qc_bins[sample.label],
            figures / "acquisition" / f"time_qc_{sample.label}.svg",
            f"{sample.label} binned event rate",
        )
        write_scatter_density_svg(
            sample.events, polygon,
            figures / "scatter" / f"scatter_{sample.label}.svg",
            f"{sample.label} {config['cell_gate']['x']} vs {config['cell_gate']['y']}",
        )
        for axis, area, height in (
            ("fsc", "FSC-A", "FSC-H"),
            ("ssc", "SSC (Violet)-A", "SSC (Violet)-H"),
        ):
            write_pulse_relationship_svg(
                sample.events[area].to_numpy(dtype=float),
                sample.events[height].to_numpy(dtype=float),
                pulse_fits[sample.label][axis],
                mad_limit,
                figures / "singlets" / f"pulse_{axis}_{sample.label}.svg",
                f"{sample.label} {area} vs {height}",
            )

    for column in (gfp_name, mcherry_name):
        pooled = np.concatenate(
            [primary[s.label][column].to_numpy(dtype=float) for s in samples]
        )
        write_reporter_density_svg(
            np.arcsinh(pooled / cofactor),
            transformed[column],
            figures / "reporters" / f"reporter_{column.replace('-', '_')}.svg",
            f"Pooled {column} with accepted threshold",
        )

    for label in SAMPLE_ROLES:
        subset = result.family_summary[result.family_summary["sample"] == label]
        write_sensitivity_delta_svg(
            subset.drop(columns=["sample"]),
            figures / "sensitivity" / f"sensitivity_delta_{label}.svg",
            f"{label} sensitivity: max |change| in LOY fraction",
        )

    # ---- montages ----
    montage_rows = []
    seed = int(config["images"]["seed"])
    maximum = int(config["images"]["maximum_per_population"])
    for sample in samples:
        manifest = image_manifests[sample.label]
        mask = masks[sample.label]
        events = sample.events
        numbers = np.arange(len(events), dtype=int) + 1
        positions = primary[sample.label].index.to_numpy()
        classes = classes_by_sample[sample.label]

        tg = np.arcsinh(events[gfp_name].to_numpy(dtype=float) / cofactor)
        tm = np.arcsinh(events[mcherry_name].to_numpy(dtype=float) / cofactor)
        d_gfp = np.abs(tg - transformed[gfp_name])
        d_mch = np.abs(tm - transformed[mcherry_name])

        groups = {
            "Scatter-kept": numbers[mask["inside"]],
            "Scatter-excluded": numbers[~mask["inside"]],
            "Pulse singlets": numbers[mask["singlet"]],
            "Pulse non-singlets": numbers[mask["pre"] & ~mask["singlet"]],
            "GFP-only": numbers[positions[classes == "GFP-only"]],
            "mCherry-only": numbers[positions[classes == "mCherry-only"]],
            "Double-positive": numbers[positions[classes == "Double-positive"]],
            "Double-negative": numbers[positions[classes == "Double-negative"]],
            "GFP boundary": numbers[positions[np.argsort(d_gfp[positions])[:maximum]]],
            "mCherry boundary": numbers[positions[np.argsort(d_mch[positions])[:maximum]]],
            "Spectral saturation": numbers[mask["spectral"]],
        }

        selections = {}
        for population in MONTAGE_POPULATIONS:
            picked = select_population_images(
                manifest, groups[population], maximum=maximum, seed=seed
            )
            if picked:
                selections[population] = picked
            for event, path in picked:
                index = event - 1
                montage_rows.append(
                    {
                        "sample": sample.label,
                        "population": population,
                        "event_number": event,
                        "image_file": path.name,
                        gfp_name: float(events[gfp_name].to_numpy(dtype=float)[index]),
                        mcherry_name: float(events[mcherry_name].to_numpy(dtype=float)[index]),
                        "gfp_threshold_distance": float(d_gfp[index]),
                        "mcherry_threshold_distance": float(d_mch[index]),
                        "any_saturated": bool(sample.saturation["any_saturated"].to_numpy()[index]),
                        "spectral_saturated": bool(
                            sample.saturation["spectral_saturated"].to_numpy()[index]
                        ),
                        "imgb3_saturated": bool(
                            sample.saturation["imgb3_saturated"].to_numpy()[index]
                        ),
                    }
                )
        if selections:
            render_montage(
                selections,
                montage_dir / f"montage_{sample.label}_full.png",
                channel_map=DEFAULT_CHANNEL_MAP,
            )
            render_montage(
                {name: rows[:8] for name, rows in selections.items()},
                montage_dir / f"montage_{sample.label}_review.png",
                channel_map=DEFAULT_CHANNEL_MAP,
            )
    pd.DataFrame(montage_rows).to_csv(montage_dir / "image_montage_manifest.csv", index=False)

    # ---- report ----
    coverage_note = ", ".join(
        f"{label} {seen}/{total}" for label, (seen, total) in result.image_coverage.items()
    )
    figure_list = [
        {"caption": path.stem.replace("_", " "), "path": str(path.relative_to(staging)).replace("\\", "/")}
        for path in sorted(figures.rglob("*.svg"))
    ]
    write_report(
        {
            "title": "A8 spectral competition analysis",
            "summary": (
                "Read-only analysis of the reciprocal dye-swap pair in "
                f"{root}. Image coverage: {coverage_note}. "
                f"Exploratory dye-swap mean log2 ratio: {dye_swap_mean:.4f}."
            ),
            "interval_note": (
                "Wilson intervals describe counting precision for the events acquired. "
                "They are not biological replicate uncertainty, because there is one "
                "acquisition per dye orientation and no independent culture replicates."
            ),
            "tables": {
                "Gate counts": result.gate_counts,
                "Competition results": result.competition,
                "Sensitivity family summary": result.family_summary,
                "Saturation composition (%)": saturation_frame,
            },
            "figures": figure_list,
            "limitations": _limitations(result, dye_swap_mean),
        },
        staging,
    )

    (staging / "README.md").write_text(_readme(result, root), encoding="utf-8")


def _limitations(result: PipelineResult, dye_swap_mean: float) -> list[str]:
    f11_seen, f11_total = result.image_coverage.get("F11", (0, 0))
    competition = result.competition
    overall = competition[competition["scope"] == "all"].set_index("sample")
    return [
        "No viability dye was acquired, so scatter-excluded events are not dead cells and "
        "the scatter gate is not a live-cell gate.",
        "The reciprocal dye swap does not reproduce: F11 gives a log2 ratio of "
        f"{overall.loc['F11', 'log2_ratio']:.4f} and G11 gives "
        f"{overall.loc['G11', 'log2_ratio']:.4f}. In both orientations the mCherry-labelled "
        "arm is larger regardless of genotype, so a single-orientation result reports "
        "reporter bias as well as biology. The exploratory dye-swap mean is "
        f"{dye_swap_mean:.4f}.",
        "Double-positive events are excluded from the assignable denominator and are "
        f"{100.0 * overall.loc['F11', 'double_positive'] / overall.loc['F11', 'events']:.1f}% of F11 and "
        f"{100.0 * overall.loc['G11', 'double_positive'] / overall.loc['G11', 'events']:.1f}% of G11 "
        "primary events. Image review shows they are a mixture of unresolved GFP/mCherry "
        "doublets and single cells with genuine low mCherry.",
        f"F11 image coverage is {f11_seen}/{f11_total} and is not missing at random: no "
        "scatter-excluded event has an image and events without images are larger. F11 "
        "montages cannot support quantitative frequency claims.",
        "The saturation sensitivity family dominates all others. A blanket "
        "Saturated == 0 rule discards roughly half of each recording because of ImgB3 (790) "
        "imaging clipping and is not a defensible primary rule for a spectral ratio.",
        "The single-stain controls do not carry the unmixed reporter parameters, so they "
        "cannot anchor or validate the reporter thresholds in unmixed space.",
        "There is one acquisition per dye orientation and no independent culture "
        "replicates, so this is a technical and QC analysis, not evidence of biological "
        "reproducibility.",
    ]


def _readme(result: PipelineResult, root: Path) -> str:
    lines = [
        "# Programmatic analysis output",
        "",
        f"Generated from `{root}` by `spectral.a8_pipeline.run_analysis`.",
        "",
        "All source FCS, CVW, TIFF and FlowJo files were opened read-only. SHA-256 hashes",
        "were recorded before the run and re-verified afterwards; the run fails if any",
        "source changed.",
        "",
        "## Contents",
        "",
        "- `analysis_config.json` — the configuration used, plus the frozen reporter thresholds and fitted mixture parameters.",
        "- `input_manifest.csv` — source files with sizes and SHA-256 hashes.",
        "- `gate_counts.csv` — nested gate counts per recording.",
        "- `event_labels.csv.gz` — one row per event with segment, gate flags, saturation flags, reporter class and image availability.",
        "- `competition_results.csv` — LOY:CTRL per recording and per acquisition segment.",
        "- `sensitivity_results.csv` — every predefined scenario.",
        "- `sensitivity_family_summary.csv` — maximum absolute change per family.",
        "- `image_manifest.csv` — every event mapped to its TIFF, with missing images explicit.",
        "- `report.md` / `report.html` — the written report.",
        "- `figures/` — acquisition, scatter, singlet, reporter and sensitivity SVGs.",
        "- `image_montages/` — stratified review montages and their companion CSV.",
        "- `logs/run.log` — the run log.",
        "",
        "## Reporter thresholds in raw units, for FlowJo cross-check",
        "",
    ]
    for column, value in result.thresholds.items():
        lines.append(f"- `{column}`: {value:.4f}")
    lines.extend(
        [
            "",
            "Image filename suffixes are ZERO-based event indices. `event_number` in every",
            "table here is one-based, so `event_number = suffix + 1`.",
            "",
        ]
    )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m spectral.a8_pipeline.run_analysis",
        description="Run the A8 spectral competition analysis read-only.",
    )
    parser.add_argument("--input", required=True, type=Path, help="input root directory")
    parser.add_argument("--output", required=True, type=Path, help="output directory to create")
    parser.add_argument("--config", required=True, type=Path, help="analysis configuration JSON")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="validate and compute everything but write no output directory",
    )
    args = parser.parse_args(argv)

    try:
        run_pipeline(args.input, args.output, args.config, dry_run=args.dry_run)
    except (OutputLocationError, PipelineValidationError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
