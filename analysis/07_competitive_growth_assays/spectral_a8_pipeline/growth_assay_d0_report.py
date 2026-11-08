r"""D0 baseline report for the 15-Oct-2026 A8 knockout competition growth assay.

The export ``20261016_171604_Export/comp_growth_assay_D0_OCT_15`` holds one
plate (Plate_004) of LOY-versus-CTRL mixtures in four knockout backgrounds
(WT, DDX3X, SREBF1, TADA2B), both reciprocal dye orientations, three replicate
wells each, plus the four single-colour SOLO controls per background and two
same-genotype GFP:mCherry dye-neutral mixtures.  ``_C1`` wells (WT and DDX3X)
are the plates that will receive the C1 DDX3X inhibitor; at D0 they are untreated.

Raw FCS, CEF and statistics files are opened read-only.  The method is the one
used for the August/September runs so that later timepoints stay comparable:
the run's own FACSChorus Scatter -> SSC Singlets -> FSC Singlets polygons,
``arcsinh(value / 1000)`` on the unmixed eGFP-A and mCherry-A parameters, and an
equal-prior two-class LDA trained on this run's SOLO controls.

Reporter-negative events are removed before the ratio is taken (the ``broad``
rule).  No unstained control was on this plate, so the rule is checked against
the 1-Sep-2026 unstained wells on the raw 625 nm detector (YG4), which is
upstream of unmixing: events removed by the rule should sit below the
unstained 99th percentile, and the retained mCherry class above it.

Run from the repository root::

    .\codex\flow_env\python.exe -m spectral.a8_pipeline.growth_assay_d0_report
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from hashlib import sha256
from html import escape
import json
import math
from pathlib import Path
import re
from typing import Iterable, Mapping

import numpy as np

from .competition_report import (
    LDAModel,
    _boundary_segment,
    _file_map,
    _primary_features,
    _sample_features,
    fit_lda,
    load_cef_gates,
)
from .fcs_io import load_fcs
from .gating import points_in_polygon
from .longitudinal_report import verify_unmixing_compatibility

BACKGROUNDS = ("WT", "DDX3X", "SREBF1", "TADA2B")
ORIENTATIONS = ("LOY_GFP", "LOY_mCherry")
GFP_COLOUR = "#1f7a4d"
MCHERRY_COLOUR = "#a8329b"
MUTED = "#6b7280"

# Reporter-negative rectangles in transformed units.  The strict rule captures
# only the unambiguous double-negative cluster; the broad rule also removes the
# GFP-negative events whose mCherry sits in the unmixing-noise band (about +5),
# where reporter-lost cells and mCherry-dim cells cannot be told apart.
NEGATIVE_RULES = {
    "strict": (6.0, 2.0),
    "broad": (6.0, 6.5),
}
# The primary analysis excludes broad negatives (chosen 2026-10-16 after the
# unstained check showed the mCherry ~+5 band is reporter-negative).  Every
# later timepoint of this assay must use the same rule.
PRIMARY_NEGATIVE_RULE = "broad"
SENSITIVITY_NEGATIVE_RULE = "strict"
# Raw detectors used for the unstained check (not affected by unmixing).
RAW_GFP, RAW_MCHERRY = "B3 (530)-A", "YG4 (625)-A"
UNSTAINED_PERCENTILE = 99.0

_SUFFIX = r"(?:_(?P<record>\d{3}))?"
_MIXTURE = re.compile(
    r"^(?P<well>[A-H]\d{2})-LOY_(?P<dye>CHERRY|GFP)_(?P<background>WT|DDX3X|SREBF1|TADA2B)"
    r"(?:_(?P<condition>C1))?_(?P<replicate>[1-3])" + _SUFFIX + "$"
)
_SOLO = re.compile(
    r"^(?P<well>[A-H]\d{2})-(?P<genotype>LOY|CTRL)_(?P<dye>CHERRY|GFP)_SOLO_"
    r"(?P<background>WT|DDX3X|SREBF1|TADA2B)" + _SUFFIX + "$"
)
_DYE_NEUTRAL = re.compile(
    r"^(?P<well>[A-H]\d{2})-(?P<genotype>LOY|CTRL)_GFP_CHERRY_(?P<replicate>[1-3])" + _SUFFIX + "$"
)


@dataclass(frozen=True)
class GrowthSample:
    name: str
    well: str
    role: str  # "mixture", "solo" or "dye_neutral"
    background: str
    condition: str  # "standard" or "C1" (plate to receive the C1 DDX3X inhibitor)
    genotype: str | None  # SOLO and dye-neutral wells only
    dye: str | None  # SOLO: the reporter; mixture: the reporter carried by LOY
    replicate: int | None

    @property
    def orientation(self) -> str | None:
        if self.role != "mixture":
            return None
        return "LOY_GFP" if self.dye == "GFP" else "LOY_mCherry"

    @property
    def group(self) -> str:
        """Mixture condition label, e.g. ``DDX3X`` or ``WT_C1``."""
        return self.background if self.condition == "standard" else f"{self.background}_{self.condition}"


def _dye(text: str) -> str:
    return "GFP" if text == "GFP" else "mCherry"


def parse_sample(name: str) -> GrowthSample:
    """Parse one Plate_004 recording name (FCS stem) into design metadata."""
    match = _MIXTURE.match(name)
    if match:
        return GrowthSample(
            name=name,
            well=match["well"],
            role="mixture",
            background=match["background"],
            condition=match["condition"] or "standard",
            genotype=None,
            dye=_dye(match["dye"]),
            replicate=int(match["replicate"]),
        )
    match = _SOLO.match(name)
    if match:
        return GrowthSample(
            name=name,
            well=match["well"],
            role="solo",
            background=match["background"],
            condition="standard",
            genotype=match["genotype"],
            dye=_dye(match["dye"]),
            replicate=None,
        )
    match = _DYE_NEUTRAL.match(name)
    if match:
        # The dye-neutral wells carry no knockout label; they sit with the WT block.
        return GrowthSample(
            name=name,
            well=match["well"],
            role="dye_neutral",
            background="WT",
            condition="standard",
            genotype=match["genotype"],
            dye=None,
            replicate=int(match["replicate"]),
        )
    raise ValueError(f"unrecognised growth-assay sample name: {name}")


def build_design(names: Iterable[str]) -> list[GrowthSample]:
    """Parse every recording and require complete replicate and control sets."""
    samples = sorted((parse_sample(name) for name in names), key=lambda sample: sample.well)
    wells = [sample.well for sample in samples]
    if len(set(wells)) != len(wells):
        raise ValueError("duplicate well identifiers in the export")
    groups: dict[tuple, list[int]] = {}
    for sample in samples:
        if sample.role == "mixture":
            groups.setdefault((sample.group, sample.orientation), []).append(sample.replicate)
        elif sample.role == "dye_neutral":
            groups.setdefault(("dye_neutral", sample.genotype), []).append(sample.replicate)
    for key, replicates in groups.items():
        if sorted(replicates) != [1, 2, 3]:
            raise ValueError(f"expected replicates 1-3 for {key}; found {sorted(replicates)}")
    for background in {sample.background for sample in samples if sample.role == "mixture"}:
        controls = {
            (sample.genotype, sample.dye)
            for sample in samples
            if sample.role == "solo" and sample.background == background
        }
        expected = {(g, d) for g in ("LOY", "CTRL") for d in ("GFP", "mCherry")}
        if controls != expected:
            raise ValueError(f"{background} is missing SOLO controls: {sorted(expected - controls)}")
    return samples


def read_instrument_singlets(path: Path) -> dict[str, int]:
    """Return the FACSChorus FSC Singlets count for every recording."""
    counts: dict[str, int] = {}
    current = None
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle):
            if row and row[0] == "Recording Name":
                current = row[1]
            elif current and row and row[0] == "FSC Singlets":
                counts[current] = int(row[1])
    return counts


def negative_mask(features: np.ndarray, rule: str) -> np.ndarray:
    gfp_limit, mcherry_limit = NEGATIVE_RULES[rule]
    values = np.asarray(features, dtype=float)
    return (values[:, 0] < gfp_limit) & (values[:, 1] < mcherry_limit)


def logit(fraction: float) -> float:
    return math.log(fraction / (1.0 - fraction))


def inverse_logit(value: float) -> float:
    return 1.0 / (1.0 + math.exp(-value))


def measure_well(
    features: np.ndarray,
    classifier: LDAModel,
    reference: LDAModel,
    sample: GrowthSample,
    all_events_classifier: LDAModel | None = None,
) -> dict:
    """Classify one well under the primary rule and every sensitivity rule.

    The primary fractions (``gfp_fraction``, ``loy_fraction``) count only
    reporter-positive events, i.e. those outside the primary negative rule.
    ``*_all_events`` keeps every event and uses ``all_events_classifier``
    (trained on every SOLO event), which is the August/September method.
    """
    gfp = classifier.predict(features).astype(bool)
    reference_gfp = reference.predict(features).astype(bool)
    all_gfp = (all_events_classifier or classifier).predict(features).astype(bool)
    negative = {rule: negative_mask(features, rule) for rule in NEGATIVE_RULES}
    positive = ~negative[PRIMARY_NEGATIVE_RULE]
    sensitivity = ~negative[SENSITIVITY_NEGATIVE_RULE]
    row: dict = {
        "primary_singlets": int(len(features)),
        "reporter_positive_singlets": int(positive.sum()),
        "gfp_count": int(gfp[positive].sum()),
        "gfp_fraction": float(gfp[positive].mean()),
        "gfp_fraction_september_classifier": float(reference_gfp[positive].mean()),
        "gfp_fraction_all_events": float(all_gfp.mean()),
        f"gfp_fraction_excluding_{SENSITIVITY_NEGATIVE_RULE}_negatives": float(gfp[sensitivity].mean()),
    }
    for rule in NEGATIVE_RULES:
        row[f"negative_{rule}_fraction"] = float(negative[rule].mean())
    if sample.role == "mixture":
        flip = sample.orientation == "LOY_mCherry"
        for key in GFP_KEYS:
            target = key.replace("gfp_fraction", "loy_fraction")
            row[target] = 1.0 - row[key] if flip else row[key]
        positives = row["reporter_positive_singlets"]
        row["loy_count"] = positives - row["gfp_count"] if flip else row["gfp_count"]
    return row


def _mean_sd(values: list[float]) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    return float(array.mean()), float(array.std(ddof=1)) if len(array) > 1 else 0.0


GFP_KEYS = (
    "gfp_fraction",
    "gfp_fraction_september_classifier",
    "gfp_fraction_all_events",
    f"gfp_fraction_excluding_{SENSITIVITY_NEGATIVE_RULE}_negatives",
)
LOY_KEYS = tuple(key.replace("gfp_fraction", "loy_fraction") for key in GFP_KEYS)


def summarise_conditions(wells: list[dict]) -> tuple[list[dict], list[dict]]:
    """Summarise replicate wells per orientation, then balance the dye swap."""
    by_orientation: dict[tuple[str, str], list[dict]] = {}
    for well in wells:
        if well["role"] == "mixture":
            by_orientation.setdefault((well["group"], well["orientation"]), []).append(well)
    orientation_rows = []
    for (group, orientation), rows in sorted(
        by_orientation.items(), key=lambda item: (_group_order(item[0][0]), item[0][1])
    ):
        summary = {
            "group": group,
            "background": rows[0]["background"],
            "condition": rows[0]["condition"],
            "orientation": orientation,
            "replicates": len(rows),
            "wells": ";".join(row["well"] for row in rows),
            "mean_primary_singlets": float(np.mean([row["primary_singlets"] for row in rows])),
            "mean_reporter_positive_singlets": float(
                np.mean([row["reporter_positive_singlets"] for row in rows])
            ),
        }
        for key in LOY_KEYS:
            mean, sd = _mean_sd([row[key] for row in rows])
            summary[f"{key}_mean"], summary[f"{key}_sd"] = mean, sd
        summary["negative_broad_fraction_mean"] = float(
            np.mean([row["negative_broad_fraction"] for row in rows])
        )
        orientation_rows.append(summary)
    balanced_rows = []
    groups = sorted({row["group"] for row in orientation_rows}, key=_group_order)
    for group in groups:
        pair = {row["orientation"]: row for row in orientation_rows if row["group"] == group}
        if set(pair) != set(ORIENTATIONS):
            continue
        balanced = {"group": group}
        for key in LOY_KEYS:
            gfp, mcherry = pair["LOY_GFP"][f"{key}_mean"], pair["LOY_mCherry"][f"{key}_mean"]
            balanced[f"{key}_balanced"] = inverse_logit(0.5 * (logit(gfp) + logit(mcherry)))
            balanced[f"{key}_dye_swap_pp"] = 100.0 * (mcherry - gfp)
        balanced_rows.append(balanced)
    return orientation_rows, balanced_rows


# Mixture batches as reported by the experimenter (2026-10-16): DDX3X and DDX3X C1
# were mixed once per orientation and dispensed to every replicate of both arms;
# WT and WT C1 were mixed as two separate batches; SREBF1 and TADA2B were one mix
# each.  Replicates within a batch are therefore splits of one mixture at D0.
# Groups absent here have no recorded batch structure and are treated as one
# batch per group, flagged "(assumed)" in the report.
MIXTURE_BATCHES = {
    "WT": "WT",
    "WT_C1": "WT_C1",
    "DDX3X": "DDX3X",
    "DDX3X_C1": "DDX3X",
    "SREBF1": "SREBF1",
    "TADA2B": "TADA2B",
}


def build_baselines(wells: list[dict]) -> list[dict]:
    """Pool every D0 well of a mixture batch into one baseline per arm and orientation."""
    mixtures = [well for well in wells if well["role"] == "mixture"]
    rows = []
    groups = sorted({well["group"] for well in mixtures}, key=_group_order)
    for group in groups:
        batch = MIXTURE_BATCHES.get(group, group)
        for orientation in ORIENTATIONS:
            members = [
                well
                for well in mixtures
                if MIXTURE_BATCHES.get(well["group"], well["group"]) == batch
                and well["orientation"] == orientation
            ]
            loy = sum(int(well["loy_count"]) for well in members)
            total = sum(int(well["reporter_positive_singlets"]) for well in members)
            _mean, sd = _mean_sd([well["loy_fraction"] for well in members])
            rows.append(
                {
                    "group": group,
                    "orientation": orientation,
                    "mixture_batch": batch,
                    "batch_recorded": group in MIXTURE_BATCHES,
                    "wells": ";".join(well["well"] for well in members),
                    "loy_count": loy,
                    "reporter_positive_singlets": total,
                    "d0_loy_fraction": loy / total,
                    "split_sd_pp": 100.0 * sd,
                }
            )
    return rows


def _group_order(group: str) -> tuple[int, int]:
    background, _, condition = group.partition("_")
    return BACKGROUNDS.index(background), 1 if condition else 0


def _sha256(path: Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


SEPTEMBER_ROOT = Path("spectral") / "20260901_125230_Mike_Run1" / "Mike_Run1"
SEPTEMBER_UNSTAINED = ("A01-unstained_CTRL", "A02-unstained_LOY")


def _gated_raw_reporters(path: Path, gates: Mapping) -> np.ndarray:
    """Raw B3 (GFP) and YG4 (mCherry) detector values of the primary singlets.

    Uses the same Scatter -> SSC Singlets -> FSC Singlets polygons as
    ``_primary_features`` so rows align event-for-event with its output.
    """
    names = ("Scatter", "SSC Singlets", "FSC Singlets")
    required = {RAW_GFP, RAW_MCHERRY, *(p for n in names for p in gates[n].parameters)}
    events = load_fcs(path, required=sorted(required)).events
    keep = np.ones(len(events), dtype=bool)
    for name in names:
        gate = gates[name]
        keep &= points_in_polygon(events.loc[:, list(gate.parameters)].to_numpy(dtype=float), gate.vertices)
    return events.loc[keep, [RAW_GFP, RAW_MCHERRY]].to_numpy(dtype=float)


def unstained_reference(repo: Path) -> dict:
    """99th-percentile raw-detector signal of the 1-Sep-2026 unstained HT1080 wells."""
    root = repo / SEPTEMBER_ROOT
    files = _file_map(root / "1")
    gates = load_cef_gates(root / "20260901_125230_Mike_Run1.cef")
    raw = np.vstack([_gated_raw_reporters(files[name], gates) for name in SEPTEMBER_UNSTAINED])
    return {
        "source": [f"{(SEPTEMBER_ROOT / '1' / name).as_posix()}.fcs" for name in SEPTEMBER_UNSTAINED],
        "events": int(len(raw)),
        "percentile": UNSTAINED_PERCENTILE,
        "gfp_detector": RAW_GFP,
        "mcherry_detector": RAW_MCHERRY,
        "gfp_median": float(np.median(raw[:, 0])),
        "mcherry_median": float(np.median(raw[:, 1])),
        "gfp_threshold": float(np.percentile(raw[:, 0], UNSTAINED_PERCENTILE)),
        "mcherry_threshold": float(np.percentile(raw[:, 1], UNSTAINED_PERCENTILE)),
    }


def unstained_check(
    features: np.ndarray, raw: np.ndarray, classifier: LDAModel, reference: Mapping
) -> dict:
    """Compare the negative rule with the unstained mCherry-detector ceiling.

    ``ambiguous`` events are the broad-but-not-strict negatives (GFP-negative,
    unmixed mCherry in the noise band).  If they are reporter-negative they sit
    below the unstained ceiling on the raw YG4 detector; the mCherry-class
    events that the rule keeps should sit above it.
    """
    threshold = reference["mcherry_threshold"]
    ambiguous = negative_mask(features, "broad") & ~negative_mask(features, "strict")
    kept_mcherry = ~negative_mask(features, PRIMARY_NEGATIVE_RULE) & ~classifier.predict(features).astype(bool)
    below = raw[:, 1] <= threshold
    return {
        "ambiguous_fraction": float(ambiguous.mean()),
        "ambiguous_below_unstained_fraction": float(below[ambiguous].mean()) if ambiguous.any() else None,
        "kept_mcherry_above_unstained_fraction": (
            float((~below)[kept_mcherry].mean()) if kept_mcherry.any() else None
        ),
        "raw_mcherry_median_strict_negatives": (
            float(np.median(raw[negative_mask(features, "strict"), 1]))
            if negative_mask(features, "strict").any() else None
        ),
    }


def _september_classifier(repo: Path) -> tuple[LDAModel, Mapping[str, str]]:
    root = repo / SEPTEMBER_ROOT
    files = _file_map(root / "1")
    gates = load_cef_gates(root / "20260901_125230_Mike_Run1.cef")
    features, labels = [], []
    for name, label in {
        "A03-CTRL_GFP": 1,
        "A05-LOY_GFP": 1,
        "A07-CTRL_cherry": 0,
        "A09-LOY_cherry": 0,
    }.items():
        values, count, _saturated = _primary_features(files[name], gates)
        features.append(values)
        labels.append(np.full(count, label, dtype=int))
    model = fit_lda(np.vstack(features), np.concatenate(labels))
    return model, load_fcs(files["A03-CTRL_GFP"]).metadata


def analyze_export(repo_root: Path, export_root: Path) -> tuple[dict, dict[str, np.ndarray]]:
    """Gate, classify and summarise every Plate_004 recording."""
    repo = Path(repo_root).resolve()
    export = Path(export_root).resolve()
    cef_path = next(export.glob("*.cef"))
    statistics_path = next(export.glob("*_Statistics.csv"))
    plate = export / "Plate_004"
    paths = {path.stem: path for path in sorted(plate.glob("*.fcs"))}
    samples = build_design(paths)
    gates = load_cef_gates(cef_path)
    instrument_singlets = read_instrument_singlets(statistics_path)
    reference_classifier, reference_metadata = _september_classifier(repo)
    unstained = unstained_reference(repo)

    features: dict[str, np.ndarray] = {}
    raw_reporters: dict[str, np.ndarray] = {}
    saturated: dict[str, float] = {}
    for sample in samples:
        verify_unmixing_compatibility(reference_metadata, load_fcs(paths[sample.name]).metadata)
        values, count, saturated_fraction = _primary_features(paths[sample.name], gates)
        recording = sample.name.split("-", 1)[1]
        if instrument_singlets.get(recording) != count:
            raise ValueError(
                f"{sample.name}: {count} gated singlets but FACSChorus reports "
                f"{instrument_singlets.get(recording)}"
            )
        features[sample.name] = values
        raw_reporters[sample.name] = _gated_raw_reporters(paths[sample.name], gates)
        if len(raw_reporters[sample.name]) != count:
            raise ValueError(f"{sample.name}: raw-detector gating disagrees with feature gating")
        saturated[sample.name] = saturated_fraction

    # The classifiers are trained on the same reporter-positive events they are
    # applied to, so unlabelled cells in the SOLO wells do not pull the class means.
    solos = [sample for sample in samples if sample.role == "solo"]
    positive = {
        s.name: features[s.name][~negative_mask(features[s.name], PRIMARY_NEGATIVE_RULE)] for s in solos
    }

    def train(controls: list[GrowthSample], events: Mapping[str, np.ndarray]) -> LDAModel:
        return fit_lda(
            np.vstack([events[s.name] for s in controls]),
            np.concatenate([np.full(len(events[s.name]), int(s.dye == "GFP"), dtype=int) for s in controls]),
        )

    classifier = train(solos, positive)
    all_events_classifier = train(solos, features)
    per_background_models = {
        background: train([s for s in solos if s.background == background], positive)
        for background in BACKGROUNDS
    }

    wells = []
    for sample in samples:
        row = {
            "well": sample.well,
            "sample": sample.name,
            "role": sample.role,
            "background": sample.background,
            "condition": sample.condition,
            "group": sample.group,
            "orientation": sample.orientation,
            "genotype": sample.genotype,
            "reporter": sample.dye,
            "replicate": sample.replicate,
            "saturated_fraction": saturated[sample.name],
        }
        row.update(
            measure_well(features[sample.name], classifier, reference_classifier, sample, all_events_classifier)
        )
        row.update(unstained_check(features[sample.name], raw_reporters[sample.name], classifier, unstained))
        if sample.role == "mixture":
            kept = features[sample.name][~negative_mask(features[sample.name], PRIMARY_NEGATIVE_RULE)]
            gfp_own = float(per_background_models[sample.background].predict(kept).mean())
            row["loy_fraction_background_classifier"] = (
                gfp_own if sample.orientation == "LOY_GFP" else 1.0 - gfp_own
            )
        if sample.role == "solo":
            expected = 1.0 if sample.dye == "GFP" else 0.0
            row["classification_error"] = abs(row["gfp_fraction"] - expected)
        wells.append(row)

    orientation_rows, balanced_rows = summarise_conditions(wells)
    dye_neutral = []
    for genotype in ("LOY", "CTRL"):
        rows = [w for w in wells if w["role"] == "dye_neutral" and w["genotype"] == genotype]
        mean, sd = _mean_sd([w["gfp_fraction"] for w in rows])
        dye_neutral.append(
            {
                "genotype": genotype,
                "wells": ";".join(w["well"] for w in rows),
                "gfp_fraction_mean": mean,
                "gfp_fraction_sd": sd,
            }
        )

    manifest = []
    for path in [cef_path, statistics_path, *sorted(plate.iterdir())]:
        if path.is_file():
            manifest.append(
                {
                    "path": path.relative_to(repo).as_posix(),
                    "size_bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                }
            )

    result = {
        "title": "Spectral competition growth assay: D0 baseline, 15 October 2026",
        "analysis_date": "2026-10-16",
        "source": export.relative_to(repo).as_posix(),
        "instrument": "BD FACSDiscover A8",
        "unmixing_matrix": reference_metadata.get("bdspectral matrix name"),
        "unmixing_matches_september": True,
        "instrument_singlet_counts_reproduced": len(samples),
        "classifier": {
            "training": "reporter-positive events of all 16 SOLO controls (GFP = 1, mCherry = 0), equal priors",
            "weights": classifier.weights.tolist(),
            "all_events_weights": all_events_classifier.weights.tolist(),
            "all_events_intercept": all_events_classifier.intercept,
            "intercept": classifier.intercept,
            "september_weights": reference_classifier.weights.tolist(),
            "september_intercept": reference_classifier.intercept,
        },
        "negative_rules": {
            rule: {"gfp_below": gfp, "mcherry_below": mcherry}
            for rule, (gfp, mcherry) in NEGATIVE_RULES.items()
        },
        "primary_negative_rule": PRIMARY_NEGATIVE_RULE,
        "unstained_reference": unstained,
        "wells": wells,
        "orientation_summaries": orientation_rows,
        "balanced_summaries": balanced_rows,
        "d0_baselines": build_baselines(wells),
        "dye_neutral": dye_neutral,
        "source_manifest": manifest,
    }
    return result, features


# ---------------------------------------------------------------- figures --


def _svg_open(width: int, height: int, title: str, description: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'font-family="sans-serif" role="img">',
        f"<title>{escape(title)}</title>",
        f"<desc>{escape(description)}</desc>",
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="20" y="28" font-size="16" fill="#111827">{escape(title)}</text>',
    ]


def _write(parts: list[str], output: Path) -> None:
    parts.append("</svg>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(parts), encoding="utf-8")


def write_well_sheet_svg(
    wells: list[dict],
    features: Mapping[str, np.ndarray],
    classifier: Mapping,
    title: str,
    output: Path,
    columns: int = 6,
    points_per_well: int = 1500,
) -> None:
    """Reporter-space panels with the LDA boundary and both negative rules."""
    x_limits, y_limits = (-8.0, 12.5), (-8.0, 12.5)
    cell_w, cell_h, plot = 230, 250, 175
    rows = max(1, math.ceil(len(wells) / columns))
    width, height = 20 + columns * cell_w, 80 + rows * cell_h
    parts = _svg_open(
        width,
        height,
        title,
        "Each panel: primary singlets on arcsinh(x/1000) reporter axes. Points are a deterministic "
        "subsample coloured by the shared LDA class; dashed line = LDA boundary; grey rectangles = "
        "strict (dark) and broad (light) reporter-negative rules. Events inside the broad rectangle are "
        "excluded; LOY / GFP-class percentages are of the remaining reporter-positive events.",
    )
    parts.append(
        f'<circle cx="24" cy="50" r="4" fill="{GFP_COLOUR}"/><text x="33" y="54" font-size="12" fill="#374151">GFP class</text>'
        f'<circle cx="114" cy="50" r="4" fill="{MCHERRY_COLOUR}"/><text x="123" y="54" font-size="12" fill="#374151">mCherry class</text>'
        f'<line x1="224" y1="50" x2="252" y2="50" stroke="#111827" stroke-width="1.5" stroke-dasharray="5,3"/>'
        f'<text x="258" y="54" font-size="12" fill="#374151">LDA boundary</text>'
        f'<rect x="352" y="44" width="14" height="12" fill="{MUTED}" fill-opacity="0.12" stroke="{MUTED}" stroke-dasharray="3,2"/>'
        f'<text x="372" y="54" font-size="12" fill="#374151">reporter-negative rules (strict / broad = excluded)</text>'
    )
    boundary = _boundary_segment(
        classifier["weights"], classifier["intercept"], x_limits, y_limits
    )
    model = LDAModel(np.asarray(classifier["weights"]), float(classifier["intercept"]))
    for index, well in enumerate(wells):
        cx = 20 + (index % columns) * cell_w
        cy = 70 + (index // columns) * cell_h
        left, top = cx + 38, cy + 36

        def px(value: float) -> float:
            return left + (value - x_limits[0]) / (x_limits[1] - x_limits[0]) * plot

        def py(value: float) -> float:
            return top + plot - (value - y_limits[0]) / (y_limits[1] - y_limits[0]) * plot

        name = well["sample"]
        if well["role"] == "mixture":
            stat = f"LOY {100 * well['loy_fraction']:.1f}% | excl. {100 * well['negative_broad_fraction']:.1f}%"
        elif well["role"] == "solo":
            stat = f"GFP class {100 * well['gfp_fraction']:.1f}% | excl. {100 * well['negative_broad_fraction']:.1f}%"
        else:
            stat = f"GFP class {100 * well['gfp_fraction']:.1f}%"
        parts.append(f'<g data-sample="{escape(name)}">')
        parts.append(f'<text x="{cx + 8}" y="{cy + 12}" font-size="11" fill="#111827">{escape(name)}</text>')
        parts.append(
            f'<text x="{cx + 8}" y="{cy + 27}" font-size="10.5" fill="#4b5563">n={well["primary_singlets"]:,} | {stat}</text>'
        )
        parts.append(
            f'<rect x="{left}" y="{top}" width="{plot}" height="{plot}" fill="#fafafa" stroke="#d1d5db"/>'
        )
        for rule, opacity in (("broad", 0.08), ("strict", 0.16)):
            gfp_limit, mcherry_limit = NEGATIVE_RULES[rule]
            parts.append(
                f'<rect x="{left}" y="{py(mcherry_limit):.1f}" width="{px(gfp_limit) - left:.1f}" '
                f'height="{top + plot - py(mcherry_limit):.1f}" fill="{MUTED}" fill-opacity="{opacity}" '
                f'stroke="{MUTED}" stroke-width="0.8" stroke-dasharray="3,2"/>'
            )
        points = np.clip(_sample_features(features[name], points_per_well), -7.9, 12.4)
        classes = model.predict(points)
        for (x, y), cls in zip(points, classes):
            colour = GFP_COLOUR if cls else MCHERRY_COLOUR
            parts.append(f'<circle cx="{px(x):.1f}" cy="{py(y):.1f}" r="1.1" fill="{colour}" fill-opacity="0.35"/>')
        if boundary is not None:
            (x1, y1), (x2, y2) = boundary
            parts.append(
                f'<line x1="{px(x1):.1f}" y1="{py(y1):.1f}" x2="{px(x2):.1f}" y2="{py(y2):.1f}" '
                f'stroke="#111827" stroke-width="1.3" stroke-dasharray="5,3"/>'
            )
        for tick in (-5, 0, 5, 10):
            parts.append(f'<text x="{px(tick):.1f}" y="{top + plot + 12}" font-size="9" text-anchor="middle" fill="#6b7280">{tick}</text>')
            parts.append(f'<text x="{left - 4}" y="{py(tick) + 3:.1f}" font-size="9" text-anchor="end" fill="#6b7280">{tick}</text>')
        parts.append(f'<text x="{left + plot / 2}" y="{top + plot + 25}" font-size="10" text-anchor="middle" fill="#4b5563">eGFP</text>')
        parts.append(
            f'<text x="{left - 24}" y="{top + plot / 2}" font-size="10" text-anchor="middle" fill="#4b5563" '
            f'transform="rotate(-90 {left - 24} {top + plot / 2})">mCherry</text>'
        )
        parts.append("</g>")
    _write(parts, output)


def write_d0_summary_svg(result: Mapping, output: Path) -> None:
    """Per-well D0 LOY fraction by condition and dye orientation."""
    groups = [row["group"] for row in result["balanced_summaries"]]
    width = 860
    row_h, top, left, right = 46, 92, 150, 700
    height = top + row_h * len(groups) + 160
    parts = _svg_open(
        width,
        height,
        "D0 LOY fraction by knockout background and dye orientation",
        "Small marks are replicate wells (circle = LOY carries GFP, triangle = LOY carries mCherry); "
        "the black tick is the dye-balanced (logit-mean) estimate. Lower panel: GFP fraction in "
        "same-genotype GFP:mCherry mixtures.",
    )
    x = lambda value: left + value * (right - left)
    parts.append(
        f'<circle cx="24" cy="54" r="5" fill="{GFP_COLOUR}"/><text x="34" y="58" font-size="12" fill="#374151">LOY = GFP, CTRL = mCherry</text>'
        f'<polygon points="214,48 208,59 220,59" fill="{MCHERRY_COLOUR}"/><text x="228" y="58" font-size="12" fill="#374151">LOY = mCherry, CTRL = GFP</text>'
        f'<line x1="430" y1="47" x2="430" y2="61" stroke="#111827" stroke-width="2.5"/><text x="438" y="58" font-size="12" fill="#374151">dye-balanced estimate</text>'
    )
    bottom = top + row_h * len(groups)
    for tick in range(0, 101, 10):
        xpos = x(tick / 100)
        parts.append(f'<line x1="{xpos:.1f}" y1="{top - 8}" x2="{xpos:.1f}" y2="{bottom}" stroke="{"#9ca3af" if tick == 50 else "#eef0f2"}"/>')
        parts.append(f'<text x="{xpos:.1f}" y="{bottom + 16}" font-size="11" text-anchor="middle" fill="#6b7280">{tick}%</text>')
    parts.append(f'<text x="{(left + right) / 2}" y="{bottom + 34}" font-size="12" text-anchor="middle" fill="#374151">LOY fraction of reporter-positive singlets (two-class LDA, broad negatives excluded)</text>')
    balanced = {row["group"]: row for row in result["balanced_summaries"]}
    for index, group in enumerate(groups):
        ycentre = top + row_h * index + row_h / 2
        parts.append(f'<text x="{left - 14}" y="{ycentre + 4:.1f}" font-size="13" text-anchor="end" fill="#111827">{escape(group.replace("_", " "))}</text>')
        for orientation, offset in (("LOY_GFP", -7), ("LOY_mCherry", 7)):
            for well in result["wells"]:
                if well["role"] != "mixture" or well["group"] != group or well["orientation"] != orientation:
                    continue
                xpos, ypos = x(well["loy_fraction"]), ycentre + offset
                tip = f"{well['sample']}: LOY {100 * well['loy_fraction']:.1f}%"
                if orientation == "LOY_GFP":
                    parts.append(f'<circle cx="{xpos:.1f}" cy="{ypos:.1f}" r="4.5" fill="{GFP_COLOUR}" fill-opacity="0.85" stroke="white" stroke-width="1"><title>{escape(tip)}</title></circle>')
                else:
                    parts.append(f'<polygon points="{xpos:.1f},{ypos - 5.5:.1f} {xpos - 5:.1f},{ypos + 4:.1f} {xpos + 5:.1f},{ypos + 4:.1f}" fill="{MCHERRY_COLOUR}" fill-opacity="0.85" stroke="white" stroke-width="1"><title>{escape(tip)}</title></polygon>')
        value = balanced[group]["loy_fraction_balanced"]
        swap = balanced[group]["loy_fraction_dye_swap_pp"]
        parts.append(f'<line x1="{x(value):.1f}" y1="{ycentre - 15:.1f}" x2="{x(value):.1f}" y2="{ycentre + 15:.1f}" stroke="#111827" stroke-width="2.5"><title>{escape(group)} balanced {100 * value:.1f}%</title></line>')
        parts.append(f'<text x="{right + 14}" y="{ycentre + 4:.1f}" font-size="11.5" fill="#374151">{100 * value:.1f}%  (swap {swap:+.1f} pp)</text>')
    # Dye-neutral panel.
    panel_top = bottom + 64
    parts.append(f'<text x="20" y="{panel_top - 10}" font-size="13" fill="#111827">Same-genotype dye-neutral mixtures (GFP fraction; 50% = no dye bias in mixing/measurement)</text>')
    for index, genotype in enumerate(("LOY", "CTRL")):
        ycentre = panel_top + 22 + index * 34
        parts.append(f'<line x1="{left}" y1="{ycentre:.1f}" x2="{right}" y2="{ycentre:.1f}" stroke="#eef0f2"/>')
        parts.append(f'<line x1="{x(0.5):.1f}" y1="{ycentre - 12:.1f}" x2="{x(0.5):.1f}" y2="{ycentre + 12:.1f}" stroke="#9ca3af"/>')
        parts.append(f'<text x="{left - 14}" y="{ycentre + 4:.1f}" font-size="13" text-anchor="end" fill="#111827">{genotype} GFP:mCherry</text>')
        values = [w["gfp_fraction"] for w in result["wells"] if w["role"] == "dye_neutral" and w["genotype"] == genotype]
        for value in values:
            parts.append(f'<circle cx="{x(value):.1f}" cy="{ycentre:.1f}" r="4.5" fill="{MUTED}" fill-opacity="0.85" stroke="white" stroke-width="1"><title>GFP {100 * value:.1f}%</title></circle>')
        parts.append(f'<text x="{right + 14}" y="{ycentre + 4:.1f}" font-size="11.5" fill="#374151">GFP {100 * float(np.mean(values)):.1f}%</text>')
    _write(parts, output)


# ---------------------------------------------------------------- report --


def _pct(value: float) -> str:
    return f"{100 * value:.1f}%"


def _c1_gap(orientation: Mapping, background: str, dye_orientation: str) -> float:
    """Percentage-point difference of the C1 (inhibitor) arm from its untreated partner at D0."""
    treated = orientation[(f"{background}_C1", dye_orientation)]["loy_fraction_mean"]
    untreated = orientation[(background, dye_orientation)]["loy_fraction_mean"]
    return 100.0 * (treated - untreated)


def _range(values: Iterable[float], digits: int = 1) -> str:
    values = list(values)
    return f"{min(values):.{digits}f}–{max(values):.{digits}f}"


def build_markdown(result: Mapping) -> str:
    solos = [w for w in result["wells"] if w["role"] == "solo"]
    max_error = max(w["classification_error"] for w in solos)
    worst = max(solos, key=lambda w: w["classification_error"])
    max_all_error = max(
        abs(w["gfp_fraction_all_events"] - (1.0 if w["reporter"] == "GFP" else 0.0)) for w in solos
    )
    orientation = {(r["group"], r["orientation"]): r for r in result["orientation_summaries"]}
    max_sd = max(
        max(r["loy_fraction_sd"] for r in result["orientation_summaries"]), 0.0
    )
    balanced = {row["group"]: row for row in result["balanced_summaries"]}
    swaps = [balanced[g]["loy_fraction_dye_swap_pp"] for g in ("WT", "DDX3X", "DDX3X_C1", "SREBF1")]
    unstained = result["unstained_reference"]
    checked = [w for w in result["wells"] if w["ambiguous_fraction"] >= 0.005]
    kept_mcherry = [
        w["kept_mcherry_above_unstained_fraction"]
        for w in result["wells"]
        if w["role"] == "mixture" or (w["role"] == "solo" and w["reporter"] == "mCherry")
    ]

    def strict_yg4(selector) -> str:
        return _range(
            (w["raw_mcherry_median_strict_negatives"] / 1000 for w in result["wells"] if selector(w)), 0
        )

    def shift(row: Mapping, key: str) -> float:
        return abs(row[f"{key}_mean"] - row["loy_fraction_mean"])

    def own(row: Mapping) -> float:
        wells = [
            w for w in result["wells"]
            if w["role"] == "mixture" and w["group"] == row["group"] and w["orientation"] == row["orientation"]
        ]
        return float(np.mean([w["loy_fraction_background_classifier"] for w in wells]))

    rows = result["orientation_summaries"]
    classifier_shift = max(
        max(shift(r, "loy_fraction_september_classifier"), abs(own(r) - r["loy_fraction_mean"])) for r in rows
    )
    all_events_shift = {
        background: max(shift(r, "loy_fraction_all_events") for r in rows if r["background"] == background)
        for background in BACKGROUNDS
    }
    lines = [
        f"# {result['title']}",
        "",
        f"**Source:** `{result['source']}` (Plate_004, 58 recordings, 10,000 events each)  ",
        f"**Instrument:** {result['instrument']}; unmixing `{result['unmixing_matrix']}` "
        "(identical spillover payload to the August/September runs)  ",
        f"**Analysed:** {result['analysis_date']}",
        "",
        "## Summary",
        "",
        "This is a baseline-only acquisition, so it says nothing yet about competitive fitness; "
        "it fixes the D0 LOY:CTRL composition that later timepoints will be compared against.",
        "",
        "LOY fraction of reporter-positive singlets (broad reporter-negative events excluded).",
        "",
        "| Condition | LOY = GFP (mean ± SD) | LOY = mCherry (mean ± SD) | Dye-balanced | Dye swap (mCh − GFP) | Excluded as reporter-negative |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in result["balanced_summaries"]:
        group = row["group"]
        gfp, mch = orientation[(group, "LOY_GFP")], orientation[(group, "LOY_mCherry")]
        neg = 0.5 * (gfp["negative_broad_fraction_mean"] + mch["negative_broad_fraction_mean"])
        lines.append(
            f"| {group.replace('_', ' ')} | {_pct(gfp['loy_fraction_mean'])} ± {100 * gfp['loy_fraction_sd']:.1f} "
            f"| {_pct(mch['loy_fraction_mean'])} ± {100 * mch['loy_fraction_sd']:.1f} "
            f"| {_pct(row['loy_fraction_balanced'])} | {row['loy_fraction_dye_swap_pp']:+.1f} pp | {_pct(neg)} |"
        )
    neutral = {row["genotype"]: row for row in result["dye_neutral"]}
    lines += [
        "",
        f"Same-genotype dye-neutral mixtures: LOY GFP:mCherry = {_pct(neutral['LOY']['gfp_fraction_mean'])} GFP "
        f"(SD {100 * neutral['LOY']['gfp_fraction_sd']:.1f} pp); CTRL GFP:mCherry = "
        f"{_pct(neutral['CTRL']['gfp_fraction_mean'])} GFP (SD {100 * neutral['CTRL']['gfp_fraction_sd']:.1f} pp).",
        "",
        "Key points:",
        "",
        f"- Replicate wells agree closely (largest within-orientation SD {100 * max_sd:.1f} pp). In every "
        "condition the replicates are splits of one mixture, so this SD measures dispensing and "
        "acquisition noise, not variation between independently prepared mixtures.",
        f"- The two dye orientations disagree by {_range(swaps)} pp in WT, DDX3X, DDX3X C1 and SREBF1 (LOY "
        f"is higher when it carries mCherry) and by {balanced['WT_C1']['loy_fraction_dye_swap_pp']:+.1f} pp in "
        "WT C1. As in the earlier runs this is a starting-composition offset. The growth readout (change "
        "in log-odds from each orientation's own D0) removes it, provided each orientation is followed "
        "separately.",
        "- **Reporter-negative events are excluded before the ratio is taken** (the broad rule; "
        "see *Unstained check*). Counting them, as the August/September reports did, puts them in "
        "the mCherry class and inflates the dye swap, most of all in TADA2B "
        f"({balanced['TADA2B']['loy_fraction_all_events_dye_swap_pp']:+.1f} pp with every event kept, "
        f"{balanced['TADA2B']['loy_fraction_dye_swap_pp']:+.1f} pp after exclusion).",
        "- **C1 (DDX3X-inhibitor) arms are untreated at D0**, so they should match their untreated "
        "partners. DDX3X C1 does: the orientation means differ by "
        f"{_c1_gap(orientation, 'DDX3X', 'LOY_GFP'):+.1f} and {_c1_gap(orientation, 'DDX3X', 'LOY_mCherry'):+.1f} pp. "
        "WT C1 does not: it differs by "
        f"{_c1_gap(orientation, 'WT', 'LOY_GFP'):+.1f} pp (LOY = GFP) and "
        f"{_c1_gap(orientation, 'WT', 'LOY_mCherry'):+.1f} pp (LOY = mCherry), far outside the replicate "
        "SD. This matches the mixing record: WT and WT C1 were two separate batches, whereas DDX3X and "
        "DDX3X C1 came from one batch. Any inhibitor "
        "effect has to be measured as each arm's change from its own D0, not as a comparison of "
        "endpoints between arms.",
        "- **TADA2B is the outlier.** About a fifth of its events are reporter-negative, in mixtures and "
        "SOLO controls alike. Unlike the other backgrounds, its double-negative events have the raw "
        "mCherry-detector signal of intact unstained HT1080 (see *Unstained check*), so they look like "
        "live cells that have lost or silenced the reporter, not debris. Its dye swap is still "
        f"{balanced['TADA2B']['loy_fraction_dye_swap_pp']:+.1f} pp after exclusion. If reporter loss "
        "continues during culture, TADA2B growth estimates will be biased, so read them with that caveat.",
        "",
        "![D0 LOY fraction](figures/d0_loy_fraction.svg)",
        "",
        "## Design as parsed from the well names",
        "",
        "- `LOY_<dye>_<KO>_<rep>`: LOY cells carry `<dye>`, CTRL cells carry the other reporter, both in "
        "the `<KO>` background (WT = no knockout). LOY fraction = fraction of events in LOY's reporter class.",
        "- `..._C1_<rep>`: the WT and DDX3X plates that will receive the C1 DDX3X inhibitor. At D0 they "
        "are untreated WT or DDX3X HT1080, but they are kept as separate arms (never pooled with the "
        "standard wells) because they are followed as separate cultures.",
        "- `<genotype>_GFP_CHERRY_<rep>`: same-genotype GFP:mCherry mixtures (dye-neutral controls).",
        "- `<genotype>_<dye>_SOLO_<KO>`: single-colour controls, four per background.",
        "- Mixing batches, as recorded by the experimenter: DDX3X and DDX3X C1 were mixed once per dye "
        "orientation and dispensed to all replicates of both arms. WT and WT C1 were mixed as two "
        "separate batches, one per arm. SREBF1 and TADA2B were each one mix per orientation.",
        "- Replicates `_1`–`_3` are therefore splits of a shared D0 mixture, not independent mixtures. "
        "Trailing `_001`–`_006` are FACSChorus re-record counters and are ignored.",
        "",
        "## D0 baselines for the growth readout",
        "",
        "Each arm's D0 is the event-pooled LOY fraction of every well from its mixing batch: the summed "
        "LOY-class events divided by the summed reporter-positive singlets. This follows the convention of the "
        "Run 1 report. DDX3X and DDX3X C1 share one baseline (six wells per orientation); WT and WT C1 "
        "each have their own. At later timepoints each replicate plate is an independently cultured "
        "unit. Its change in log-odds is taken against this shared baseline, so the interval will "
        "reflect plate-to-plate culture variation, not mixture-preparation variation, as in Run 1.",
        "",
        "| Arm | Orientation | Mixing batch | Wells pooled | Reporter-positive singlets | D0 LOY | Split SD |",
        "|---|---|---|---|---:|---:|---:|",
        *[
            f"| {row['group'].replace('_', ' ')} | {row['orientation']} | "
            f"{row['mixture_batch'].replace('_', ' ')}{'' if row['batch_recorded'] else ' (assumed)'} | "
            f"{row['wells'].replace(';', ', ')} | {row['reporter_positive_singlets']:,} | "
            f"{_pct(row['d0_loy_fraction'])} | {row['split_sd_pp']:.1f} pp |"
            for row in result["d0_baselines"]
        ],
        "",
        "Comparing the inhibitor arm with its untreated arm is cleaner for DDX3X than for WT. The DDX3X "
        "arms share a mixture, so only the treatment differs. For WT, the treatment is confounded with "
        "one mixing batch per arm. Within-arm changes from D0 remove the ~5 pp starting offset, but any "
        "batch effect on subsequent growth cannot be separated from the inhibitor effect.",
        "",
        "## Methods",
        "",
        "1. Raw FCS files were read without modification. The run's own CEF polygons "
        "(Scatter → SSC Singlets → FSC Singlets) were applied in raw units. The gated singlet count "
        f"reproduces the FACSChorus statistics export exactly in all {result['instrument_singlet_counts_reproduced']} recordings.",
        "2. Unmixed `eGFP-A` and `mCherry-A` were transformed as `arcsinh(value / 1000)`.",
        "3. Reporter-negative events were removed with the broad rule (eGFP < "
        f"{NEGATIVE_RULES['broad'][0]:g} and mCherry < {NEGATIVE_RULES['broad'][1]:g} in transformed units). "
        "This removes the double-negative cluster and the GFP-negative events whose unmixed mCherry sits "
        "in the unmixing-noise band (about +5). The rule was checked against unstained cells from the "
        "1 September run (step 4). Every later timepoint of this assay must use the same rule.",
        "4. Unstained check. The raw `YG4 (625)-A` detector is upstream of unmixing, so its scale is "
        "comparable between days with the same YG settings. The September unstained CTRL and LOY wells "
        f"({unstained['events']:,} primary singlets) give a {unstained['percentile']:g}th-percentile ceiling "
        f"of {unstained['mcherry_threshold']:,.0f} (median {unstained['mcherry_median']:,.0f}).",
        "5. An equal-prior two-class LDA (GFP vs mCherry) was trained on the reporter-positive events of "
        "all 16 SOLO controls from this plate and applied to the reporter-positive events of every well. "
        "LOY fraction = LOY-class events / reporter-positive events.",
        "6. Sensitivity analyses: (a) every event kept, with an LDA trained on every SOLO event (the "
        "August/September method); (b) only the strict double-negative cluster excluded (eGFP < "
        f"{NEGATIVE_RULES['strict'][0]:g} and mCherry < {NEGATIVE_RULES['strict'][1]:g}); (c) the September "
        "classifier; (d) a separate LDA per knockout background.",
        "7. Dye-balanced D0 = inverse-logit of the mean of the two orientation logits.",
        "8. Saturated events are retained, as before (the flag is dominated by the ImgB3 imaging detector).",
        "",
        "## Unstained check",
        "",
        "No unstained control was on this plate. On the unmixed axes the September unstained cells cannot "
        "be transferred: this plate's reporter-positive cells are about one arcsinh unit brighter, the "
        "violet detector settings are about 10% higher, and the double-negative cluster moves between "
        "wells. On the raw 625 nm detector the comparison does work:",
        "",
        "- Events removed by the broad rule but not the strict rule (GFP-negative, unmixed mCherry about +5) "
        f"are below the unstained ceiling in {_range((100 * w['ambiguous_below_unstained_fraction'] for w in checked), 0)}% "
        f"of events, across all {len(checked)} wells where they make up at least 0.5% of events. They carry "
        "no mCherry and are correctly excluded.",
        "- mCherry-class events that the rule keeps are above the ceiling in "
        f"{_range((100 * v for v in kept_mcherry), 0)}% of events (mixtures and mCherry SOLO controls).",
        "- The double-negative cluster differs by background. In TADA2B mixtures its median YG4 signal is "
        f"{strict_yg4(lambda w: w['role'] == 'mixture' and w['background'] == 'TADA2B')}k, the same as unstained HT1080 "
        f"({unstained['mcherry_median'] / 1000:.0f}k), so these are intact, unlabelled cells. In WT and DDX3X "
        f"mixtures it is {strict_yg4(lambda w: w['role'] == 'mixture' and w['background'] in {'WT', 'DDX3X'})}k, "
        "well below an unstained cell, which is more like debris or dead cells than reporter loss.",
        "- Limits: the September plate had no knockout backgrounds, so the ceiling is the autofluorescence "
        "of unlabelled CTRL and LOY HT1080. Raw signal is only comparable while the YG detector settings stay the same. "
        "A ceiling-only rule would also discard the "
        f"{_range((100 - 100 * v for v in kept_mcherry), 0)}% of mCherry-class cells that sit below it "
        "(and 18% of the September mCherry control), so the ceiling is used to check the rule, not to define it.",
        "",
        "## Controls",
        "",
        f"After exclusion, the largest SOLO misclassification is {_pct(max_error)} ({worst['sample']}). "
        f"With every event kept it reached {_pct(max_all_error)}, because unlabelled cells are forced into "
        "the mCherry class.",
        "",
        "| Control | Primary singlets | Excluded (broad) | Strict negatives | GFP class (reporter-positive) | Error | GFP class, all events |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for w in solos:
        lines.append(
            f"| {w['sample']} | {w['primary_singlets']:,} | {_pct(w['negative_broad_fraction'])} "
            f"| {_pct(w['negative_strict_fraction'])} | {_pct(w['gfp_fraction'])} | {_pct(w['classification_error'])} "
            f"| {_pct(w['gfp_fraction_all_events'])} |"
        )
    lines += [
        "",
        "Notable controls: the DDX3X LOY-GFP control has a clear unlabelled cluster, the SREBF1 "
        "CTRL-mCherry control has about a quarter of its events excluded, and all four TADA2B controls "
        "except CTRL-GFP have 16–28% excluded. In TADA2B these are live unlabelled cells (see *Unstained "
        "check*), so the TADA2B reporter lines themselves are not fully labelled.",
        "",
        "## Sensitivity of D0 LOY fraction",
        "",
        "Orientation means (%).",
        "",
        "| Condition | Orientation | Primary (excl. broad neg.) | All events (Aug/Sep method) | Excl. strict neg. only | September LDA | Per-background LDA |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in result["orientation_summaries"]:
        wells = [
            w for w in result["wells"]
            if w["role"] == "mixture" and w["group"] == row["group"] and w["orientation"] == row["orientation"]
        ]
        lines.append(
            f"| {row['group'].replace('_', ' ')} | {row['orientation']} | {100 * row['loy_fraction_mean']:.1f} "
            f"| {100 * row['loy_fraction_all_events_mean']:.1f} "
            f"| {100 * row['loy_fraction_excluding_strict_negatives_mean']:.1f} "
            f"| {100 * row['loy_fraction_september_classifier_mean']:.1f} | {100 * own(row):.1f} |"
        )
    lines += [
        "",
        f"Once negatives are excluded, classifier choice moves orientation means by at most {100 * classifier_shift:.1f} pp. "
        "Keeping every event instead shifts them by up to "
        + ", ".join(f"{100 * all_events_shift[b]:.1f} pp ({b})" for b in BACKGROUNDS)
        + ". The values here are therefore not directly comparable with the August/September reports, which "
        "kept every event. Within this assay that does not matter, provided every timepoint uses the same rule.",
        "",
        "## Per-well reporter plots",
        "",
        "- [SOLO controls](figures/wells_controls.svg)",
        "- [WT, WT C1 (inhibitor arm) and dye-neutral mixtures](figures/wells_WT.svg)",
        "- [DDX3X and DDX3X C1 (inhibitor arm) mixtures](figures/wells_DDX3X.svg)",
        "- [SREBF1 and TADA2B mixtures](figures/wells_SREBF1_TADA2B.svg)",
        "",
        "## Recommendations for the next timepoints",
        "",
        "- Acquire an unstained control for each knockout background on each day: the knockout lines "
        "without reporter, or at least parental HT1080. That lets reporter-negative cells be gated on the "
        "day, rather than borrowed from September.",
        "- Keep the broad negative rule, and keep this plate's SOLO-trained classifier (or retrain on each "
        "day's SOLO controls and report both). Follow each dye orientation separately from its own D0.",
        "- For future inhibitor experiments, mix the treated and untreated arms of each background as "
        "one batch (as was done for DDX3X), so that treatment is not confounded with mixing.",
        "- Prepare at least two independent mixtures per orientation if the uncertainty should cover "
        "mixture preparation rather than only post-split culture variability.",
        "- Re-sort or re-check the TADA2B reporter lines. A fifth of their cells are unlabelled live cells, "
        "and further reporter loss during culture would bias the growth estimate.",
        "",
        "## Files",
        "",
        "`wells.csv` (every well, every rule, unstained check) · `conditions.csv` (orientation means/SD) · "
        "`balanced.csv` (dye-balanced D0 and dye-swap) · `d0_baselines.csv` (event-pooled D0 per arm "
        "and orientation, for the growth readout) · `source_manifest.csv` (size + SHA-256 of all "
        f"{len(result['source_manifest'])} source files) · `results.json` (everything above, machine-readable).",
        "",
        "## Reproduction",
        "",
        "```powershell",
        ".\\codex\\flow_env\\python.exe -m spectral.a8_pipeline.growth_assay_d0_report",
        "```",
        "",
        "Refuses to write into an output directory that already contains files.",
        "",
    ]
    return "\n".join(lines)


def _write_csv(path: Path, rows: list[dict]) -> None:
    fields: list[str] = []
    for row in rows:
        fields.extend(key for key in row if key not in fields)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(result: Mapping, features: Mapping[str, np.ndarray], output: Path) -> None:
    destination = Path(output)
    # Empty directories are tolerated because Dropbox can hold them open after a delete.
    if destination.exists() and any(path.is_file() for path in destination.rglob("*")):
        raise FileExistsError(f"output already contains files: {destination}")
    figures = destination / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    wells = list(result["wells"])
    classifier = result["classifier"]
    write_d0_summary_svg(result, figures / "d0_loy_fraction.svg")
    write_well_sheet_svg(
        [w for w in wells if w["role"] == "solo"], features, classifier,
        "SOLO controls (4 per knockout background)", figures / "wells_controls.svg", columns=4,
    )
    write_well_sheet_svg(
        [w for w in wells if w["background"] == "WT" and w["role"] != "solo"], features, classifier,
        "WT, WT C1 (DDX3X-inhibitor arm, untreated at D0) and dye-neutral mixtures", figures / "wells_WT.svg",
    )
    write_well_sheet_svg(
        [w for w in wells if w["background"] == "DDX3X" and w["role"] == "mixture"], features, classifier,
        "DDX3X and DDX3X C1 (DDX3X-inhibitor arm, untreated at D0) mixtures", figures / "wells_DDX3X.svg",
    )
    write_well_sheet_svg(
        [w for w in wells if w["background"] in {"SREBF1", "TADA2B"} and w["role"] == "mixture"],
        features, classifier, "SREBF1 and TADA2B mixtures", figures / "wells_SREBF1_TADA2B.svg",
    )
    (destination / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    (destination / "report.md").write_text(build_markdown(result), encoding="utf-8")
    _write_csv(destination / "wells.csv", wells)
    _write_csv(destination / "conditions.csv", list(result["orientation_summaries"]))
    _write_csv(destination / "balanced.csv", list(result["balanced_summaries"]))
    _write_csv(destination / "d0_baselines.csv", list(result["d0_baselines"]))
    _write_csv(destination / "source_manifest.csv", list(result["source_manifest"]))


def main(argv: list[str] | None = None) -> int:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Generate the 15-Oct-2026 D0 growth-assay report")
    parser.add_argument(
        "--input",
        type=Path,
        default=repo / "20261016_171604_Export" / "comp_growth_assay_D0_OCT_15",
    )
    parser.add_argument(
        "--output", type=Path, default=repo / "analysis" / "spectral_growth_assay_D0_20261015"
    )
    args = parser.parse_args(argv)
    result, features = analyze_export(repo, args.input)
    write_outputs(result, features, args.output.resolve())
    print(args.output.resolve() / "report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
