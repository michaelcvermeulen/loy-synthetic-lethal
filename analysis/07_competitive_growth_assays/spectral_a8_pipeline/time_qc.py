from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


QC_SIGNALS = (
    "FSC-A",
    "SSC (Violet)-A",
    "eGFP-A",
    "mCherry-A",
    "Autofluorescence_001-A",
)


@dataclass(frozen=True)
class Segment:
    start: int
    stop: int
    label: str

    @property
    def size(self) -> int:
        return self.stop - self.start


@dataclass(frozen=True)
class TimeQCResult:
    good_events: np.ndarray
    bins: pd.DataFrame


def detect_segments(time: np.ndarray) -> list[Segment]:
    values = np.asarray(time, dtype=float)
    if values.ndim != 1:
        raise ValueError("time must be one-dimensional")
    if values.size == 0:
        return []
    delta = np.diff(values)
    positive = delta[delta > 0]
    typical = float(np.median(positive)) if positive.size else 0.0
    boundaries = np.flatnonzero((delta < 0) | ((delta > 1.0) & (delta > 20.0 * typical))) + 1
    edges = np.r_[0, boundaries, values.size]
    return [
        Segment(int(start), int(stop), f"Run {index + 1}")
        for index, (start, stop) in enumerate(zip(edges[:-1], edges[1:]))
    ]


def robust_outliers(values: np.ndarray, z_limit: float = 3.5) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    median = float(np.median(array))
    mad = float(1.4826 * np.median(np.abs(array - median)))
    if mad == 0.0:
        return ~np.isclose(array, median, rtol=1e-12, atol=1e-12)
    return np.abs(array - median) / mad > z_limit


def time_qc(
    events: pd.DataFrame,
    segments: list[Segment],
    min_events: int = 200,
    z_limit: float = 3.5,
) -> TimeQCResult:
    required = {"Time", *QC_SIGNALS}
    missing = sorted(required - set(events.columns))
    if missing:
        raise KeyError(f"missing time-QC parameters: {missing}")
    if min_events < 2:
        raise ValueError("min_events must be at least 2")

    rows: list[dict[str, float | int | str | bool]] = []
    event_bins: list[np.ndarray] = []
    for segment in segments:
        count = max(1, segment.size // min_events)
        for local_bin, indices in enumerate(np.array_split(np.arange(segment.start, segment.stop), count)):
            time_values = events.iloc[indices]["Time"].to_numpy(dtype=float)
            duration = float(time_values[-1] - time_values[0]) if indices.size > 1 else 0.0
            rate = float((indices.size - 1) / duration) if duration > 0 else float("inf")
            row: dict[str, float | int | str | bool] = {
                "segment": segment.label,
                "bin": local_bin + 1,
                "start": int(indices[0]),
                "stop": int(indices[-1] + 1),
                "events": int(indices.size),
                "time_start": float(time_values[0]),
                "time_stop": float(time_values[-1]),
                "event_rate": rate,
            }
            for signal in QC_SIGNALS:
                row[signal] = float(np.median(events.iloc[indices][signal].to_numpy(dtype=float)))
            rows.append(row)
            event_bins.append(indices)

    table = pd.DataFrame(rows)
    table["bad"] = False
    for segment_label, positions in table.groupby("segment", sort=False).groups.items():
        positions_array = np.asarray(list(positions), dtype=int)
        subset = table.loc[positions_array]
        rate_bad = robust_outliers(subset["event_rate"].to_numpy(dtype=float), z_limit)
        signal_bad = np.column_stack(
            [robust_outliers(subset[signal].to_numpy(dtype=float), z_limit) for signal in QC_SIGNALS]
        )
        table.loc[positions_array, "bad"] = rate_bad | (signal_bad.sum(axis=1) >= 2)

    good = np.ones(len(events), dtype=bool)
    for is_bad, indices in zip(table["bad"].to_numpy(dtype=bool), event_bins):
        if is_bad:
            good[indices] = False
    return TimeQCResult(good_events=good, bins=table)
