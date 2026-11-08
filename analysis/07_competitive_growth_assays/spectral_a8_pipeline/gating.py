from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import theilslopes


@dataclass(frozen=True)
class LineFit:
    intercept: float
    slope: float
    residual_center: float
    residual_scale: float
    valid_count: int


def points_in_polygon(points: np.ndarray, polygon: np.ndarray) -> np.ndarray:
    xy = np.asarray(points, dtype=float)
    vertices = np.asarray(polygon, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError("points must have shape (n, 2)")
    if vertices.ndim != 2 or vertices.shape[1] != 2 or len(vertices) < 3:
        raise ValueError("polygon must contain at least three 2D vertices")
    x, y = xy[:, 0], xy[:, 1]
    inside = np.zeros(len(xy), dtype=bool)
    xj, yj = vertices[-1]
    for xi, yi in vertices:
        crosses = ((yi > y) != (yj > y)) & (
            x < (xj - xi) * (y - yi) / ((yj - yi) + np.finfo(float).eps) + xi
        )
        inside ^= crosses
        xj, yj = xi, yi
    return inside


def scale_polygon(polygon: np.ndarray, factor: float) -> np.ndarray:
    vertices = np.asarray(polygon, dtype=float)
    if factor <= 0:
        raise ValueError("polygon scale factor must be positive")
    centroid = vertices.mean(axis=0)
    return centroid + factor * (vertices - centroid)


def _scaled_mad(values: np.ndarray) -> float:
    center = float(np.median(values))
    return float(1.4826 * np.median(np.abs(values - center)))


def pulse_relationship_mask(
    area: np.ndarray | pd.Series,
    height: np.ndarray | pd.Series,
    mad_limit: float = 3.5,
) -> tuple[np.ndarray, LineFit]:
    area_values = np.asarray(area, dtype=float)
    height_values = np.asarray(height, dtype=float)
    if area_values.shape != height_values.shape:
        raise ValueError("pulse area and height must have identical shapes")
    valid = np.isfinite(area_values) & np.isfinite(height_values) & (area_values > 0) & (height_values > 0)
    if valid.sum() < 3:
        raise ValueError("at least three positive finite pulse pairs are required")
    x = np.log10(area_values[valid])
    y = np.log10(height_values[valid])
    order = np.argsort(x)
    sample_positions = np.linspace(0, len(order) - 1, min(len(order), 2000), dtype=int)
    sample = order[sample_positions]
    robust = theilslopes(y[sample], x[sample])
    intercept, slope = float(robust.intercept), float(robust.slope)
    residual = y - (intercept + slope * x)
    center = float(np.median(residual))
    scale = _scaled_mad(residual)
    tolerance = mad_limit * scale if scale > 0 else 1e-8
    valid_keep = np.abs(residual - center) <= tolerance
    keep = np.zeros(area_values.shape, dtype=bool)
    keep[np.flatnonzero(valid)] = valid_keep
    return keep, LineFit(intercept, slope, center, scale, int(valid.sum()))


def pulse_singlet_mask(
    events: pd.DataFrame,
    mad_limit: float = 3.5,
) -> tuple[np.ndarray, dict[str, LineFit]]:
    required = {"FSC-A", "FSC-H", "SSC (Violet)-A", "SSC (Violet)-H"}
    missing = sorted(required - set(events.columns))
    if missing:
        raise KeyError(f"missing pulse parameters: {missing}")
    fsc, fsc_fit = pulse_relationship_mask(events["FSC-A"], events["FSC-H"], mad_limit)
    ssc, ssc_fit = pulse_relationship_mask(events["SSC (Violet)-A"], events["SSC (Violet)-H"], mad_limit)
    return fsc & ssc, {"fsc": fsc_fit, "ssc": ssc_fit}
