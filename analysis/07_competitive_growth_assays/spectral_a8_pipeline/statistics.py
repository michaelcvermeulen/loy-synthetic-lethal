from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class CompetitionResult:
    loy: int
    ctrl: int
    assignable: int
    gfp_only: int
    mcherry_only: int
    double_positive: int
    double_negative: int
    loy_fraction: float
    ratio: float
    log2_ratio: float
    wilson_low: float
    wilson_high: float


def wilson_interval(successes: int, total: int, z: float = 1.96) -> tuple[float, float]:
    if total < 0 or successes < 0 or successes > total:
        raise ValueError("Wilson counts must satisfy 0 <= successes <= total")
    if total == 0:
        return math.nan, math.nan
    proportion = successes / total
    denominator = 1.0 + z**2 / total
    center = (proportion + z**2 / (2.0 * total)) / denominator
    half = z * math.sqrt(proportion * (1.0 - proportion) / total + z**2 / (4.0 * total**2)) / denominator
    return center - half, center + half


def competition_metrics(classes: np.ndarray, loy_reporter: str) -> CompetitionResult:
    labels = np.asarray(classes, dtype=object)
    gfp = int(np.sum(labels == "GFP-only"))
    mcherry = int(np.sum(labels == "mCherry-only"))
    double_positive = int(np.sum(labels == "Double-positive"))
    double_negative = int(np.sum(labels == "Double-negative"))
    if loy_reporter == "GFP":
        loy, ctrl = gfp, mcherry
    elif loy_reporter == "mCherry":
        loy, ctrl = mcherry, gfp
    else:
        raise ValueError("loy_reporter must be 'GFP' or 'mCherry'")
    assignable = loy + ctrl
    fraction = loy / assignable if assignable else math.nan
    ratio = loy / ctrl if ctrl else (math.inf if loy else math.nan)
    log2_ratio = math.log2(ratio) if ratio > 0 else (-math.inf if ratio == 0 else math.nan)
    low, high = wilson_interval(loy, assignable)
    return CompetitionResult(
        loy=loy,
        ctrl=ctrl,
        assignable=assignable,
        gfp_only=gfp,
        mcherry_only=mcherry,
        double_positive=double_positive,
        double_negative=double_negative,
        loy_fraction=fraction,
        ratio=ratio,
        log2_ratio=log2_ratio,
        wilson_low=low,
        wilson_high=high,
    )
