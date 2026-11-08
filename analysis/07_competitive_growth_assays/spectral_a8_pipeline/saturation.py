from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


MASK_COLUMNS = ("SaturatedChannels", "SaturatedChannels1", "SaturatedChannels2")


def _spectral_names(prefix: str, wavelengths: Sequence[int]) -> list[str]:
    return [f"{prefix}{index} ({wavelength})" for index, wavelength in enumerate(wavelengths, start=1)]


DETECTOR_NAMES = tuple(
    ["B-SSC", "ImgB1 (535)", "ImgB2 (600)", "ImgB3 (790)", "B-ALL", "B-FSC"]
    + _spectral_names("B", [500, 515, 530, 545, 575, 590, 605, 625, 655, 675, 700, 725, 750, 780, 810, 845])
    + _spectral_names("V", [420, 440, 460, 475, 500, 515, 530, 545, 575, 590, 605, 625, 655, 675, 700, 725, 750, 780, 810, 845])
    + ["V-SSC", "V-ALL"]
    + _spectral_names("YG", [575, 590, 605, 625, 655, 675, 700, 725, 750, 780, 810, 845])
    + _spectral_names("R", [655, 675, 700, 725, 750, 780, 810, 845])
    + _spectral_names("UV", [375, 390, 420, 440, 460, 475, 500, 515, 530, 545, 575, 590, 605, 625, 655, 675, 700, 725, 750, 780, 810, 845])
)

IMAGING_SCATTER_INDICES = np.array([0, 1, 2, 3, 4, 5, 42, 43], dtype=int)
SCATTER_LIGHTLOSS_INDICES = np.array([0, 4, 5, 42, 43], dtype=int)
SPECTRAL_INDICES = np.array([i for i in range(len(DETECTOR_NAMES)) if i not in set(IMAGING_SCATTER_INDICES)], dtype=int)


def reinterpret_u32(values: np.ndarray) -> np.ndarray:
    array = np.ascontiguousarray(values, dtype=np.float32)
    return array.view(np.uint32)


def decode_mask_matrix(stored: np.ndarray) -> np.ndarray:
    array = np.asarray(stored)
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError("saturation storage must have shape (events, 3)")
    words = reinterpret_u32(array).reshape(-1, 3)
    bit = np.arange(32, dtype=np.uint32)
    return np.concatenate(
        [((words[:, word, None] >> bit) & np.uint32(1)).astype(bool) for word in range(3)],
        axis=1,
    )


def decode_saturation(
    events: pd.DataFrame,
    detector_names: Sequence[str] = DETECTOR_NAMES,
) -> pd.DataFrame:
    missing = sorted(set(MASK_COLUMNS) - set(events.columns))
    if missing:
        raise KeyError(f"missing saturation parameters: {missing}")
    if len(detector_names) > 96:
        raise ValueError("detector mapping exceeds 96 saturation bits")

    decoded = decode_mask_matrix(events.loc[:, MASK_COLUMNS].to_numpy())
    active = decoded[:, : len(detector_names)]
    any_saturated = decoded.any(axis=1)
    detector_count = active.sum(axis=1)
    return pd.DataFrame(
        {
            "any_saturated": any_saturated,
            "spectral_saturated": active[:, SPECTRAL_INDICES].any(axis=1),
            "imgb1_saturated": active[:, 1],
            "imgb2_saturated": active[:, 2],
            "imgb3_saturated": active[:, 3],
            "imgb3_only_saturated": active[:, 3] & (detector_count == 1),
            "scatter_lightloss_saturated": active[:, SCATTER_LIGHTLOSS_INDICES].any(axis=1),
            "imaging_scatter_saturated": active[:, IMAGING_SCATTER_INDICES].any(axis=1),
        },
        index=events.index,
    )
