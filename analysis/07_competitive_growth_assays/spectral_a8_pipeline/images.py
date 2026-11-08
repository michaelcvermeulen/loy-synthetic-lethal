"""Map extracted A8 TIFFs to FCS events and render deterministic review montages.

Two properties of this export differ from the original plan assumptions and are
handled explicitly here:

1. The eight-digit filename suffix is a ZERO-based event index, not a one-based
   event number. G11 contains exactly 20,000 files numbered 00000000-00019999
   for 20,000 events, so a one-based reading is impossible. parse_event_number
   returns the literal suffix; build_image_manifest converts it to the one-based
   event number used everywhere else in the pipeline.

2. PIL.ImageSequence.Iterator yields stale pixel data for these multi-frame
   float TIFFs. Frames must be read with an explicit seek followed by load,
   which read_tiff_channels does.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

EVENT_SUFFIX = re.compile(r"_(\d{8})\.tiff$", re.IGNORECASE)

# Frame index of each imaging channel inside the extracted TIFF stacks. The TIFFs
# carry no channel metadata, so this was established empirically by correlating
# each frame's peak against the fully populated per-channel "Max Intensity"
# parameters in the FCS; see the Task 6 entry in implementation_log.md.
#
#   frame 1 -> FSC                       Spearman 0.939
#   frame 2 -> SSC (Imaging)             Spearman 0.967
#   frame 3 -> GFP_reporter eGFP         Spearman 0.944
#   frame 5 -> mCherry_reporter mCherry  Spearman 0.950
#
# Frame 5 is independent of GFP (rho -0.058), so it is a usable second opinion on
# the spectral mCherry axis. Frames 0 and 4 are not confidently assigned.
DEFAULT_CHANNEL_MAP: Mapping[str, int] = {"gray": 1, "green": 3, "red": 5}


def parse_event_number(path: Path) -> int:
    """Return the literal zero-based event index encoded in a TIFF filename."""
    match = EVENT_SUFFIX.search(Path(path).name)
    if match is None:
        raise ValueError(f"filename does not encode an event index: {path}")
    return int(match.group(1))


def build_image_manifest(root: Path, event_count: int) -> pd.DataFrame:
    """Map every FCS event to its TIFF, leaving missing images explicit.

    event_number is one-based. image_index is the zero-based filename suffix.
    Events without an extracted image get a null path and remain valid rows so
    that image coverage is reported, never used as an exclusion.
    """
    resolved = Path(root)
    if event_count < 1:
        raise ValueError("event_count must be positive")

    found: dict[int, Path] = {}
    for path in sorted(resolved.rglob("*.tiff")):
        index = parse_event_number(path)
        if index in found:
            raise ValueError(
                f"duplicate image event {index}: {found[index].name} and {path.name}"
            )
        if not 0 <= index < event_count:
            raise ValueError(
                f"image event {index} is out of range for {event_count} events: {path.name}"
            )
        found[index] = path

    indices = np.arange(event_count, dtype=int)
    paths = [found.get(int(index)) for index in indices]
    return pd.DataFrame(
        {
            "event_number": indices + 1,
            "image_index": indices,
            "path": paths,
            "available": [path is not None for path in paths],
        }
    )


def read_tiff_channels(path: Path) -> np.ndarray:
    """Read every frame of a multi-frame TIFF as (frames, height, width)."""
    with Image.open(Path(path)) as image:
        frames = []
        for index in range(getattr(image, "n_frames", 1)):
            image.seek(index)
            image.load()
            frames.append(np.array(image, dtype=np.float32))
    shapes = {frame.shape for frame in frames}
    if len(shapes) != 1:
        raise ValueError(f"frames disagree in shape for {path}: {sorted(shapes)}")
    return np.stack(frames)


def _normalize(values: np.ndarray, low: float, high: float) -> np.ndarray:
    if high <= low:
        return np.zeros_like(values, dtype=np.float32)
    return np.clip((values - low) / (high - low), 0.0, 1.0).astype(np.float32)


def _pad(tile: np.ndarray, height: int, width: int) -> np.ndarray:
    padded = np.zeros((height, width, tile.shape[2]), dtype=tile.dtype)
    top = (height - tile.shape[0]) // 2
    left = (width - tile.shape[1]) // 2
    padded[top : top + tile.shape[0], left : left + tile.shape[1]] = tile
    return padded


def render_montage(
    groups: Mapping[str, Sequence[tuple[int, Path]]],
    output: Path,
    channel_map: Mapping[str, int] = DEFAULT_CHANNEL_MAP,
    percentiles: tuple[float, float] = (50.0, 99.5),  # median black point: >50% of each frame is background
    label_height: int = 14,
    margin: int = 110,
) -> Path:
    """Render a labelled montage with one row per population.

    Display normalization uses fixed percentiles computed once across every tile
    in the montage, so tiles are comparable within a montage. Normalization is a
    display choice only and changes no quantitative FCS measurement.
    """
    if not groups:
        raise ValueError("at least one population is required")

    stacks: dict[str, list[tuple[int, np.ndarray]]] = {}
    for label, entries in groups.items():
        stacks[label] = [(event, read_tiff_channels(path)) for event, path in entries]

    populated = [row for row in stacks.values() if row]
    if not populated:
        raise ValueError("no images were supplied for any population")

    used = sorted(set(channel_map.values()))
    available = min(stack.shape[0] for row in populated for _, stack in row)
    if used and max(used) >= available:
        raise ValueError(
            f"channel_map references frame {max(used)} but stacks have {available} frames"
        )

    limits: dict[int, tuple[float, float]] = {}
    for frame_index in used:
        pooled = np.concatenate(
            [stack[frame_index].ravel() for row in populated for _, stack in row]
        )
        limits[frame_index] = (
            float(np.percentile(pooled, percentiles[0])),
            float(np.percentile(pooled, percentiles[1])),
        )

    tiles: dict[str, list[tuple[int, np.ndarray]]] = {}
    for label, row in stacks.items():
        rendered = []
        for event, stack in row:
            height, width = stack.shape[1], stack.shape[2]
            rgb = np.zeros((height, width, 3), dtype=np.float32)
            for name, frame_index in channel_map.items():
                low, high = limits[frame_index]
                layer = _normalize(stack[frame_index], low, high)
                if name == "gray":
                    rgb += layer[..., None] * np.array([0.55, 0.55, 0.55], dtype=np.float32)
                elif name == "green":
                    rgb[..., 1] += layer
                elif name == "red":
                    rgb[..., 0] += layer
                elif name == "blue":
                    rgb[..., 2] += layer
                else:
                    raise ValueError(f"unsupported montage channel: {name}")
            rendered.append((event, np.clip(rgb, 0.0, 1.0)))
        tiles[label] = rendered

    tile_height = max(tile.shape[0] for row in tiles.values() for _, tile in row)
    tile_width = max(tile.shape[1] for row in tiles.values() for _, tile in row)
    columns = max(len(row) for row in tiles.values())
    canvas_width = margin + columns * (tile_width + 4)
    canvas_height = 18 + len(tiles) * (tile_height + label_height + 8)

    canvas = Image.new("RGB", (canvas_width, canvas_height), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    draw.text((6, 4), f"channels {dict(channel_map)}", fill=(0, 0, 0))

    y = 18
    for label, row in tiles.items():
        draw.text((6, y + label_height + tile_height // 2), label, fill=(0, 0, 0))
        x = margin
        for event, tile in row:
            padded = _pad(tile, tile_height, tile_width)
            canvas.paste(
                Image.fromarray((padded * 255).astype(np.uint8), mode="RGB"),
                (x, y + label_height),
            )
            draw.text((x, y), f"event {event}", fill=(0, 0, 0))
            x += tile_width + 4
        y += tile_height + label_height + 8

    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(destination)
    return destination


def select_population_images(
    manifest: pd.DataFrame,
    event_numbers: Iterable[int],
    maximum: int = 50,
    seed: int = 20260828,
) -> list[tuple[int, Path]]:
    """Deterministically sample available images for one population."""
    wanted = np.asarray(list(event_numbers), dtype=int)
    lookup = manifest.set_index("event_number")
    present = [
        int(event)
        for event in wanted
        if event in lookup.index and bool(lookup.loc[event, "available"])
    ]
    rng = np.random.default_rng(seed)
    if len(present) > maximum:
        chosen = rng.choice(np.asarray(present, dtype=int), size=maximum, replace=False)
        present = sorted(int(event) for event in chosen)
    return [(event, Path(lookup.loc[event, "path"])) for event in present]
