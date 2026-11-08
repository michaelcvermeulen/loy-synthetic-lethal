from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import pandas as pd


@dataclass(frozen=True)
class FCSData:
    """Events and metadata loaded without modifying the source FCS file."""

    path: Path
    events: pd.DataFrame
    metadata: Mapping[str, str]
    matrix_name: str | None

