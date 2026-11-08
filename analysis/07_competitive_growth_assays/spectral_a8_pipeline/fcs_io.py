from __future__ import annotations

from pathlib import Path
from typing import Sequence

import flowio
import numpy as np
import pandas as pd

from .models import FCSData


def load_fcs(path: Path, required: Sequence[str] = ()) -> FCSData:
    source = Path(path).resolve()
    fd = flowio.FlowData(str(source))
    metadata = dict(fd.text)
    names = [metadata[f"p{index}n"] for index in range(1, fd.channel_count + 1)]
    missing = sorted(set(required) - set(names))
    if missing:
        raise KeyError(f"missing required parameters: {missing}")
    values = np.asarray(fd.events, dtype=np.float64).reshape(fd.event_count, fd.channel_count)
    events = pd.DataFrame(values, columns=names)
    return FCSData(
        path=source,
        events=events,
        metadata=metadata,
        matrix_name=metadata.get("bdspectral matrix name"),
    )
