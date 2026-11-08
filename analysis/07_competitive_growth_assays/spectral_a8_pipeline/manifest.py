from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Iterable


SOURCE_SUFFIXES = {".fcs", ".cvw", ".tiff", ".flowjo"}
OUTPUT_DIR_NAMES = {"programmatic_analysis", ".programmatic_analysis_staging"}


def hash_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class SourceRecord:
    path: Path
    size: int
    mtime_ns: int
    sha256: str

    @classmethod
    def from_path(cls, path: Path) -> "SourceRecord":
        resolved = path.resolve()
        stat = resolved.stat()
        return cls(resolved, stat.st_size, stat.st_mtime_ns, hash_file(resolved))


def build_manifest(root: Path) -> list[SourceRecord]:
    resolved = root.resolve()
    paths = (
        path
        for path in resolved.rglob("*")
        if path.is_file()
        and path.suffix.lower() in SOURCE_SUFFIXES
        and not any(part in OUTPUT_DIR_NAMES for part in path.relative_to(resolved).parts)
    )
    return [SourceRecord.from_path(path) for path in sorted(paths)]


def verify_manifest(records: Iterable[SourceRecord]) -> None:
    for record in records:
        if not record.path.exists():
            raise RuntimeError(f"source changed: missing {record.path}")
        stat = record.path.stat()
        if stat.st_size != record.size or stat.st_mtime_ns != record.mtime_ns:
            raise RuntimeError(f"source changed: metadata differs for {record.path}")
        if hash_file(record.path) != record.sha256:
            raise RuntimeError(f"source changed: content differs for {record.path}")

