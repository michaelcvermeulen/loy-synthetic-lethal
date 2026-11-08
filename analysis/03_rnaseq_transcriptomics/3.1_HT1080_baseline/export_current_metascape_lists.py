"""Export current HT1080 differential-expression lists for Metascape.

The foreground lists use the thesis primary gene-level rule:
overall BH FDR < 0.05, absolute log2 fold change >= 1, and matching
effect directions in August 2025 and May 2026. The fold-change filter is
post hoc; the FDR pertains to the non-zero-effect test, not a TREAT test.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

import pandas as pd


ROOT = Path(__file__).resolve().parent
TABLES = ROOT / "outputs" / "programme_analysis_v1" / "tables"
OUT = ROOT / "outputs" / "metascape_current_de"
SOURCE_SELECTED = TABLES / "directionally_consistent_genes.csv"
SOURCE_ALL = TABLES / "DE_overall.csv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_symbols(path: Path, values: pd.Series) -> None:
    path.write_text("\n".join(values.astype(str)) + "\n", encoding="utf-8")


def main() -> None:
    selected = pd.read_csv(SOURCE_SELECTED)
    tested = pd.read_csv(SOURCE_ALL)

    required = {
        "gene_id", "gene_name", "logFC", "FDR", "AUG_logFC", "MAY_logFC",
        "passes_main_cutoffs", "directionally_consistent", "direction",
    }
    missing = required.difference(selected.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    selected = selected.loc[
        selected["passes_main_cutoffs"].eq(True)
        & selected["directionally_consistent"].eq(True)
    ].copy()
    selected["metascape_direction"] = selected["logFC"].map(
        lambda value: "UP" if value > 0 else "DOWN"
    )

    assert selected["gene_name"].notna().all()
    assert selected["gene_name"].is_unique
    assert ((selected["logFC"] > 0) == selected["metascape_direction"].eq("UP")).all()
    assert ((selected["AUG_logFC"] > 0) == (selected["MAY_logFC"] > 0)).all()
    assert (selected["FDR"] < 0.05).all()
    assert (selected["logFC"].abs() >= 1).all()

    up = selected.loc[selected["metascape_direction"].eq("UP")].sort_values(
        ["FDR", "gene_name"], kind="stable"
    )
    down = selected.loc[selected["metascape_direction"].eq("DOWN")].sort_values(
        ["FDR", "gene_name"], kind="stable"
    )
    background = tested.loc[tested["FDR"].notna(), "gene_name"].dropna().drop_duplicates()
    background = background.sort_values(kind="stable")

    assert len(up) == 229
    assert len(down) == 433
    assert len(selected) == 662
    assert len(background) == 25461
    assert set(selected["gene_name"]).issubset(set(background))

    OUT.mkdir(parents=True, exist_ok=True)
    files = {
        "up": OUT / "HT1080_current_UP_in_Y_depleted.txt",
        "down": OUT / "HT1080_current_DOWN_in_Y_depleted.txt",
        "background": OUT / "HT1080_current_TESTED_BACKGROUND.txt",
        "audit": OUT / "HT1080_current_selected_DE_genes_with_statistics.csv",
        "two_column": OUT / "HT1080_current_UP_and_DOWN_gene_lists.csv",
    }
    write_symbols(files["up"], up["gene_name"])
    write_symbols(files["down"], down["gene_name"])
    write_symbols(files["background"], background)

    audit_columns = [
        "gene_id", "gene_name", "gene_type", "chrom", "start", "end",
        "mean_CPM", "logFC", "CI.L", "CI.R", "PValue", "FDR",
        "AUG_logFC", "AUG_FDR", "MAY_logFC", "MAY_FDR",
        "interaction_logFC", "interaction_FDR", "both_periods_FDR05",
        "metascape_direction", "statistical_scope",
    ]
    selected.sort_values(["metascape_direction", "FDR", "gene_name"])[audit_columns].to_csv(
        files["audit"], index=False
    )
    width = max(len(up), len(down))
    pd.DataFrame(
        {
            "UP_in_Y_depleted": up["gene_name"].tolist() + [""] * (width - len(up)),
            "DOWN_in_Y_depleted": down["gene_name"].tolist() + [""] * (width - len(down)),
        }
    ).to_csv(files["two_column"], index=False)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "contrast": "Y-depleted minus sgCTRL",
        "model": "limma-voom batch + condition; current programme-analysis workflow",
        "selection": {
            "overall_BH_FDR": "< 0.05",
            "absolute_log2FC": ">= 1",
            "period_direction": "matching signs in August 2025 and May 2026",
            "magnitude_test_warning": (
                "The fold-change threshold is a post-hoc filter. BH FDR controls the "
                "non-zero-effect test, not the claim that absolute log2FC exceeds 1."
            ),
        },
        "scope": "autosomal and non-PAR X genes; Y-linked and PAR genes excluded upstream",
        "counts": {
            "up": len(up),
            "down": len(down),
            "selected_total": len(selected),
            "tested_rows_with_finite_FDR": int(tested["FDR"].notna().sum()),
            "tested_unique_symbols": len(background),
            "duplicate_symbol_rows_collapsed_in_background": int(
                tested.loc[tested["FDR"].notna(), "gene_name"].duplicated().sum()
            ),
        },
        "sources": {
            str(SOURCE_SELECTED.relative_to(ROOT)): sha256(SOURCE_SELECTED),
            str(SOURCE_ALL.relative_to(ROOT)): sha256(SOURCE_ALL),
        },
        "outputs": {},
    }
    for label, path in files.items():
        manifest["outputs"][label] = {
            "file": path.name,
            "sha256": sha256(path),
        }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    readme = f"""# HT1080 current differential-expression lists for Metascape

These lists come from the current HT1080 limma-voom batch + condition model.
The contrast is Y-depleted minus sgCTRL. Genes were selected at overall BH FDR
< 0.05 and absolute log2 fold change >= 1, with matching effect directions in
August 2025 and May 2026.

Counts: {len(up)} UP, {len(down)} DOWN, {len(selected)} total. The custom
background contains {len(background)} unique tested gene symbols. Upload the UP
and DOWN text files as separate lists and use the tested-background text file as
the custom background. Each text file contains one HGNC-style symbol per line
with no header.

The fold-change threshold is a post-hoc filter. The reported FDR applies to the
non-zero-effect test, not a formal test that absolute log2FC exceeds 1. The
analysis is restricted to autosomal and non-PAR X genes; Y-linked and PAR genes
were excluded upstream. These samples are serial passages of one established
population per condition, so the list does not demonstrate reproducibility
across independently generated Y-depleted populations.

`HT1080_current_selected_DE_genes_with_statistics.csv` preserves the model
statistics and period effects. `manifest.json` records input/output hashes,
counts and the exact selection rule.

Reproduce from this analysis directory:

```powershell
python -X utf8 export_current_metascape_lists.py
```
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")

    print(json.dumps(manifest["counts"], indent=2))


if __name__ == "__main__":
    main()
