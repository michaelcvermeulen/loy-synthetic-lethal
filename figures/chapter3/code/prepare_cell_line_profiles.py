"""Extract release-specific DepMap figure data with explicit, externally supplied Y calls.

This script does not infer LOY calls or borrow statistics from another release.
The output directory should be empty (or a new staging directory).
"""
from pathlib import Path
import argparse
import hashlib
import json
import re

import numpy as np
import pandas as pd

Y_GENES = ["DDX3Y", "EIF1AY", "RPS4Y1", "UTY", "KDM5D", "USP9Y", "ZFY", "PRKY", "TMSB4Y", "NLGN4Y"]
X_GENES = ["DDX3X", "EIF1AX", "RPS4X", "ZFX"]
MODELS = {"HT1080": "HT1080", "A549": "A549", "DU145": "DU145", "NCIH2452": "H2452"}


def symbol(column):
    return re.sub(r" \([^)]*\)$", "", str(column))


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def default_rows(frame, label):
    if "IsDefaultEntryForModel" in frame:
        keep = frame.IsDefaultEntryForModel.astype(str).str.lower().isin(["yes", "true", "1"])
        frame = frame.loc[keep].copy()
    if frame.ModelID.duplicated().any():
        raise ValueError(f"{label}: duplicate models after selecting default entries")
    return frame


def locate(root, name, section, required=True):
    candidates = [root / name, root / 'raw/depmap_public_26q1' / section / name]
    found = [path for path in candidates if path.is_file()]
    if len(found) > 1:
        raise ValueError(f'Ambiguous input {name}: {found}')
    if not found and required:
        raise FileNotFoundError(f'{name} has not arrived in {root}')
    return found[0] if found else None


def prepare(args):
    paths = {
        "metadata": locate(args.data_dir, "Model.csv", "metadata"),
        "gene_effect": locate(args.data_dir, "CRISPRGeneEffect.csv", "crispr"),
        "expression": locate(args.data_dir, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", "expression"),
    }
    readme_path = locate(args.data_dir, "README.txt", "documentation", required=False)
    source_manifest_path = args.data_dir / 'docs/data_file_manifest.csv'
    if readme_path:
        paths['readme'] = readme_path
    elif source_manifest_path.is_file():
        paths['source_manifest'] = source_manifest_path
    else:
        raise FileNotFoundError('Neither the release README nor the HPC file manifest has arrived')
    if args.calls:
        paths['calls'] = args.calls
    controls_path = locate(args.data_dir, "AchillesCommonEssentialControls.csv", "crispr", required=False)
    if controls_path:
        paths['common_essential'] = controls_path
    for label, path in paths.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"{label} is absent or empty: {path}")
    before = {label: (path.stat().st_size, path.stat().st_mtime_ns) for label, path in paths.items()}
    if readme_path:
        readme = readme_path.read_text(encoding="utf-8-sig")
        if not re.search(re.escape(args.release), readme[:2000], re.I):
            raise ValueError(f"README does not identify the requested {args.release} release")
    source_manifest = pd.read_csv(source_manifest_path) if source_manifest_path.is_file() else None

    model = pd.read_csv(paths["metadata"])
    if model.ModelID.duplicated().any():
        raise ValueError("Model.csv has duplicate ModelIDs")
    stripped = model.CellLineName.str.replace(r"[^A-Za-z0-9]", "", regex=True).str.upper()
    model["cell_line"] = stripped.map(MODELS).fillna(model.CellLineName)
    selected = model.loc[stripped.isin(MODELS)].copy()
    if len(selected) != 4 or set(selected.cell_line) != set(MODELS.values()):
        raise ValueError("Cannot uniquely resolve HT1080, A549, DU145 and NCI-H2452")
    ht_id = selected.loc[selected.cell_line.eq("HT1080"), "ModelID"].item()

    expr_header = pd.read_csv(paths["expression"], nrows=0).columns
    expr_columns = [c for c in expr_header if symbol(c) in Y_GENES]
    absent_y = sorted(set(Y_GENES) - {symbol(c) for c in expr_columns})
    id_columns = [c for c in ["ModelID", "IsDefaultEntryForModel", "SequencingID"] if c in expr_header]
    expr = default_rows(pd.read_csv(paths["expression"], usecols=id_columns + expr_columns), "Expression")
    expr = expr.rename(columns={c: symbol(c) for c in expr_columns})
    if expr.columns.duplicated().any():
        raise ValueError("Ambiguous expression gene symbols")
    ht_expr = expr.loc[expr.ModelID.eq(ht_id)]
    if len(ht_expr) != 1:
        raise ValueError("No unique default HT1080 expression profile")
    ht_expr = ht_expr.melt(id_vars=id_columns, value_vars=[g for g in Y_GENES if g in expr],
                           var_name="gene", value_name="expression")
    ht_expr["units"] = "log2(TPM + 1)"
    selected_expr = expr.loc[expr.ModelID.isin(selected.ModelID)].merge(
        selected[['ModelID', 'cell_line']], on='ModelID', validate='one_to_one')
    selected_expr = selected_expr.melt(id_vars=id_columns+['cell_line'],
                                      value_vars=[g for g in Y_GENES if g in expr],
                                      var_name='gene', value_name='expression')
    selected_expr['units'] = 'log2(TPM + 1)'

    effects = pd.read_csv(paths["gene_effect"], index_col=0)
    effects.index.name = "ModelID"
    if effects.index.duplicated().any():
        raise ValueError("Duplicate models in gene-effect matrix")
    effects.columns = [symbol(c) for c in effects.columns]
    if effects.columns.duplicated().any():
        raise ValueError("Ambiguous gene-effect gene symbols")
    if set(X_GENES) - set(effects.columns):
        raise ValueError("Requested X-linked genes missing from gene-effect matrix")
    if set(selected.ModelID) - set(effects.index):
        raise ValueError("At least one selected model lacks a gene-effect profile")
    if not np.isfinite(effects.loc[selected.ModelID, X_GENES].to_numpy()).all():
        raise ValueError("Selected-model gene effects contain missing/nonfinite values")

    # Descriptive selection of strong, relatively selective baseline dependencies.
    # The full ranking is retained so this working display choice can be reviewed.
    controls = pd.read_csv(controls_path, sep="\t").iloc[:, 0].map(symbol) if controls_path else []
    reference = effects.drop(index=ht_id)
    deps = pd.DataFrame({
        "gene_effect": effects.loc[ht_id],
        "reference_median": reference.median(),
        "reference_q25": reference.quantile(.25),
        "reference_q75": reference.quantile(.75),
        "n_reference": reference.count(),
        "reference_dependency_fraction": reference.lt(-.5).sum() / reference.count(),
    })
    deps.index.name = "gene"
    deps["difference_from_reference_median"] = deps.gene_effect - deps.reference_median
    deps["common_essential"] = deps.index.isin(controls)
    deps = deps.sort_values("difference_from_reference_median")
    shortlist = deps.loc[(~deps.common_essential) & deps.gene_effect.lt(-.5)
                         & deps.reference_dependency_fraction.le(.20)].head(16)
    if len(shortlist) == 0:
        raise ValueError("No HT1080 genes meet the stated display criteria")
    if not controls_path:
        deps['common_essential'] = pd.NA
        shortlist = shortlist.copy()
        shortlist['common_essential'] = pd.NA

    cohort, comparison, selected_calls = None, None, None
    counts = {}
    if args.calls:
        calls = pd.read_csv(paths["calls"])
        if args.call_column not in calls or "ModelID" not in calls:
            raise ValueError("Calls require ModelID and the requested call column")
        if calls.ModelID.duplicated().any():
            raise ValueError("Duplicate ModelIDs in Y calls")
        calls = calls[["ModelID", args.call_column]].rename(columns={args.call_column: "LOY_status"})
        calls["LOY_status"] = calls.LOY_status.astype("string").str.strip()
        cohort = model.merge(calls, on="ModelID", how="left", validate="one_to_one")
        cohort["has_gene_effect"] = cohort.ModelID.isin(effects.index)
        cohort["included"] = cohort.Sex.eq("Male") & cohort.LOY_status.isin(["XY", "LOY"]) & cohort.has_gene_effect
        retained = cohort.loc[cohort.included, ["ModelID", "cell_line", "LOY_status", "OncotreeLineage"]]
        selected_calls = cohort.loc[cohort.ModelID.isin(selected.ModelID)]
        if not selected_calls.included.all():
            raise ValueError("A selected model lacks a male XY/LOY call; inspect it before plotting")
        values = effects[X_GENES].reset_index().melt(id_vars="ModelID", var_name="gene", value_name="gene_effect")
        comparison = retained.merge(values, on="ModelID", validate="one_to_many")
        comparison["highlight"] = comparison.cell_line.isin(MODELS.values())
        counts = retained.LOY_status.value_counts().to_dict()

    manifest = {
        "release_label": f"DepMap Public {args.release}",
        "comparison_label": args.calls_description,
        "calls_description": args.calls_description,
        "call_column": args.call_column,
        "population": "Male-derived models with supplied XY/LOY calls and CRISPR gene effects",
        "reference_population_for_ht1080_dependencies": "All other screened models, both sexes and all lineages",
        "dependency_display_selection": "Exclude supplied common-essential controls, if available; HT1080 gene effect < -0.5; effect < -0.5 in <=20% of other models; show up to 16 most negative differences from the other-model median. Descriptive ranking, not an inferential discovery test.",
        "common_essential_controls_applied": bool(controls_path),
        "absent_y_expression_genes": absent_y,
        "statistics": "Not calculated by this extraction or plotting script",
        "cohort_counts": counts,
        "sources": {},
    }
    if not controls_path:
        manifest['dependency_display_selection'] = manifest['dependency_display_selection'].replace('Exclude common-essential controls; ', '')
    for label, path in paths.items():
        sha = digest(path)
        after = (path.stat().st_size, path.stat().st_mtime_ns)
        if after != before[label]:
            raise RuntimeError(f"Input changed during reading: {path}. Wait for transfer completion.")
        manifest["sources"][label] = {"path": str(path.resolve()), "bytes": after[0], "sha256": sha}
        if source_manifest is not None and label in ['metadata', 'gene_effect', 'expression', 'readme']:
            expected = source_manifest.loc[source_manifest.new_relative_path.eq(path.relative_to(args.data_dir).as_posix())]
            if len(expected) != 1 or sha != expected.iloc[0].sha256:
                raise RuntimeError(f'Input does not match the HPC manifest: {path}')
            manifest['sources'][label]['matches_hpc_manifest'] = True

    args.output_dir.mkdir(parents=True, exist_ok=True)
    targets = ["ht1080_y_expression.csv", "ht1080_dependency_ranking.csv", "ht1080_dependency_context.csv",
               "dependency_comparison.csv", "cohort_inclusion.csv", "selected_model_metadata.csv", "figure_config.json"]
    if any((args.output_dir / filename).exists() for filename in targets):
        raise FileExistsError("Use a new staging directory to preserve existing figure source data")
    ht_expr.to_csv(args.output_dir / targets[0], index=False)
    deps.to_csv(args.output_dir / targets[1])
    shortlist.to_csv(args.output_dir / targets[2])
    if comparison is not None:
        comparison.to_csv(args.output_dir / targets[3], index=False)
        cohort.to_csv(args.output_dir / targets[4], index=False)
    selected.to_csv(args.output_dir / targets[5], index=False)
    selected_expr.to_csv(args.output_dir / 'selected_models_y_expression.csv', index=False)
    (args.output_dir / targets[6]).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"release": args.release, "cohort": manifest["cohort_counts"],
                      "selected": selected[['cell_line', 'ModelID']].to_dict('records'),
                      "missing_y_genes": absent_y, "dependencies": shortlist.index.tolist()}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--release", default="26Q1")
    parser.add_argument("--calls", type=Path, help='Omit to prepare only baseline profiles while calls are pending')
    parser.add_argument("--call-column", default="LOY_binary")
    parser.add_argument("--calls-description", default='Y calls pending; baseline profiles only')
    parser.add_argument("--output-dir", type=Path, required=True)
    prepare(parser.parse_args())
