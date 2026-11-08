"""Transparent guide filtering from saved tables (Python standard library only).

Every candidate receives a decision. Missing mismatch counts remain unknown.
No genome search, CFD recalculation or historical selection claim is made.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys

VERSION = "1.0.0"
PROJECT = Path(__file__).resolve().parents[1]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrMT"]
COUNT_COLS = CHROMS + ["rev_" + c for c in CHROMS]
MISSING = {"", "NA", "N/A", "NAN", "NULL", "NONE"}
FILTER_KEYS = ("minimum_exact_y_matches", "maximum_exact_non_y_matches",
               "maximum_non_y_1_mismatch_sites", "maximum_non_y_2_mismatch_sites")
RANKING = "exact_y_matches_descending_then_spacer_alphabetical"
OUTPUT_COLUMNS = ["rank", "guide_id", "guide", "region", "exact_y_matches",
                  "exact_y_forward", "exact_y_reverse", "exact_non_y_matches",
                  "non_y_mm1", "non_y_mm2", "non_y_mm3", "non_y_mm4",
                  "legacy_non_y_cfd_mean", "legacy_non_y_cfd_max", "legacy_non_y_cfd_sum",
                  "annotation_row_present", "required_mismatch_counts_available",
                  "decision", "reasons"]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def number(value, field, optional=False, integer=False):
    if value is None or str(value).strip().upper() in MISSING:
        if optional:
            return None
        raise ValueError(f"Missing required value: {field}")
    try:
        value = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Invalid number for {field}") from exc
    if not math.isfinite(value) or value < 0 or (integer and value != int(value)):
        raise ValueError(f"Invalid non-negative {'integer' if integer else 'number'}: {field}")
    return int(value) if integer else value


def read_table(path, required):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: missing columns {sorted(set(required) - set(reader.fieldnames or []))}")
        rows = list(reader)
    seen = set()
    for row in rows:
        guide = row["guide"].strip().upper()
        if len(guide) != 20 or set(guide) - set("ACGT"):
            raise ValueError(f"{path}: invalid 20-nt spacer {guide!r}")
        if guide in seen:
            raise ValueError(f"{path}: duplicate spacer {guide}")
        seen.add(guide)
        row["guide"] = guide
    return rows


def load_candidates(path):
    rows = read_table(path, ["guide", "sumY", "sumOther"] + COUNT_COLS)
    if not rows:
        raise ValueError("Candidate table is empty")
    for row in rows:
        for field in COUNT_COLS + ["sumY", "sumOther"]:
            row[field] = number(row[field], field, integer=True)
        if row["sumY"] != row["chrY"] + row["rev_chrY"]:
            raise ValueError(f"Inconsistent Y count: {row['guide']}")
        non_y = sum(row[c] for c in COUNT_COLS if c not in ("chrY", "rev_chrY"))
        if row["sumOther"] != non_y:
            raise ValueError(f"Inconsistent non-Y count: {row['guide']}")
    return rows


def load_annotations(path):
    rows = read_table(path, ["guide", "mismatch_1", "mismatch_2"])
    out = {}
    for row in rows:
        if "chrY_rows" in row and number(row["chrY_rows"], "chrY_rows", integer=True) != 0:
            raise ValueError("Non-Y annotations contain chrY site rows")
        parsed = {f"non_y_mm{i}": number(row.get(f"mismatch_{i}"), f"mismatch_{i}",
                                          optional=True, integer=True) for i in range(1, 5)}
        for metric in ("mean", "max", "sum"):
            parsed[f"legacy_non_y_cfd_{metric}"] = number(row.get(f"mismatch_{metric}_CFD"),
                                                          f"mismatch_{metric}_CFD", optional=True)
        out[row["guide"]] = parsed
    return out


def evaluate(candidates, annotations, labels, filters):
    decisions = []
    for candidate in candidates:
        guide = candidate["guide"]
        annotation = annotations.get(guide, {})
        label = labels.get(guide, {})
        row = dict.fromkeys(OUTPUT_COLUMNS)
        row.update(guide=guide, guide_id=label.get("guide_id", ""), region=label.get("region", ""),
                   exact_y_matches=candidate["sumY"], exact_y_forward=candidate["chrY"],
                   exact_y_reverse=candidate["rev_chrY"], exact_non_y_matches=candidate["sumOther"],
                   annotation_row_present=guide in annotations)
        row.update(annotation)
        reasons = []
        if row["exact_non_y_matches"] > filters["maximum_exact_non_y_matches"]:
            reasons.append("exact_non_y_matches_above_limit")
        if row["exact_y_matches"] < filters["minimum_exact_y_matches"]:
            reasons.append("exact_y_matches_below_minimum")
        row["required_mismatch_counts_available"] = all(row[f"non_y_mm{i}"] is not None for i in (1, 2))
        if reasons:
            row["decision"] = "rejected_exact_match_filters"
        else:
            for i in (1, 2):
                count = row[f"non_y_mm{i}"]
                if count is not None and count > filters[f"maximum_non_y_{i}_mismatch_sites"]:
                    reasons.append(f"non_y_{i}_mismatch_sites_above_limit")
            if reasons:
                row["decision"] = "rejected_mismatch_filters"
            elif not row["required_mismatch_counts_available"]:
                row["decision"] = "unassessed_missing_mismatch_counts"
                reasons = [f"missing_non_y_mm{i}" for i in (1, 2) if row[f"non_y_mm{i}"] is None]
            else:
                row["decision"] = "retained"
                reasons = ["passes_all_configured_filters"]
        row["reasons"] = ";".join(reasons)
        decisions.append(row)
    # Explicit deterministic tie-break; no score-derived or manually forced ranks.
    decisions.sort(key=lambda r: (-r["exact_y_matches"], r["guide"]))
    retained = [r for r in decisions if r["decision"] == "retained"]
    for rank, row in enumerate(retained, 1):
        row["rank"] = rank
    return decisions


def write_csv(path, rows, columns):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def stage_counts(rows, filters):
    exact = [r for r in rows if r["exact_non_y_matches"] <= filters["maximum_exact_non_y_matches"]]
    eligible = [r for r in exact if r["exact_y_matches"] >= filters["minimum_exact_y_matches"]]
    assessed = [r for r in eligible if r["required_mismatch_counts_available"]]
    mm1 = [r for r in assessed if r["non_y_mm1"] <= filters["maximum_non_y_1_mismatch_sites"]]
    mm2 = [r for r in mm1 if r["non_y_mm2"] <= filters["maximum_non_y_2_mismatch_sites"]]
    stages = [
        ("input", "input", len(rows), "Distinct 20-nt guide sequences"),
        ("exact_non_y_filter", "filter", len(exact), f"Exact non-Y matches <= {filters['maximum_exact_non_y_matches']}"),
        ("exact_y_abundance_filter", "filter", len(eligible), f"Exact Y matches >= {filters['minimum_exact_y_matches']}"),
        ("required_mismatch_counts_available", "data_coverage", len(assessed), "Eligible guides with both one- and two-mismatch counts"),
        ("required_mismatch_counts_missing", "data_coverage", len(eligible)-len(assessed), "Eligible guides lacking either count; missing is not zero"),
        ("one_mismatch_filter", "filter_within_assessed_guides", len(mm1), f"Non-Y one-mismatch sites <= {filters['maximum_non_y_1_mismatch_sites']}"),
        ("two_mismatch_filter", "filter_within_assessed_guides", len(mm2), f"Non-Y two-mismatch sites <= {filters['maximum_non_y_2_mismatch_sites']}"),
    ]
    return [dict(stage=s, kind=k, guides=n, description=d) for s, k, n, d in stages]


def run(config_path, output=None, allow_partial=False):
    config_path = Path(config_path).resolve()
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    if config.get("pipeline_version") != VERSION or config.get("ranking") != RANKING:
        raise ValueError("Unsupported pipeline version or ranking policy")
    if config.get("annotation_scope") != "non_Y":
        raise ValueError("Mismatch annotations must be explicitly scoped to non-Y sites")
    if set(config["filters"]) != set(FILTER_KEYS):
        raise ValueError("Unrecognised or missing filter setting")
    filters = {k: number(v, k, integer=True) for k, v in config["filters"].items()}
    paths = {k: (config_path.parent / v).resolve() for k, v in config["inputs"].items()}
    candidates = load_candidates(paths["candidates"])
    annotations = load_annotations(paths["non_y_mismatch_annotations"])
    labels = {r["guide"]: r for r in read_table(paths["guide_labels"], ["guide", "guide_id", "region"])}
    historical = read_table(paths["historical_shortlist_for_comparison_only"], ["guide"])
    guide_set = {r["guide"] for r in candidates}
    if annotations.keys() - guide_set:
        raise ValueError("Mismatch annotations include guides absent from the candidate table")
    if labels.keys() - guide_set:
        raise ValueError("Guide labels include spacers absent from the candidate table")
    decisions = evaluate(candidates, annotations, labels, filters)
    counts = stage_counts(decisions, filters)
    counts_by_stage = {r["stage"]: r["guides"] for r in counts}
    partial = counts_by_stage["required_mismatch_counts_missing"] > 0
    created = datetime.now(timezone.utc)
    destination = Path(output).resolve() if output else (config_path.parent / config["output_root"] / created.strftime("%Y%m%dT%H%M%S_%fZ"))
    if destination.exists():
        raise ValueError(f"Output folder already exists; choose a new run folder: {destination}")
    destination.mkdir(parents=True)
    retained = [r for r in decisions if r["decision"] == "retained"]
    unassessed = [r for r in decisions if r["decision"].startswith("unassessed")]
    rejected = [r for r in decisions if r["decision"].startswith("rejected")]
    for name, subset in [("all_guide_decisions", decisions), ("retained_guides", retained),
                         ("unassessed_guides", unassessed), ("rejected_guides", rejected),
                         ("named_guide_decisions", [r for r in decisions if r["guide_id"]])]:
        write_csv(destination / f"{name}.csv", subset, OUTPUT_COLUMNS)
    write_csv(destination / "filtering_counts.csv", counts, ["stage", "kind", "guides", "description"])
    with (destination / "unassessed_guides.fasta").open("w", encoding="ascii") as handle:
        for row in unassessed:
            handle.write(f">{row['guide']}\n{row['guide']}\n")
    historical_set = {r["guide"] for r in historical}
    retained_set = {r["guide"] for r in retained}
    comparison = [{"guide": g, "in_new_retained_set": g in retained_set, "in_historical_top_hits": g in historical_set}
                  for g in sorted(retained_set | historical_set)]
    write_csv(destination / "historical_comparison.csv", comparison, ["guide", "in_new_retained_set", "in_historical_top_hits"])
    snapshot = json.dumps(config, indent=2) + "\n"
    (destination / "config_snapshot.json").write_text(snapshot, encoding="utf-8")
    summary = {
        "pipeline_version": VERSION, "created_utc": created.isoformat(),
        "status": "partial_mismatch_coverage" if partial else "complete_mismatch_coverage",
        "partial_run_explicitly_allowed": allow_partial,
        "counts": {"input_guides": len(decisions), "exact_match_eligible": counts_by_stage["exact_y_abundance_filter"],
                   "eligible_with_required_mismatch_counts": counts_by_stage["required_mismatch_counts_available"],
                   "eligible_missing_required_mismatch_counts": counts_by_stage["required_mismatch_counts_missing"],
                   "retained": len(retained), "unassessed": len(unassessed), "rejected": len(rejected)},
        "historical_comparison": {"shared": len(retained_set & historical_set),
                                  "new_only": len(retained_set - historical_set),
                                  "historical_only": len(historical_set - retained_set)},
        "historical_comparison_interpretation": "Comparison only: available annotations derive from the historical shortlist, so agreement is not independent reproduction of historical discovery.",
        "filters": filters, "ranking": RANKING, "config_path": str(config_path),
        "config_sha256": sha256(config_path), "pipeline_sha256": sha256(__file__),
        "python_version": sys.version,
        "inputs": {k: {"path": str(p), "sha256": sha256(p), "bytes": p.stat().st_size} for k, p in paths.items()},
        "annotation_provenance": config["annotation_provenance"],
        "score_interpretation": config["score_interpretation"],
    }
    exact_rule = ("no exact matches to non-Y chromosomes" if filters["maximum_exact_non_y_matches"] == 0
                  else f"at most {filters['maximum_exact_non_y_matches']} exact matches to non-Y chromosomes")
    mismatch_rule = ("no non-Y sites with one or two mismatches"
                     if filters["maximum_non_y_1_mismatch_sites"] == filters["maximum_non_y_2_mismatch_sites"] == 0
                     else f"at most {filters['maximum_non_y_1_mismatch_sites']} non-Y one-mismatch sites and {filters['maximum_non_y_2_mismatch_sites']} non-Y two-mismatch sites")
    methods = f"""# Guide filtering from saved reference-match tables

## Procedure for this reanalysis

Candidate 20-nt guide sequences were filtered using saved reference-match counts, combining matches in both orientations. Guides were required to have at least {filters['minimum_exact_y_matches']} exact Y-chromosome matches and {exact_rule}. Among guides with available mismatch annotations, we retained those with {mismatch_rule}. Guides lacking a required mismatch count were classified as unassessed unless an observed criterion established rejection. Three- and four-mismatch counts and saved CFD summaries were retained as annotations. Retained guides were ranked by decreasing exact Y-match count, with alphabetical spacer order used to resolve ties. The same filters were applied to all candidates, including the named historical guides.

## Results and coverage of this run

The input contained **{len(decisions):,} guides**. Exact-match specificity and Y-abundance filters retained **{counts_by_stage['exact_y_abundance_filter']:,}**. Required mismatch counts were available for **{counts_by_stage['required_mismatch_counts_available']:,}** of these guides and missing for **{counts_by_stage['required_mismatch_counts_missing']:,}**. The final decisions were **{len(retained):,} retained**, **{len(unassessed):,} unassessed** and **{len(rejected):,} rejected** across the input table.

Run status: **{summary['status']}**. The number of guides with available annotations is a data-coverage count, not an additional biological filter. Historical shortlist overlap was {len(retained_set & historical_set):,} guides. Available annotations originate from that historical shortlist; overlap therefore does not independently reproduce the original selection from the full candidate universe.

## Scope for thesis writing

This text describes the new table-based reanalysis. It should not be presented as a recovered record of the original experimental guide-selection procedure. No genome search or new activity prediction was performed. {config['score_interpretation']}

Reference context: {config['reference_description']}
"""
    (destination / "methods_and_run_summary.md").write_text(methods, encoding="utf-8")
    summary["output_sha256"] = {p.name: sha256(p) for p in sorted(destination.iterdir()) if p.is_file()}
    (destination / "run_manifest.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(destination), "status": summary["status"], **summary["counts"]}, indent=2))
    if partial and not allow_partial:
        print("Incomplete mismatch coverage: diagnostic outputs written. Use --allow-partial to explicitly accept an annotated-subset run.", file=sys.stderr)
        return 2
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=PROJECT / "filtering_config.json")
    parser.add_argument("--out", type=Path, help="New output folder; existing folders are never overwritten")
    parser.add_argument("--allow-partial", action="store_true", help="Allow a run with incomplete mismatch coverage; unassessed guides remain separate")
    args = parser.parse_args()
    try:
        return run(args.config, args.out, args.allow_partial)
    except (ValueError, KeyError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
