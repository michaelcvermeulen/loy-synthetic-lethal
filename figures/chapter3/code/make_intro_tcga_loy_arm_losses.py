"""Recreate Qi et al. Figure 2E from the frozen source workbooks.

Run from any directory with numpy, pandas, scipy, openpyxl and matplotlib
available. No downloads or changes to the original inputs.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox, offset_copy

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "figures/chapter3"
SOURCE = OUT / "source_data/qi_loy_arm_aneuploidy"
PREVALENCE = OUT / "source_data/qi_tcga_loy_prevalence"
DATA = OUT / "source_data/intro_tcga_loy_arm_losses"
DATA.mkdir(parents=True, exist_ok=True)
STEM = "ch03_main_intro_tcga_loy_arm_losses"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


# Rebuild prevalence from individual calls and verify the published summaries.
cases = pd.read_excel(PREVALENCE / "Qi_S1_Table.xlsx", sheet_name="Table S1A")
published = pd.read_excel(PREVALENCE / "Qi_S1_Table.xlsx", sheet_name="Table S1C")
assert cases.case_id.is_unique and cases.Gender.eq("male").all()
called = cases[cases.Y_status.notna()]
assert len(called) == 5014 and called.Y_status.eq("LOY").sum() == 1390
counts = called.groupby("Cohort").agg(
    n_y_status=("case_id", "size"),
    n_complete_loy=("Y_status", lambda s: s.eq("LOY").sum()),
)
types = published.loc[
    published.Male_sample_count.notna() & published.Fraction_LOY.notna(),
    ["Tumor_type", "Male_sample_count", "Fraction_LOY"],
]
table = types.merge(counts, left_on="Tumor_type", right_index=True,
                    validate="one_to_one")
assert len(table) == 28 and table.n_y_status.sum() == 5003
assert np.array_equal(table.Male_sample_count, table.n_y_status)
table["complete_loy_percent"] = 100 * table.n_complete_loy / table.n_y_status
np.testing.assert_allclose(table.complete_loy_percent / 100, table.Fraction_LOY,
                           atol=1e-12, rtol=0)

# Match the authors' notebook: use all male arm-score sample records, including
# repeat records from a case and cases without a Qi Y-status call.
clinical = pd.read_excel(SOURCE / "TCGA-CDR-SupplementalTableS1.xlsx",
                         sheet_name="TCGA-CDR",
                         usecols=["bcr_patient_barcode", "type", "gender"])
assert clinical.bcr_patient_barcode.is_unique
arms = pd.read_excel(SOURCE / "Taylor_2018_Table_S2.xlsx", skiprows=1)
arms["case_id"] = arms.Sample.str[:12]
arms = arms.merge(clinical, left_on="case_id", right_on="bcr_patient_barcode",
                  how="left", validate="many_to_one", indicator=True)
matched = arms._merge.eq("both")
assert arms.loc[matched, "Type"].eq(arms.loc[matched, "type"]).all()
male = arms.loc[arms.gender.eq("MALE")].copy()
assert len(male) == 5039 and male.case_id.nunique() == 5005
assert male.AS_del.notna().all()
np.testing.assert_allclose(male["AneuploidyScore(AS)"], male.AS_del + male.AS_amp)
means = male.groupby("type").agg(
    mean_autosomal_arm_losses=("AS_del", "mean"),
    n_arm_records=("AS_del", "size"), n_arm_cases=("case_id", "nunique"),
)
table = table.merge(means, left_on="Tumor_type", right_index=True,
                    validate="one_to_one")
names = pd.read_csv(PREVALENCE / "prevalence_by_tumour_type.csv")
table = table.merge(names[["Tumor_type", "tumour_name"]], on="Tumor_type",
                    validate="one_to_one")
table = table.sort_values("Tumor_type").reset_index(drop=True)
assert len(table) == 28 and table.n_arm_records.sum() == 5028
assert table.mean_autosomal_arm_losses.notna().all()

# Independent comparison with the previous source reconstruction.
earlier = pd.read_csv(SOURCE / "prevalence_vs_mean_arm_loss.csv").set_index("Tumor_type")
np.testing.assert_allclose(table.complete_loy_percent,
                           earlier.loc[table.Tumor_type, "percentage_complete_loy"], atol=1e-12)
np.testing.assert_allclose(table.mean_autosomal_arm_losses,
                           earlier.loc[table.Tumor_type, "mean_arm_losses"], atol=1e-12)
r, p = pearsonr(table.complete_loy_percent, table.mean_autosomal_arm_losses)
assert round(float(r), 2) == 0.33 and round(float(p), 4) == 0.0833

columns = ["Tumor_type", "tumour_name", "complete_loy_percent", "n_complete_loy",
           "n_y_status", "mean_autosomal_arm_losses", "n_arm_records", "n_arm_cases"]
table[columns].to_csv(DATA / "plotted_values.csv", index=False)
plotted_records = male[male.type.isin(table.Tumor_type)]
plotted_records[["Sample", "case_id", "type", "AS_del"]].to_csv(
    DATA / "male_arm_loss_records.csv", index=False)

font = apply_style(9)
focus = {"KIRP", "KIRC", "UVM"}
red, grey = "#E64B35", "#737B84"
fig, ax = plt.subplots(figsize=(160 / 25.4, 135 / 25.4), dpi=180)
fig.subplots_adjust(left=.135, right=.985, bottom=.12, top=.985)
ax.set(xlim=(-2, 89), ylim=(-.4, 14.0),
       xlabel="Male tumours with complete Y-chromosome loss (%)",
       ylabel="Mean number of autosomal chromosome arms lost")
ax.set_xticks([0, 20, 40, 60, 80])
ax.set_yticks([0, 2, 4, 6, 8, 10, 12])
ax.tick_params(length=3, pad=3)
ax.xaxis.labelpad = 7
ax.yaxis.labelpad = 7
texts = []
# Label offsets in physical points; only labels move, never data coordinates.
offsets = {
    "PCPG": (35, -14, "left"), "GBM": (8, 4, "left"),
    "THYM": (-4, 30, "left"), "PRAD": (23, 5, "left"),
    "THCA": (6, -2, "left"), "LAML": (6, 3, "left"),
    "LGG": (6, 5, "left"), "SARC": (-8, -13, "right"),
    "LIHC": (6, 6, "left"), "COAD": (6, -3, "left"),
    "SKCM": (6, 5, "left"), "HNSC": (6, -4, "left"),
    "CHOL": (0, 14, "left"), "STAD": (10, 6, "left"),
    "PAAD": (8, -8, "left"),
}
for row in table.itertuples():
    marked = row.Tumor_type in focus
    colour = red if marked else grey
    ax.scatter(row.complete_loy_percent, row.mean_autosomal_arm_losses,
               s=34 if marked else 27, c=colour, edgecolors="white",
               linewidths=.6, zorder=3)
    dx, dy, ha = offsets.get(row.Tumor_type, (6, 4, "left"))
    texts.append(ax.text(row.complete_loy_percent,
                         row.mean_autosomal_arm_losses, row.Tumor_type,
                         transform=offset_copy(ax.transData, fig=fig, x=dx, y=dy, units="points"),
                         ha=ha, va="center",
                         fontsize=8.5, color=red if marked else "#3F454B",
                         weight="semibold" if marked else "normal", zorder=4))
stat = ax.text(.985, .98, f"Pearson r = {r:.3f}; P = {p:.4f}\nn = 28 cancer types",
               transform=ax.transAxes, va="top", ha="right", fontsize=9)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
boxes = [t.get_window_extent(renderer) for t in texts]
for row, box in zip(table.itertuples(), boxes):
    point = ax.transData.transform((row.complete_loy_percent, row.mean_autosomal_arm_losses))
    edge = np.array([np.clip(point[0], box.x0, box.x1), np.clip(point[1], box.y0, box.y1)])
    distance = np.linalg.norm(edge - point)
    if distance > 10 * fig.dpi / 72:
        direction = (edge - point) / distance
        ends = np.array([point + direction * 4 * fig.dpi / 72,
                         edge - direction * 2 * fig.dpi / 72])
        line = ax.transData.inverted().transform(ends)
        ax.plot(line[:, 0], line[:, 1], color="#9A9FA5", lw=.45, zorder=2)
overlaps = [(texts[i].get_text(), texts[j].get_text())
            for i in range(len(texts)) for j in range(i + 1, len(texts))
            if boxes[i].overlaps(boxes[j])]
assert not overlaps, f"Overlapping cancer labels: {overlaps}"
radius = 3.5 * fig.dpi / 72
point_overlaps = []
for row in table.itertuples():
    px, py = ax.transData.transform((row.complete_loy_percent, row.mean_autosomal_arm_losses))
    point_box = Bbox.from_bounds(px - radius, py - radius, 2 * radius, 2 * radius)
    point_overlaps.extend((row.Tumor_type, text.get_text())
                          for text, box in zip(texts, boxes) if point_box.overlaps(box))
assert not point_overlaps, f"Labels obscure data points: {point_overlaps}"
assert not any(b.overlaps(stat.get_window_extent(renderer)) for b in boxes)
figure_box = fig.bbox
assert all(figure_box.contains(b.x0, b.y0) and figure_box.contains(b.x1, b.y1)
           for b in boxes)
pd.DataFrame([{"Tumor_type": t.get_text(),
               "offset_x_points": offsets.get(t.get_text(), (6, 4, "left"))[0],
               "offset_y_points": offsets.get(t.get_text(), (6, 4, "left"))[1],
               "horizontal_alignment": offsets.get(t.get_text(), (6, 4, "left"))[2]} for t in texts]).to_csv(
                   DATA / "label_positions.csv", index=False)
for ext in ("pdf", "svg", "png"):
    fig.savefig(OUT / f"{STEM}.{ext}", dpi=600)
fig.savefig(OUT / "previews" / f"{STEM}.png", dpi=180)
plt.close(fig)

source_paths = [PREVALENCE / "Qi_S1_Table.xlsx", SOURCE / "Taylor_2018_Table_S2.xlsx",
                SOURCE / "TCGA-CDR-SupplementalTableS1.xlsx",
                SOURCE / "LOY_LOX WEX Final analysis.ipynb",
                SOURCE / "prevalence_vs_mean_arm_loss.csv"]
meta = {
    "figure": STEM, "created_for": "Chapter 3, shorter 3.1 Introduction draft",
    "source_figure": "Qi et al. (2023), Cell, Figure 2E",
    "doi": "https://doi.org/10.1016/j.cell.2023.06.006",
    "notebook_url": "https://github.com/rheinbaylab/LOY_in_primary_tumors/blob/master/Notebooks/LOY_LOX%20WEX%20Final%20analysis.ipynb",
    "notebook_cells_zero_based": [41, 43],
    "pearson_r": float(r), "pearson_p_two_sided": float(p), "n_cancer_types": 28,
    "n_y_status_cases_plotted": int(table.n_y_status.sum()),
    "n_arm_loss_records_plotted": len(plotted_records),
    "n_arm_loss_unique_cases_plotted": int(plotted_records.case_id.nunique()),
    "unmatched_arm_records_excluded_before_sex_selection": int((~matched).sum()),
    "mean_population": "All available male Taylor sample records per cancer type, as in Qi notebook; repeat case records retained; not restricted to Qi-called cases",
    "LOY_definition": "Complete genomic LOY only; partial and relative losses excluded from numerator",
    "arm_endpoint": "Taylor AS_del, autosomal arm losses only; gains and Y excluded",
    "highlighted_types": sorted(focus),
    "highlight_reason": "Cancer types discussed in Introduction; no significance classification",
    "restyling": "Original robust regression line/confidence band omitted; no new trend or selection model",
    "font": font, "size_mm": [160, 135], "raster_dpi": 600,
    "label_overlap_check": "passed", "missing_value_imputation": "none",
    "versions": {"numpy": np.__version__, "pandas": pd.__version__,
                 "matplotlib": matplotlib.__version__},
    "sources": {str(p.relative_to(ROOT)): sha256(p) for p in source_paths},
    "script_sha256": sha256(Path(__file__)),
    "exports": {ext: sha256(OUT / f"{STEM}.{ext}") for ext in ("pdf", "svg", "png")},
}
(DATA / "provenance.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
print(json.dumps({k: meta[k] for k in ["figure", "pearson_r", "pearson_p_two_sided",
      "n_y_status_cases_plotted", "n_arm_loss_records_plotted",
      "n_arm_loss_unique_cases_plotted", "label_overlap_check"]}, indent=2))
