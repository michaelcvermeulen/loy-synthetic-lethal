"""Supplementary full-label HT1080 ATAC volcano, including LOC-annotated regions.

Same data, test and drawing as main panel c (shared atac_volcano helper); only the
label rule differs: 20 lowest P values and the 5 most extreme log2FC per direction,
LOC genes included, plus the five main-figure loci.
"""
from pathlib import Path
import sys, json
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "figures/chapter3"
DATA = FIG / "source_data/ht1080_atac_accessibility"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from thesis_style import apply_style
from atac_volcano import select_labels, draw_volcano, DIR
apply_style(8)
STEM = "ch03_supp_ht1080_atac_volcano_labelled"

da = pd.read_csv(DATA / "differential_peaks.csv")
loci = pd.read_csv(DATA / "locus_definitions.csv")
lab = select_labels(da, n_p=20, n_ext=5, forced=dict(zip(loci.peak_id, loci.label)), include_loc=True)
lab[["peak_id", "text", "direction", "logFC", "PValue", "FDR", "annotation", "nearest_tss_distance_signed"]].to_csv(
    DATA / "volcano_labels_supp.csv", index=False)

fig = plt.figure(figsize=(183 / 25.4, 130 / 25.4))
ax = fig.add_axes([.09, .11, .88, .74])
n_treat, cut = draw_volcano(ax, da, lab, fontsize=6.4, point_scale=1.5, xlim=(-6.2, 3.6), label_gap=.27,
                            headroom=1.12, lost_split=-1.95, right_lo=4.6)
h = [Line2D([], [], ls="", marker="o", markersize=4, color=DIR["gained"], label="Gained (FDR < 0.05)"),
     Line2D([], [], ls="", marker="o", markersize=4, color=DIR["lost"], label="Lost (FDR < 0.05)"),
     Line2D([], [], ls="", marker="o", markersize=5, markerfacecolor="none", markeredgecolor="#222222",
            label="TREAT |log₂FC| > 0.5")]
fig.legend(handles=h, loc="upper center", bbox_to_anchor=(.5, .995), ncol=3, fontsize=7, columnspacing=1.8)
fig.text(.5, .925, "Labels name the nearest TSS gene; bold labels are promoter-overlapping regions",
         ha="center", fontsize=6.6, color="#555555")

(FIG / "previews").mkdir(exist_ok=True)
for ext in ["pdf", "svg", "png"]:
    target, stage = FIG / f"{STEM}.{ext}", FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600, facecolor="white"); stage.replace(target)
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=180)
plt.close(fig)
json.dump({"labels": len(lab), "loc_labels": int(lab.text.str.startswith("LOC").sum()), "treat": n_treat},
          open(DATA / "volcano_supp_checks.json", "w"), indent=1)
print("written", STEM, len(lab), "labels")
