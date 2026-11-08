"""Supplementary figure: promoter vs distal motif enrichment at closing and opening chromatin (HT1080, ATAC FDR < 0.10).

Sources (ATACseq/SEPT_08_2026_REPROCESS; code/84_homer_all_significant.py and 83_polycomb_matched.R, no refit here):
  motifs/TF_extension/HOMER/primary_TMM_allFDR05_promoter_fdr010/known_motif_all_results.csv
  motifs/TF_extension/HOMER/primary_TMM_allFDR05_distal_fdr010/known_motif_all_results.csv
  integration/thesis_RNA_SEPT11/polycomb_matched/polycomb_matched_tests.csv (and _K1)
Each HOMER arm compares targets with GC/accessibility-matched non-differential regions of the same class
(promoter vs promoter; distal vs distal); BH across 879 motifs x 2 directions within each arm.
Motif rows (fixed rule): the 8 lowest-FDR distal closing archetypes, the 3 lowest-FDR distal opening archetypes,
and the 3 lowest-P promoter closing archetypes not already included; one JASPAR motif per archetype.
Added on request (22 Sep 2026): FOXC1 (MA0032.2) and TEAD1 (MA0090.4), appended to the distal closing group
because both families led the all-region closing analysis. TEAD1 represents the TEAD archetype (AC0201); TEAD2
(MA1121.2, unmapped to an archetype) has a lower distal FDR (1.6e-7 vs 8.6e-6).
"""
from pathlib import Path
import sys, shutil, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "figures/chapter3"
DATA = FIG / "source_data/ht1080_atac_promoter_vs_distal"
RUN = ROOT.parent / "ATACseq/SEPT_08_2026_REPROCESS"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(8)
STEM = "ch03_supp_ht1080_atac_promoter_vs_distal"
LOWER, HIGHER, GREY = "#3C5488", "#E64B35", "#555555"
DIVERGE = LinearSegmentedColormap.from_list("div", [LOWER, "#FFFFFF", HIGHER])
FDR = "BH_FDR_both_directions"

DATA.mkdir(parents=True, exist_ok=True)
srcs = {"HOMER_promoter_fdr010.csv": RUN / "motifs/TF_extension/HOMER/primary_TMM_allFDR05_promoter_fdr010/known_motif_all_results.csv",
        "HOMER_distal_fdr010.csv": RUN / "motifs/TF_extension/HOMER/primary_TMM_allFDR05_distal_fdr010/known_motif_all_results.csv",
        "polycomb_matched_tests_K3.csv": RUN / "integration/thesis_RNA_SEPT11/polycomb_matched/polycomb_matched_tests.csv"}
for k, v in srcs.items():
    shutil.copyfile(v, DATA / k)
pd.DataFrame([dict(file=k, source=str(v), sha256=hashlib.sha256(v.read_bytes()).hexdigest()) for k, v in srcs.items()]
             ).to_csv(DATA / "source_manifest.csv", index=False)
pro = pd.read_csv(DATA / "HOMER_promoter_fdr010.csv"); dis = pd.read_csv(DATA / "HOMER_distal_fdr010.csv")
pc = pd.read_csv(DATA / "polycomb_matched_tests_K3.csv")

# ---- rows ----
def top(d, direction, n, by=FDR, exclude=()):
    s = d[(d.direction == direction) & ~d.archetype.isin(exclude)].sort_values([by, "PValue", "motif_id"])
    return s.drop_duplicates("archetype").head(n)
r1 = top(dis, "lost", 8)
EXTRA = ["MA0032.2", "MA0090.4"]  # FOXC1, TEAD1 (requested)
r1 = pd.concat([r1, dis[(dis.direction == "lost") & dis.motif_id.isin(EXTRA)].set_index("motif_id").loc[EXTRA].reset_index()])
r2 = top(dis, "gained", 3)
r3 = top(pro, "lost", 3, by="PValue", exclude=set(r1.archetype) | set(r2.archetype))
rows = pd.concat([r1.assign(group="Distal closing"), r3.assign(group="Promoter closing"), r2.assign(group="Distal opening")])
rows = rows[["motif_id", "motif_name", "archetype", "group"]].reset_index(drop=True)
COLS = [("Closing\npromoter", pro, "lost"), ("Closing\ndistal", dis, "lost"),
        ("Opening\npromoter", pro, "gained"), ("Opening\ndistal", dis, "gained")]
grid = []
for i, r in rows.iterrows():
    for j, (name, d, direction) in enumerate(COLS):
        x = d[(d.motif_id == r.motif_id) & (d.direction == direction)].iloc[0]
        grid.append(dict(row=i, col=j, motif_id=r.motif_id, motif_name=r.motif_name, set=name.replace("\n", " "),
                         target_total=int(x.target_total), fold_enrichment=x.fold_enrichment, PValue=x.PValue, BH_FDR=x[FDR]))
grid = pd.DataFrame(grid); grid.to_csv(DATA / "panel_a_grid.csv", index=False)

fig = plt.figure(figsize=(150 / 25.4, 132 / 25.4))


def letter(x, y, ch):
    fig.text(x, y, ch, fontsize=12, weight="semibold", va="bottom")


# a  dot grid
a = fig.add_axes([.17, .135, .38, .71])
LIM = 1.6
norm = Normalize(-LIM, LIM)
for g in grid.itertuples():
    fe = g.fold_enrichment
    l2 = np.log2(fe) if np.isfinite(fe) and fe > 0 else -LIM
    col = DIVERGE(norm(np.clip(l2, -LIM, LIM)))
    size = 8 + 9 * min(-np.log10(max(g.BH_FDR, 1e-12)), 10)
    filled = g.BH_FDR < .05
    faint = abs(norm(np.clip(l2, -LIM, LIM)) - .5) < .1
    a.scatter(g.col, g.row, s=size, color=col if filled else "white",
              edgecolors="#333333" if filled else ("#BBBBBB" if faint else col), linewidths=.35 if filled else .9, zorder=3)
n = len(rows)
a.set(xlim=(-.6, 3.6), ylim=(n - .5, -1.6))
a.set_yticks(range(n), [m.upper() for m in rows.motif_name], fontsize=6.3)  # human symbol style, roman
for t_, grp in zip(a.get_yticklabels(), rows.group):
    t_.set_color(HIGHER if "opening" in grp else LOWER)
a.set_xticks([])
for j, (name, d, direction) in enumerate(COLS):
    tot = int(d[d.direction == direction].target_total.iloc[0])
    a.text(j, -1.45, name, ha="center", va="bottom", fontsize=6.4, weight="semibold", linespacing=1.05,
           color=LOWER if direction == "lost" else HIGHER)
    a.text(j, -.85, f"n = {tot:,}", ha="center", va="bottom", fontsize=5.6, color=GREY)
a.axvline(1.5, color="#CCCCCC", lw=.8, zorder=0)
for j in (0.5, 2.5):
    a.axvline(j, color="#EEEEEE", lw=.6, zorder=0)
bounds = rows.groupby("group", sort=False).size().cumsum()
for b_ in bounds.iloc[:-1]:
    a.axhline(b_ - .5, color="#DDDDDD", lw=.6, zorder=0)
for grp, end in bounds.items():
    start = end - (rows.group == grp).sum()
    a.text(3.75, (start + end - 1) / 2, grp.replace(" ", "\n"), fontsize=5.6, color=GREY, va="center", linespacing=1.1)
a.tick_params(axis="y", length=0, pad=3)
for sp in a.spines.values(): sp.set_visible(False)
letter(.02, .92, "a")
cbx = fig.add_axes([.19, .075, .16, .012])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=DIVERGE), cax=cbx, orientation="horizontal")
cb.set_ticks([-1.5, 0, 1.5]); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.4, length=2)
cb.set_label("log₂ fold enrichment vs matched background", fontsize=5.6, labelpad=1)
lh = [Line2D([], [], ls="", marker="o", markersize=np.sqrt(8 + 9 * v), markerfacecolor="#999999", markeredgecolor="#333333",
             markeredgewidth=.35, label=l) for v, l in [(2, "FDR 1e-2"), (6, "1e-6"), (10, "≤ 1e-10")]]
lh.append(Line2D([], [], ls="", marker="o", markersize=4, markerfacecolor="white", markeredgecolor="#999999", label="FDR ≥ 0.05"))
fig.legend(handles=lh, loc="lower left", bbox_to_anchor=(.40, .045), ncol=2, fontsize=5.4, handletextpad=.3,
           columnspacing=.8, borderpad=0, frameon=False)

# b  significant motif counts
b = fig.add_axes([.74, .40, .23, .32])
cnt = []
for name, d in [("Promoter", pro), ("Distal", dis)]:
    for direction, col in [("lost", LOWER), ("gained", HIGHER)]:
        s = d[d.direction == direction]
        cnt.append(dict(set=name, direction=direction, regions=int(s.target_total.iloc[0]), sig=int((s[FDR] < .05).sum())))
cnt = pd.DataFrame(cnt); cnt.to_csv(DATA / "panel_b_counts.csv", index=False)
for k, r in cnt.iterrows():
    x = (0 if r.set == "Promoter" else 1) + (-.18 if r.direction == "lost" else .18)
    b.bar(x, r.sig, width=.34, color=LOWER if r.direction == "lost" else HIGHER)
    b.text(x, r.sig + 2, f"{r.sig}", ha="center", va="bottom", fontsize=5.8)
b.set_xticks([0, 1], ["Promoter", "Distal"], fontsize=6.4)
b.set_ylabel("Motifs at BH FDR < 0.05", fontsize=6.3); b.set_ylim(0, 110); b.tick_params(axis="y", labelsize=6)
b.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=LOWER, label="Closing"), plt.Rectangle((0, 0), 1, 1, color=HIGHER, label="Opening")],
         loc="upper left", fontsize=5.6, handlelength=1, borderpad=0, frameon=False)
letter(.66, .78, "b")

(FIG / "previews").mkdir(exist_ok=True)
for ext in ["pdf", "svg", "png"]:
    target, stage = FIG / f"{STEM}.{ext}", FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600, facecolor="white"); stage.replace(target)
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=200)
plt.close(fig)
print("written", STEM, len(rows), "rows"); print(cnt.to_string(index=False))
