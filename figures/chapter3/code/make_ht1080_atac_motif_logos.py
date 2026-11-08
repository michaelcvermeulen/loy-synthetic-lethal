"""Compact HOMER known-motif table with sequence logos (HT1080 ATAC, closing and opening peaks).

Sources (ATACseq/SEPT_08_2026_REPROCESS; TF extension 2026-09-08, no refit):
  motifs/TF_extension/HOMER/known_motif_all_results.csv (primary_TMM arm)
  config/TF_extension/JASPAR2024_CORE_vertebrates.jaspar (the matrices HOMER was run with)
Rows: the 10 lowest-FDR archetypes enriched in closing (lost) peaks and every motif enriched in opening (gained)
peaks (BH FDR < 0.05 across 879 motifs x 2 directions), one motif per Vierstra archetype.
Logos: information content (bits) from JASPAR counts with a 0.25 pseudocount per base; drawn directly here.
"""
from pathlib import Path
import sys, re, shutil, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "figures/chapter3"
# Optional argument "allFDR05": plot the comparison HOMER arm (all FDR < 0.05 peaks, no fold-change filter;
# code/84_homer_all_significant.py) instead of the primary FDR < 0.05 and |log2FC| > 1 arm.
# Variants (author decisions 22 Sep 2026): "allFDR05" = all FDR < 0.05 peaks (the supplementary figure);
# "primary" = FDR < 0.05 and |log2FC| > 1 (exploratory); "promoter" = all FDR < 0.05 promoter peaks vs promoter
# background (exploratory; no motif passes FDR 0.05, so rows are the lowest-P archetypes).
VARIANT = sys.argv[1] if len(sys.argv) > 1 else "allFDR05"
assert VARIANT in ("primary", "allFDR05", "promoter", "promoter010", "distal010")
ARM_OF = {"primary": "primary_TMM", "allFDR05": "primary_TMM_allFDR05", "promoter": "primary_TMM_allFDR05_promoter",
          "promoter010": "primary_TMM_allFDR05_promoter_fdr010", "distal010": "primary_TMM_allFDR05_distal_fdr010"}  # promoter010 = promoter peaks at target FDR < 0.10
DATA = FIG / {"primary": "source_data/ht1080_atac_motif_logos_FC1", "allFDR05": "source_data/ht1080_atac_motif_logos_allFDR05",
              "promoter": "source_data/ht1080_atac_motif_logos_promoter",
              "promoter010": "source_data/ht1080_atac_motif_logos_promoter_fdr010",
              "distal010": "source_data/ht1080_atac_motif_logos_distal_fdr010"}[VARIANT]
RUN = ROOT.parent / "ATACseq/SEPT_08_2026_REPROCESS"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style, SKILL
apply_style(8)
STEM = {"allFDR05": "ch03_main_ht1080_atac_motif_logos",
        "primary": "ch03_exploratory_ht1080_atac_motif_logos_fold_change_filtered",
        "promoter": "ch03_exploratory_ht1080_atac_motif_logos_promoter",
        "promoter010": "ch03_exploratory_ht1080_atac_motif_logos_promoter_fdr010",
        "distal010": "ch03_exploratory_ht1080_atac_motif_logos_distal_fdr010"}[VARIANT]
LOWER, HIGHER, GREY = "#3C5488", "#E64B35", "#555555"
BASE_COL = {"A": "#1A9850", "C": "#2B6CB0", "G": "#F2A900", "T": "#D7301F"}
FONT = FontProperties(fname=SKILL / "assets/fonts/OpenSans-SemiBold.ttf")
N_LOST = 10

DATA.mkdir(parents=True, exist_ok=True)
HOMER_RES = RUN / ("motifs/TF_extension/HOMER/known_motif_all_results.csv" if VARIANT == "primary"
                   else f"motifs/TF_extension/HOMER/{ARM_OF[VARIANT]}/known_motif_all_results.csv")
srcs = {"known_motif_all_results.csv": HOMER_RES,
        "JASPAR2024_CORE_vertebrates.jaspar": RUN / "config/TF_extension/JASPAR2024_CORE_vertebrates.jaspar"}
for k, v in srcs.items():
    shutil.copyfile(v, DATA / k)
pd.DataFrame([dict(file=k, source=str(v), sha256=hashlib.sha256(v.read_bytes()).hexdigest()) for k, v in srcs.items()]
             ).to_csv(DATA / "source_manifest.csv", index=False)

# ---- JASPAR matrices ----
pfm, cur = {}, None
for line in open(DATA / "JASPAR2024_CORE_vertebrates.jaspar"):
    line = line.strip()
    if line.startswith(">"):
        cur = line[1:].split()[0]; pfm[cur] = {}
    elif line and cur:
        base = line[0]
        pfm[cur][base] = [float(x) for x in re.findall(r"[-\d.eE]+", line[1:])]
def info_matrix(mid):
    m = np.array([pfm[mid][b] for b in "ACGT"]).T + .25
    p = m / m.sum(axis=1, keepdims=True)
    ic = 2 + (p * np.log2(p)).sum(axis=1)
    return p * ic[:, None]

# ---- rows ----
h = pd.read_csv(DATA / "known_motif_all_results.csv")
h = h[h.arm == ARM_OF[VARIANT]].copy()
FDR = "BH_FDR_both_directions"
if VARIANT in ("promoter", "promoter010"):   # (almost) nothing passes FDR 0.05: show the lowest-P archetypes, flagged as not significant
    sig = h.sort_values(["PValue", "motif_id"])
    lost = sig[sig.direction == "lost"].drop_duplicates("archetype").head(N_LOST)
    gained = sig[sig.direction == "gained"].drop_duplicates("archetype").head(4)
else:
    sig = h[h[FDR] < .05].sort_values([FDR, "motif_id"])
    lost = sig[sig.direction == "lost"].drop_duplicates("archetype").head(N_LOST)
    gained = sig[sig.direction == "gained"].drop_duplicates("archetype")
    if VARIANT == "distal010":  # 20 opening motifs, mostly redundant SOX models: show the top 6 archetypes
        gained = gained.head(6)
rows = pd.concat([lost.assign(section="Closing peaks"), gained.assign(section="Opening peaks")]).reset_index(drop=True)
rows[["section", "motif_id", "motif_name", "archetype", "target_motif_count", "target_total", "target_fraction",
      "background_fraction", "fold_enrichment", "PValue", FDR]].to_csv(DATA / "motif_rows.csv", index=False)
assert all(m in pfm for m in rows.motif_id)

# ---- layout ----
ROW_H = .062
n = len(rows); n_sec = rows.section.nunique()
H_MM = 22 + (n + n_sec) * 6.4 + (6 if VARIANT.startswith(("promoter", "distal")) else 0)
fig = plt.figure(figsize=(120 / 25.4, H_MM / 25.4))
top = 1 - 11 / H_MM
row_h = 6.4 / H_MM
X = {"name": .02, "logo": .20, "logo_w": .30, "bar": .56, "bar_w": .20, "p": .80}
fig.text(X["name"], top + .012, "Motif (JASPAR)", fontsize=6.6, weight="semibold")
fig.text(X["logo"], top + .012, "Sequence logo", fontsize=6.6, weight="semibold")
fig.text(X["bar"], top + .012, "Target / background (%)", fontsize=6.6, weight="semibold")
fig.text(X["p"], top + .012, "P / BH FDR", fontsize=6.6, weight="semibold")
max_pct = 100 * max(rows.target_fraction.max(), rows.background_fraction.max())
y = top
for sec, colour in [("Closing peaks", LOWER), ("Opening peaks", HIGHER)]:
    sub = rows[rows.section == sec]
    if sub.empty:
        continue
    y -= row_h
    tot = int(sub.target_total.iloc[0])
    if VARIANT == "allFDR05":
        fig.text(.008, y + row_h * .3, "a" if sec == "Closing peaks" else "b", fontsize=10, weight="semibold")
    fig.text(X["name"] + (.018 if VARIANT == "allFDR05" else 0), y + row_h * .3, f"{sec} (n = {tot} regions)", fontsize=6.6, weight="semibold", color=colour)
    fig.add_artist(plt.Line2D([X["name"], .98], [y + row_h * .15] * 2, color="#DDDDDD", lw=.5, transform=fig.transFigure))
    for r in sub.itertuples():
        y -= row_h
        fig.text(X["name"], y + row_h * .5, r.motif_name.upper(), fontsize=6.6, va="center", color=colour)  # human symbol style, roman (author, 22 Sep 2026)
        fig.text(X["name"], y + row_h * .12, r.motif_id, fontsize=5.2, va="center", color=GREY)
        # logo
        ic = info_matrix(r.motif_id)
        lw = min(X["logo_w"], X["logo_w"] * len(ic) / 20)
        ax = fig.add_axes([X["logo"], y + row_h * .08, lw, row_h * .84])
        ax.set_xlim(0, len(ic)); ax.set_ylim(0, 2); ax.axis("off")
        for i, col in enumerate(ic):
            base_y = 0.0
            for bi in np.argsort(col):
                hgt = col[bi]
                if hgt < 1e-3:
                    continue
                b = "ACGT"[bi]
                tp = TextPath((0, 0), b, size=1, prop=FONT)
                bb = tp.get_extents()
                tr = Affine2D().translate(-bb.x0, -bb.y0).scale(.92 / bb.width, hgt / bb.height).translate(i + .04, base_y)
                ax.add_patch(PathPatch(tr.transform_path(tp), color=BASE_COL[b], lw=0))
                base_y += hgt
        ax.plot([0, 0], [0, 2], color="#BBBBBB", lw=.4)
        # enrichment bars (target vs background %)
        bx = fig.add_axes([X["bar"], y + row_h * .12, X["bar_w"], row_h * .76])
        bx.barh([1], [100 * r.target_fraction], color=colour, height=.8)
        bx.barh([0], [100 * r.background_fraction], color="#BBBBBB", height=.8)
        bx.set_xlim(0, max_pct * 1.3); bx.set_ylim(-.6, 1.6); bx.axis("off")
        bx.text(100 * r.target_fraction + max_pct * .02, 1, f"{100 * r.target_fraction:.0f}", fontsize=5.2, va="center")
        bx.text(100 * r.background_fraction + max_pct * .02, 0, f"{100 * r.background_fraction:.0f}", fontsize=5.2,
                va="center", color=GREY)
        fig.text(X["bar"] + X["bar_w"] + .005, y + row_h * .5, f"{r.fold_enrichment:.1f}×", fontsize=5.8,
                 va="center", weight="semibold", color=colour)
        # P and FDR
        fig.text(X["p"] + .03, y + row_h * .5, f"{r.PValue:.0e}  /  {getattr(r, FDR):.0e}".replace("e-0", "e-"),
                 fontsize=5.8, va="center")
fig.text(X["logo"], y - row_h * .55, "Letter height = information content (0–2 bits)", fontsize=5.2, color=GREY)
if VARIANT == "distal010":
    ns = h[h[FDR] < .05].groupby("direction").size()
    fig.text(X["name"], y - row_h * 1.25, f"Distal peaks (ATAC FDR < 0.10) vs distal background: {ns.get('lost', 0)} closing and "
             f"{ns.get('gained', 0)} opening motifs at BH FDR < 0.05; top archetypes shown", fontsize=5.2, color=GREY, weight="semibold")
if VARIANT.startswith("promoter"):
    n_sig = int((h[FDR] < .05).sum()); thr = "0.10" if VARIANT == "promoter010" else "0.05"
    fig.text(X["name"], y - row_h * 1.25, f"Promoter peaks (ATAC FDR < {thr}) vs promoter background: {n_sig} motif(s) at "
             f"BH FDR < 0.05; other rows are the lowest-P archetypes", fontsize=5.2, color=HIGHER, weight="semibold")
fig.text(X["bar"], y - row_h * .55, "Coloured: target regions; grey: matched background", fontsize=5.2, color=GREY)

(FIG / "previews").mkdir(exist_ok=True)
for ext in ["pdf", "svg", "png"]:
    target, stage = FIG / f"{STEM}.{ext}", FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600, facecolor="white"); stage.replace(target)
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=220)
plt.close(fig)
print("written", STEM, n, "rows;", f"{H_MM:.0f} mm tall")
