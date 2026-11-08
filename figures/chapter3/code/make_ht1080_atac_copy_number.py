"""Main figure: copy-number context of HT1080 differential accessibility.

Source: ATACseq/SEPT_08_2026_REPROCESS/code/85_atac_copy_number_context.py (outputs in
integration/thesis_RNA_SEPT11/copy_number_context/) and the D0 lpWGS hg38 100 kb bins from
72_hg38_WGS_copy_number.R. No biological models are refitted; display-only linear trend lines use the frozen chromosome summaries. Tables are frozen into source_data/ht1080_atac_copy_number/.
a = per-chromosome shifts (ATAC per passage, RNA, D0 WGS) with differential-region rates;
b = arm-level RNA/lpWGS, ATAC/lpWGS, and combined/August/May RNA/ATAC scatters;
c = chr18 and chr13 position tracks; d = matched motif sensitivity;
e = RNA-ATAC promoter concordance in all vs copy-number-neutral promoters.
"""
from pathlib import Path
import sys, shutil, hashlib, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "figures/chapter3"
DATA = FIG / "source_data/ht1080_atac_copy_number"
RUN = ROOT.parent / "ATACseq/SEPT_08_2026_REPROCESS"
CNO = RUN / "integration/thesis_RNA_SEPT11/copy_number_context"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(8)
STEM = "ch03_main_ht1080_atac_copy_number"
LOWER, HIGHER, GREY = "#3C5488", "#E64B35", "#555555"
ATAC_C, RNA_C, WGS_C = "#333333", "#00A087", "#8C6BB1"

DATA.mkdir(parents=True, exist_ok=True)
srcs = {f: CNO / f for f in ["chromosome_shift_by_passage.csv", "motif_CMH_chromosome_exclusion.csv", "gene_CN_context.csv",
                             "summary.json", "region_CN_context.csv"]}
srcs["WGS_hg38_100kb_bins.csv"] = RUN / "integration/RNA_interaction_CN/WGS_hg38_100kb_bins.csv"
for k, v in srcs.items():
    shutil.copyfile(v, DATA / k)
pd.DataFrame([dict(file=k, source=str(v), sha256=hashlib.sha256(v.read_bytes()).hexdigest()) for k, v in srcs.items()]
             ).to_csv(DATA / "source_manifest.csv", index=False)
sh = pd.read_csv(DATA / "chromosome_shift_by_passage.csv")
mot = pd.read_csv(DATA / "motif_CMH_chromosome_exclusion.csv")
gen = pd.read_csv(DATA / "gene_CN_context.csv")
reg = pd.read_csv(DATA / "region_CN_context.csv", usecols=["chrom", "start", "end", "logFC", "call"])
bins = pd.read_csv(DATA / "WGS_hg38_100kb_bins.csv", usecols=["chrom", "start0", "end0", "use", "delta_segment_log2"])
summ = json.load(open(DATA / "summary.json"))

fig = plt.figure(figsize=(183 / 25.4, 215 / 25.4))


def letter(x, y, ch):
    fig.text(x, y, ch, fontsize=12, weight="semibold", va="bottom")


# ---- a: per-chromosome shifts + differential-region rates ----
order = [f"chr{i}" for i in range(1, 23)]
s = sh.set_index("chrom").loc[order]
n_reg = reg.groupby("chrom").size().reindex(order)
x = np.arange(len(order))
a = fig.add_axes([.075, .795, .425, .17])
a.axhline(0, color="#BBBBBB", lw=.6)
for k, (ps, off) in enumerate([("11", -.18), ("12", 0), ("21", .18)]):
    a.scatter(x + off, s[f"ATAC_p{ps}"], s=9, color=ATAC_C, lw=0, zorder=3, label="ATAC (each passage)" if k == 0 else None)
a.scatter(x, s.RNA_median_logFC, s=16, marker="D", facecolor="white", edgecolor=RNA_C, lw=.8, zorder=4, label="RNA")
a.hlines(s.WGS_D0_median_segment_delta, x - .38, x + .38, color=WGS_C, lw=1.6, zorder=2, label="D0 lpWGS")
a.set_xlim(-.6, len(order) - .4); a.set_ylim(-.42, .48)
a.set_xticks(x, [c.replace("chr", "") for c in order], fontsize=6.3)
a.set_ylabel("Y-depleted − sgCTRL\n(log2, chromosome median)", fontsize=6.5)
a.tick_params(axis="y", labelsize=6)
for c in ("chr13", "chr18"):
    a.axvspan(order.index(c) - .5, order.index(c) + .5, color="#F2F2F2", zorder=0)
a.legend(loc="upper left", ncol=3, fontsize=5.8, frameon=False, handletextpad=.3, columnspacing=1, borderpad=0)
cs = summ["chromosome_shift_spearman"]
a.spines["bottom"].set_visible(False); a.tick_params(axis="x", length=0, labelbottom=False)
# rate strip
r_ = fig.add_axes([.075, .695, .425, .075], sharex=a)
w = .36
r_.bar(x - w / 2, 100 * s.lost_FDR05 / n_reg, width=w, color=LOWER, lw=0, label="Lost")
r_.bar(x + w / 2, 100 * s.gained_FDR05 / n_reg, width=w, color=HIGHER, lw=0, label="Gained")
for c in ("chr13", "chr18"):
    r_.axvspan(order.index(c) - .5, order.index(c) + .5, color="#F2F2F2", zorder=0)
i18 = order.index("chr18")
r_.text(i18 + w / 2, 100 * s.loc["chr18", "gained_FDR05"] / n_reg["chr18"] + .15,
        f"{int(s.loc['chr18', 'gained_FDR05'])} of {summ['chr18']['gained_total']} gained", fontsize=5.6, ha="center", va="bottom", color=HIGHER)
r_.set_ylim(0, 2.9)
r_.set_ylabel("% regions\nFDR < 0.05", fontsize=6.2); r_.tick_params(labelsize=6)
r_.set_xticks(x, [c.replace("chr", "") for c in order], fontsize=6.3); r_.set_xlabel("Chromosome", fontsize=6.5)
r_.legend(loc="upper right", ncol=2, fontsize=5.6, frameon=False, handlelength=.9, borderpad=0, bbox_to_anchor=(1, 1.05))
letter(.01, .955, "a")

letter(.525, .955, "b")

# Five arm comparisons: DNA comparisons above, joint and period-specific RNA below.
from scipy.stats import spearmanr, pearsonr
from adjustText import adjust_text
scatter_data=pd.read_csv(FIG/'source_data/ht1080_rna_period_atac_arms/arm_values.csv').set_index('arm')
plots=[
 ([.58,.735,.16,.22],'WGS_D0_median_segment_delta','RNA_median_logFC','lpWGS','RNA',RNA_C),
 ([.81,.735,.16,.22],'WGS_D0_median_segment_delta','ATAC_mean_shift','lpWGS','ATAC',ATAC_C),
 ([.075,.47,.255,.155],'RNA_median_logFC','ATAC_mean_shift','RNA, combined','ATAC','#8C6BB1'),
 ([.395,.47,.255,.155],'August_RNA','ATAC_mean_shift','RNA, August','ATAC','#8C6BB1'),
 ([.715,.47,.255,.155],'May_RNA','ATAC_mean_shift','RNA, May','ATAC','#8C6BB1')]
for rect,xfield,field,xlabel,ylabel,color in plots:
 ax=fig.add_axes(rect);xx=scatter_data[xfield];yy=scatter_data[field]
 check=spearmanr(xx,yy);pr=pearsonr(xx,yy)
 ax.axhline(0,color='#DDDDDD',lw=.5);ax.axvline(0,color='#DDDDDD',lw=.5)
 ax.scatter(xx,yy,s=10,color=color,edgecolor='white',lw=.25,zorder=3)
 slope,intercept=np.polyfit(xx,yy,1);lx=np.array([xx.min(),xx.max()]);ax.plot(lx,intercept+slope*lx,color=color,lw=.7)
 ax.set_xlim(-.38,.55);ax.set_ylim(-.29,.58)
 ax.tick_params(labelsize=6,pad=2);ax.set_xticks([-.3,0,.3]);ax.set_yticks([-.2,0,.2,.4])
 ax.set_xlabel(xlabel+' difference (log₂)',fontsize=6.5,labelpad=2)
 ax.set_ylabel(ylabel+' difference (log₂)',fontsize=6.5,labelpad=2)
 pf=f'{check.pvalue:.1e}' if check.pvalue<.001 else f'{check.pvalue:.3f}'
 txt=f'ρ = {check.statistic:.2f}; P = {pf}'
 if xfield!='WGS_D0_median_segment_delta':
  txt+=f'\nr = {pr.statistic:.2f}; P = {pr.pvalue:.3f}'
 ax.text(.03,.98,txt,transform=ax.transAxes,fontsize=5.7,va='top')
 texts=[]
 for arm in ['chr13q','chr18p','chr18q','chr20p']:
  if arm=='chr20p' and xfield=='WGS_D0_median_segment_delta':continue
  row=scatter_data.loc[arm]
  texts.append(ax.text(row[xfield],row[field],arm.replace('chr',''),fontsize=5.5))
 adjust_text(texts,x=xx.to_numpy(),y=yy.to_numpy(),ax=ax,arrowprops=dict(arrowstyle='-',color='.5',lw=.4))
fig.text(.54,.674,'39 autosomal arms; ρ, Spearman; r, Pearson',fontsize=6)

# ---- c: chr18 and chr13 position tracks ----
def track(ax, chrom, title):
    d = reg[reg.chrom == chrom].sort_values("start")
    mid = (d.start + d.end) / 2e6
    ax.scatter(mid, d.logFC.clip(-3, 3), s=1.2, color="#C8C8C8", lw=0, rasterized=True, zorder=1)
    edges = np.arange(0, mid.max() + 5, 5)
    cat = pd.cut(mid, edges)
    med = d.logFC.groupby(cat, observed=False).median()
    ctr = (edges[:-1] + edges[1:]) / 2
    ok = d.groupby(cat, observed=False).size().to_numpy() >= 20
    ax.plot(ctr[ok], med.to_numpy()[ok], color=ATAC_C, lw=1.1, zorder=3, label="ATAC log2FC, 5 Mb median")
    b = bins[(bins.chrom == chrom) & bins.use & bins.delta_segment_log2.notna()]
    ax.scatter((b.start0 + b.end0) / 2e6, b.delta_segment_log2, s=1.5, color=WGS_C, lw=0, zorder=2, rasterized=True,
               label="D0 lpWGS segment Δ")
    for call, col, y in [("lost", LOWER, -1.9), ("gained", HIGHER, 1.9)]:
        q = d[d.call == call]
        ax.vlines((q.start + q.end) / 2e6, y - .18, y + .18, color=col, lw=.6, zorder=3)
    ax.axhline(0, color="#999999", lw=.5)
    ax.set_ylim(-2.2, 2.2); ax.set_xlim(0, mid.max() + 1)
    ax.set_title(title, fontsize=6.8, loc="left", pad=2)
    ax.tick_params(labelsize=6)


b1 = fig.add_axes([.075, .245, .405, .14]); b2 = fig.add_axes([.075, .05, .405, .14])
track(b1, "chr18", f"chr18  (gained {summ['chr18']['gained_FDR05']}; lost 4)")
track(b2, "chr13", "chr13  (gained 1; lost 36)")
b1.set_ylabel("log2 (Y-depleted / sgCTRL)", fontsize=6.3); b2.set_ylabel("log2 (Y-depleted / sgCTRL)", fontsize=6.3)
b2.set_xlabel("Position (Mb, hg38)", fontsize=6.5)
b1.legend(handles=[Line2D([], [], color=ATAC_C, lw=1.1, label="ATAC, 5 Mb median"),
                   Line2D([], [], ls="", marker="o", ms=2.5, color=WGS_C, label="lpWGS Δ"),
                   Line2D([], [], ls="", marker="o", ms=2.5, color="#C8C8C8", label="Region"),
                   Line2D([], [], color=HIGHER, lw=.8, label="Gained / lost")],
          loc="lower left", bbox_to_anchor=(0, 1.11), fontsize=5.3, frameon=False, ncol=4, handlelength=1.2, columnspacing=.8, borderpad=0)
letter(.01, .405, "c")

# ---- d: motif CMH, all vs excluding chr13/chr18 targets ----
c = fig.add_axes([.62, .245, .35, .15])
fams = ["TFAP2", "ZIC", "TEAD", "GATA", "RUNX", "EBF", "SOX6"]
for i, m in enumerate(fams):
    for dy, lab, mk in [(-.16, "all targets", "o"), (.16, "excluding", "s")]:
        r = mot[(mot.motif == m) & mot.target_set.str.startswith(lab)].iloc[0]
        col = HIGHER if r.direction == "gained" else LOWER
        c.plot([r.CI_low, r.CI_high], [i + dy] * 2, color=col, lw=.8)
        c.scatter(r.MH_OR, i + dy, s=13, marker=mk, color=col if mk == "o" else "white", edgecolor=col, lw=.8, zorder=3)
c.axvline(1, color="#888888", lw=.6, ls=(0, (4, 3)))
c.set_xscale("log"); c.set_xlim(.8, 10); c.set_ylim(len(fams) - .5, -.7)
c.set_xticks([1, 2, 4, 8], ["1", "2", "4", "8"], fontsize=6); c.minorticks_off()
c.set_yticks(range(len(fams)), fams, fontsize=6.3)
for t_, m in zip(c.get_yticklabels(), fams):
    t_.set_color(HIGHER if m == "SOX6" else LOWER)
c.tick_params(axis="y", length=0); c.spines["left"].set_visible(False)
c.set_xlabel("Mantel–Haenszel OR vs matched controls (95% CI)", fontsize=6.2)
c.legend(handles=[Line2D([], [], ls="", marker="o", ms=3.5, color=GREY, label="All targets"),
                  Line2D([], [], ls="", marker="s", ms=3.5, mfc="white", mec=GREY, label="Chromosome exclusion")],
         loc="lower center", bbox_to_anchor=(.45, 1.0), ncol=2, fontsize=5.6, frameon=False, handletextpad=.2, borderpad=0)
letter(.55, .405, "d")

# ---- e: concordance in CN-neutral promoters ----
dax = fig.add_axes([.62, .05, .35, .13])
sens = summ["lost_peak_rate_sensitivity"]
bars = [("Lower-expression\nDE genes", sens["all"]["pct"], sens["CN_neutral_0.20"]["pct"],
         (sens["all"]["with_lost_peak_010"], sens["all"]["DE_lower"]), (sens["CN_neutral_0.20"]["with_lost_peak_010"], sens["CN_neutral_0.20"]["DE_lower"])),
        ("Non-DE genes", sens["nonDE_all"]["pct"], sens["nonDE_CN_neutral_0.20"]["pct"],
         (sens["nonDE_all"]["with_lost_peak_010"], sens["nonDE_all"]["genes"]), (sens["nonDE_CN_neutral_0.20"]["with_lost_peak_010"], sens["nonDE_CN_neutral_0.20"]["genes"]))]
out = []
for i, (lab, p_all, p_neu, n_all, n_neu) in enumerate(bars):
    dax.barh(i - .19, p_all, height=.36, color=LOWER, lw=0)
    dax.barh(i + .19, p_neu, height=.36, color="white", edgecolor=LOWER, lw=.8)
    dax.text(p_all + .6, i - .19, f"{n_all[0]}/{n_all[1]:,}", va="center", fontsize=5.6)
    dax.text(p_neu + .6, i + .19, f"{n_neu[0]}/{n_neu[1]:,}", va="center", fontsize=5.6)
    out += [dict(genes=lab.replace("\n", " "), promoters="all", pct=p_all, n=n_all[0], total=n_all[1]),
            dict(genes=lab.replace("\n", " "), promoters="CN-neutral", pct=p_neu, n=n_neu[0], total=n_neu[1])]
pd.DataFrame(out).to_csv(DATA / "panel_d_concordance.csv", index=False)
dax.set_yticks([0, 1], [b[0] for b in bars], fontsize=6.2); dax.set_ylim(1.6, -.6)
dax.set_xlim(0, 36); dax.tick_params(axis="x", labelsize=6); dax.tick_params(axis="y", length=0)
dax.set_xlabel("% with a lost promoter peak (ATAC FDR < 0.10)", fontsize=6.2)
dax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=LOWER, label="All promoters"),
                    plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=LOWER, label="lpWGS |Δ| < 0.2")],
           loc="lower right", fontsize=5.6, frameon=False, handlelength=1, borderpad=0)
letter(.55, .19, "e")

(FIG / "previews").mkdir(exist_ok=True)
for ext in ["pdf", "svg", "png"]:
    target, stage = FIG / f"{STEM}.{ext}", FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600, facecolor="white"); stage.replace(target)
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=200)
plt.close(fig)
print("written", STEM)
