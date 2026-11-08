"""TSS-centred accessibility at promoters of RNA-defined gene groups, with matched controls.

Companion to the region-centred S08 heatmaps. Genes: the 13,300-gene promoter
integration universe (ATACseq/SEPT_08_2026_REPROCESS/integration/thesis_RNA_SEPT11,
RNA_ATAC_promoter_genes.csv). Coverage: tss_transcript_matrix.npz (extract_ht1080_atac_transcript_tss_matrix.py;
knownGene transcript TSSs, +/-3 kb, 50 bp bins, strand-oriented, fragment CPM of non-Y fragments).

One TSS per gene, chosen at extraction: the gene's annotated transcript TSS with the highest central
(+/-100 bp) coverage pooled over all six libraries, so the choice does not depend on condition.
(The pipeline's gene-level TSS, the most upstream transcript start, often lies upstream of the active
promoter and was not used.)
Target groups:
  down, up        thesis RNA classes (independent of ATAC calls);
  down_no_lost_peak  down genes without a called lost promoter peak (is closing wider than the calls?);
  promoter_lost   genes with >=1 lost promoter peak (primary FDR < 0.05); selected on ATAC,
                  shown only as a positive reference (circular by construction).
Controls: not-DE genes with no differential promoter peak, K = 3 per target, sampled without
replacement within deciles of pooled central TSS coverage (log scale); drawn independently per
group (seed 20260923 + group index). 200 further draws per group give the draw-to-draw spread
(tss_groups_control_draws.csv).
Centre statistic: mean over the central 200 bp (4 bins). Group ratio = mean Y-depleted /
mean sgCTRL (ratio of group means, as in S08). Per-gene log2 ratio uses pseudocount 0.1.
Test: two-sided Mann-Whitney U of per-gene log2 ratio, target vs its controls; BH over 4 tests.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT.parent / "ATACseq/SEPT_08_2026_REPROCESS"
DATA = ROOT / "figures/chapter3/source_data/ht1080_atac_accessibility"
SAMPLES = ["CTRL_11", "CTRL_12", "CTRL_21", "LOY_11", "LOY_12", "LOY_21"]
K, SEED, PSEUDO, NDEC, NDRAW = 3, 20260923, .1, 10, 200

z = np.load(DATA / "tss_transcript_matrix.npz")
nb = z[SAMPLES[0]].shape[1]
C = slice(nb // 2 - 2, nb // 2 + 2)
sites = pd.DataFrame([s.split(":") for s in z["sites"]], columns=["chrom", "tss", "strand", "entrez_id"])
sites["row"] = np.arange(len(sites))
centre = {s: z[s][:, C].mean(1) for s in SAMPLES}
sites["pooled_centre"] = np.mean([centre[s] for s in SAMPLES], axis=0)
assert not sites.entrez_id.duplicated().any()

g = pd.read_csv(RUN / "integration/thesis_RNA_SEPT11/RNA_ATAC_promoter_genes.csv", dtype={"entrez_id": str})
assert len(g) == 13300
g = g.merge(sites, on="entrez_id", how="inner", suffixes=("", "_tss"))
n_no_tss = 13300 - len(g)
for s in SAMPLES:
    g[s] = centre[s][g.row]
g["ctrl_centre"] = g[SAMPLES[:3]].mean(1)
g["loy_centre"] = g[SAMPLES[3:]].mean(1)
g["log2_ratio"] = np.log2((g.loy_centre + PSEUDO) / (g.ctrl_centre + PSEUDO))
g["stratum"] = pd.qcut(np.log2(g.pooled_centre + PSEUDO), NDEC, labels=False)

targets = {"down": g.rna_class == "down", "down_no_lost_peak": (g.rna_class == "down") & (g.lost_peaks == 0),
           "up": g.rna_class == "up", "promoter_lost": g.lost_peaks > 0}
pool = g[(g.rna_class == "not_DE") & (g.lost_peaks == 0) & (g.gained_peaks == 0)]


def draw(tg, rng):
    """K controls per target within each stratum, without replacement; targets subsampled if a stratum is short."""
    keep, pick, acc = [], [], []
    for s in sorted(tg.stratum.unique()):
        t_s = tg[tg.stratum == s].sort_values("entrez_id")
        c_s = pool[pool.stratum == s].sort_values("entrez_id")
        n_t = min(len(t_s), len(c_s) // K)
        keep.append(t_s if n_t == len(t_s) else t_s.iloc[np.sort(rng.choice(len(t_s), n_t, replace=False))])
        pick.append(c_s.iloc[np.sort(rng.choice(len(c_s), n_t * K, replace=False))])
        acc.append(dict(stratum=int(s), targets=len(t_s), matched_targets=n_t, controls=n_t * K))
    return pd.concat(keep), pd.concat(pick), acc


def summarise(t_, c_):
    return dict(target_ratio=t_.loy_centre.mean() / t_.ctrl_centre.mean(),
                control_ratio=c_.loy_centre.mean() / c_.ctrl_centre.mean(),
                median_difference=t_.log2_ratio.median() - c_.log2_ratio.median(),
                MW_P=mannwhitneyu(t_.log2_ratio, c_.log2_ratio, alternative="two-sided").pvalue)


# Controls are drawn independently for each group (seed + group index), so a group's result does not depend on
# which other groups are analysed. The plotted draw is the first; NDRAW further draws measure draw-to-draw spread.
rows, acct, reps = [], [], []
for i, (name, mask) in enumerate(targets.items()):
    tg = g[mask]
    t_, c_, acc = draw(tg, np.random.default_rng(SEED + i))
    rows += [t_.assign(group=name), c_.assign(group=f"{name}_ref")]
    acct += [dict(group=name, **a) for a in acc]
    rng = np.random.default_rng(SEED + 1000 + i)
    for r in range(NDRAW):
        reps.append(dict(group=name, draw=r, **summarise(*draw(tg, rng)[:2])))
sel = pd.concat(rows, ignore_index=True)
acct = pd.DataFrame(acct)
reps = pd.DataFrame(reps)
robust = reps.groupby("group").agg(
    draws=("draw", "size"), control_ratio_median=("control_ratio", "median"),
    control_ratio_q025=("control_ratio", lambda x: x.quantile(.025)), control_ratio_q975=("control_ratio", lambda x: x.quantile(.975)),
    median_difference_median=("median_difference", "median"),
    median_difference_q025=("median_difference", lambda x: x.quantile(.025)),
    median_difference_q975=("median_difference", lambda x: x.quantile(.975)),
    MW_P_median=("MW_P", "median"), MW_P_max=("MW_P", "max"),
    pct_draws_P_below_005=("MW_P", lambda x: 100 * (x < .05).mean())).reset_index()

def smd(a, b):
    return (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2)

bal, tests, passage = [], [], []
for name in targets:
    t_, c_ = sel[sel.group == name], sel[sel.group == f"{name}_ref"]
    bal.append(dict(group=name, targets=len(t_), controls=len(c_),
                    SMD_log2_pooled_centre=smd(np.log2(c_.pooled_centre + PSEUDO), np.log2(t_.pooled_centre + PSEUDO)),
                    target_median_rna_logFC=t_.rna_logFC.median(), control_median_rna_logFC=c_.rna_logFC.median()))
    u = mannwhitneyu(t_.log2_ratio, c_.log2_ratio, alternative="two-sided")
    tests.append(dict(group=name, targets=len(t_), controls=len(c_),
                      target_ratio=t_.loy_centre.mean() / t_.ctrl_centre.mean(),
                      control_ratio=c_.loy_centre.mean() / c_.ctrl_centre.mean(),
                      target_median_log2=t_.log2_ratio.median(), control_median_log2=c_.log2_ratio.median(),
                      median_difference=t_.log2_ratio.median() - c_.log2_ratio.median(),
                      target_pct_below_control_median=100 * (t_.log2_ratio < c_.log2_ratio.median()).mean(),
                      MW_P=u.pvalue))
    for grp, d in [(name, t_), (f"{name}_ref", c_)]:
        for p in ["11", "12", "21"]:
            passage.append(dict(group=grp, passage=f"p{p}", ratio=d[f"LOY_{p}"].mean() / d[f"CTRL_{p}"].mean()))
        noCN = d[~d.chrom.isin(["chr13", "chr18"])]
        passage.append(dict(group=grp, passage="all, excluding chr13/chr18",
                            ratio=noCN.loy_centre.mean() / noCN.ctrl_centre.mean()))
tests = pd.DataFrame(tests)
tests["BH_FDR"] = multipletests(tests.MW_P, method="fdr_bh")[1]
bal, passage = pd.DataFrame(bal), pd.DataFrame(passage)
allg = dict(ratio_all_universe=float(g.loy_centre.mean() / g.ctrl_centre.mean()),
            median_log2_all_universe=float(g.log2_ratio.median()))

cols = ["group", "entrez_id", "gene_name", "chrom", "tss", "strand", "row", "rna_class", "rna_logFC", "rna_FDR",
        "lost_peaks", "gained_peaks", "atac_median_logFC", "stratum", "pooled_centre", *SAMPLES,
        "ctrl_centre", "loy_centre", "log2_ratio"]
sel[cols].to_csv(DATA / "tss_groups_genes.csv", index=False)
acct.to_csv(DATA / "tss_groups_matching_accounting.csv", index=False)
bal.to_csv(DATA / "tss_groups_matching_balance.csv", index=False)
tests.to_csv(DATA / "tss_groups_tests.csv", index=False)
passage.to_csv(DATA / "tss_groups_ratio_by_passage.csv", index=False)
robust.to_csv(DATA / "tss_groups_control_draws.csv", index=False)
json.dump(dict(universe=13300, genes_without_tss=n_no_tss, K=K, seed=SEED, robustness_draws=NDRAW, strata=f"deciles of log2 pooled central coverage",
               centre_bp=200, pseudocount=PSEUDO, **allg), open(DATA / "tss_groups_summary.json", "w"), indent=1)
print(f"genes without TSS row: {n_no_tss}")
print(acct.groupby("group")[["targets", "matched_targets", "controls"]].sum())
print(bal.round(3).to_string())
print(tests.round(4).to_string())
print(passage.pivot(index="group", columns="passage", values="ratio").round(3).to_string())
print(robust.round(4).to_string())
print(allg)
