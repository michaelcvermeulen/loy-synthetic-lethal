"""Copy-number context of HT1080 differential accessibility (all tested regions and RNA-ATAC joint genes).

Inputs (no refit):
  differential/primary_TMM_all.csv                              84,510 tested regions, primary edgeR model
  integration/RNA_interaction_CN/WGS_hg38_100kb_bins.csv        lpWGS D0 bins, delta = log2(LOY seg / CTRL seg)
  SEPT_11_2026_THESIS/.../ht1080_rna_atac_integration/integration_genes.csv   13,300-gene universe (joint_010)
Each region gets the length-weighted segment delta of the usable bins it overlaps (same rule as
73_interaction_CN_report.interval_context; filtered bins are left missing, not neutral).
CN classes: loss (delta <= -T), neutral (|delta| < T), gain (delta >= T), uncovered; T = 0.20 primary, 0.10/0.30 sensitivity.
One WGS sample per condition: deltas are descriptive relative dosage flags, not tested CN calls.
Outputs: integration/thesis_RNA_SEPT11/copy_number_context/
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "integration/thesis_RNA_SEPT11/copy_number_context"
OUT.mkdir(parents=True, exist_ok=True)
GENES = ROOT.parents[1] / "SEPT_11_2026_THESIS/figures/chapter3/source_data/ht1080_rna_atac_integration/integration_genes.csv"


def interval_context(chrom, start, end, bychrom):
    """Length-weighted segment delta over usable bins (copied from 73_interaction_CN_report.py)."""
    res = dict(covered_fraction=0., delta_mean=np.nan, bins=0)
    if chrom not in bychrom or end <= start:
        return res
    b = bychrom[chrom]
    i = np.searchsorted(b["end0"], start, side="right"); j = np.searchsorted(b["start0"], end, side="left")
    if j <= i:
        return res
    w = np.minimum(b["end0"][i:j], end) - np.maximum(b["start0"][i:j], start); v = b["delta_segment_log2"][i:j]
    good = (w > 0) & b["use"][i:j] & np.isfinite(v); w = w[good]; v = v[good]
    if not len(v):
        return res
    res.update(covered_fraction=float(w.sum() / (end - start)), delta_mean=float(np.average(v, weights=w)), bins=len(v))
    return res


def cn_class(d, t):
    return np.select([d.isna(), d <= -t, d >= t], ["uncovered", "loss", "gain"], "neutral")


bins = pd.read_csv(ROOT / "integration/RNA_interaction_CN/WGS_hg38_100kb_bins.csv",
                   usecols=["chrom", "start0", "end0", "use", "delta_segment_log2"])
bychrom = {c: {k: f[k].to_numpy() for k in ["start0", "end0", "delta_segment_log2", "use"]} for c, f in bins.groupby("chrom", sort=False)}
peaks = pd.read_csv(ROOT / "differential/primary_TMM_all.csv")
ctx = pd.DataFrame([interval_context(r.chrom, r.start, r.end, bychrom) for r in peaks.itertuples()])
p = pd.concat([peaks[["peak_id", "chrom", "start", "end", "logFC", "logCPM", "FDR", "direction"]], ctx.add_prefix("CN_")], axis=1)
p["call"] = np.where(p.FDR < .05, np.where(p.logFC < 0, "lost", "gained"), "ns")
p["call_010"] = np.where(p.FDR < .10, np.where(p.logFC < 0, "lost", "gained"), "ns")
for t in (.1, .2, .3):
    p[f"CN_class_{t:.2f}"] = cn_class(p.CN_delta_mean, t)
p.to_csv(OUT / "region_CN_context.csv", index=False)

# 1. call rates by CN class, with Fisher tests of each class vs neutral
rows = []
for t in (.1, .2, .3):
    col = f"CN_class_{t:.2f}"
    for fdr, callcol in [(.05, "call"), (.10, "call_010")]:
        neu = p[p[col] == "neutral"]
        for cls in ["loss", "neutral", "gain", "uncovered"]:
            s = p[p[col] == cls]
            r = dict(threshold=t, ATAC_FDR=fdr, CN_class=cls, regions=len(s), lost=int((s[callcol] == "lost").sum()),
                     gained=int((s[callcol] == "gained").sum()), median_logFC=s.logFC.median(),
                     median_CN_delta=s.CN_delta_mean.median())
            for d in ("lost", "gained"):
                r[f"pct_{d}"] = 100 * r[d] / max(len(s), 1)
                if cls in ("loss", "gain") and len(s):
                    a, b = r[d], len(s) - r[d]; c = int((neu[callcol] == d).sum()); e = len(neu) - c
                    orr, pv = stats.fisher_exact([[a, b], [c, e]])
                    r[f"OR_{d}_vs_neutral"], r[f"P_{d}_vs_neutral"] = orr, pv
            rows.append(r)
rates = pd.DataFrame(rows)
prim = rates[(rates.threshold == .2)]
pcols = [c for c in rates.columns if c.startswith("P_")]
flat = prim[pcols].to_numpy().ravel(); ok = ~np.isnan(flat)
adj = np.full(flat.shape, np.nan); adj[ok] = stats.false_discovery_control(flat[ok])
adj = adj.reshape(prim[pcols].shape)
for k, c in enumerate(pcols):
    rates.loc[prim.index, c.replace("P_", "BH_")] = adj[:, k]
rates.to_csv(OUT / "call_rates_by_CN_class.csv", index=False)

# 2. continuous association: region logFC vs CN delta; chromosome medians
cov = p[p.CN_delta_mean.notna()]
rho_all = stats.spearmanr(cov.logFC, cov.CN_delta_mean)
chrom = (cov.groupby("chrom").agg(regions=("peak_id", "size"), median_logFC=("logFC", "median"),
                                  median_CN_delta=("CN_delta_mean", "median"),
                                  lost=("call", lambda s: (s == "lost").sum()), gained=("call", lambda s: (s == "gained").sum()))
         .reset_index())
chrom["pct_lost"] = 100 * chrom.lost / chrom.regions; chrom["pct_gained"] = 100 * chrom.gained / chrom.regions
chrom.to_csv(OUT / "chromosome_summary.csv", index=False)
rho_chr = stats.spearmanr(chrom.median_logFC, chrom.median_CN_delta)
rho_chr_lost = stats.spearmanr(chrom.pct_lost, chrom.median_CN_delta)

# 3. RNA-ATAC joint genes (ATAC FDR 0.10): CN at their promoter peaks
g = pd.read_csv(GENES)
pk = p.set_index("peak_id").CN_delta_mean
g["CN_promoter_delta"] = [np.nanmedian(pk.reindex(str(s).split(";")).to_numpy()) if pk.reindex(str(s).split(";")).notna().any() else np.nan
                          for s in g.peak_ids]
g["CN_class_0.20"] = cn_class(g.CN_promoter_delta, .2)
g["is_joint_010"] = g.joint_010.str.contains("concordant|opposite", case=False, na=False)
g[["entrez_id", "gene_name", "chrom", "rna_class", "atac_class_010", "joint_010", "is_joint_010", "CN_promoter_delta", "CN_class_0.20"]
  ].to_csv(OUT / "gene_CN_context.csv", index=False)
jt = pd.crosstab(g.joint_010, g["CN_class_0.20"]).reset_index()
jt.to_csv(OUT / "joint_genes_by_CN_class.csv", index=False)
j = g[g.is_joint_010]
chr19 = j[j.chrom == "chr19"]

# Sensitivity: concordant enrichment among DE-down genes restricted to CN-neutral promoters
down = g[g.rna_class.astype(str).str.contains("lower|down", case=False)]
sens = {}
for name, sub in [("all", down), ("CN_neutral_0.20", down[down["CN_class_0.20"] == "neutral"])]:
    sens[name] = dict(DE_lower=len(sub), with_lost_peak_010=int((sub.lost_010 > 0).sum()),
                      pct=100 * (sub.lost_010 > 0).mean())
nonde = g[g.rna_class.astype(str).str.contains("not", case=False)]
for name, sub in [("all", nonde), ("CN_neutral_0.20", nonde[nonde["CN_class_0.20"] == "neutral"])]:
    sens["nonDE_" + name] = dict(genes=len(sub), with_lost_peak_010=int((sub.lost_010 > 0).sum()), pct=100 * (sub.lost_010 > 0).mean())

sig = p[p.call != "ns"]
summary = dict(
    regions=len(p), regions_CN_covered=int(p.CN_delta_mean.notna().sum()),
    CN_class_counts_0_20=p["CN_class_0.20"].value_counts().to_dict(),
    significant_FDR05=len(sig), significant_by_CN_class_0_20=pd.crosstab(sig.call, sig["CN_class_0.20"]).to_dict(),
    pct_significant_in_neutral_0_20=100 * (sig["CN_class_0.20"] == "neutral").mean(),
    spearman_logFC_vs_CN_delta=dict(rho=rho_all.statistic, P=rho_all.pvalue, n=len(cov)),
    spearman_chromosome_median_logFC_vs_CN=dict(rho=rho_chr.statistic, P=rho_chr.pvalue, n=len(chrom)),
    spearman_chromosome_pct_lost_vs_CN=dict(rho=rho_chr_lost.statistic, P=rho_chr_lost.pvalue),
    joint_010_genes=len(j), joint_010_by_CN_class=j["CN_class_0.20"].value_counts().to_dict(),
    joint_010_chr19=len(chr19), joint_010_chr19_CN_delta_range=[float(chr19.CN_promoter_delta.min()), float(chr19.CN_promoter_delta.max())],
    lost_peak_rate_sensitivity=sens,
)
json.dump(summary, open(OUT / "summary.json", "w"), indent=2, default=float)
print(json.dumps(summary, indent=2, default=float))
print(rates[rates.threshold == .2].round(4).to_string(index=False))
print(chrom.round(3).to_string(index=False))

# 4. Chromosome-level shifts within each passage (ATAC LOY - CTRL logCPM, centred on the genome median),
#    RNA median logFC per chromosome (13,300-gene universe, centred) and WGS D0 median segment delta.
lc = pd.read_csv(ROOT / "counts/logCPM_primary.tsv", sep="\t")
lc["chrom"] = lc.peak_id.str.split(":").str[0]
shift = {}
for ps in ("11", "12", "21"):
    d = lc[f"LOY_{ps}"] - lc[f"CTRL_{ps}"]
    shift[f"ATAC_p{ps}"] = d.groupby(lc.chrom).median() - d.median()
shift = pd.DataFrame(shift)
shift["RNA_median_logFC"] = g.groupby("chrom").rna_logFC.median() - g.rna_logFC.median()
wgs = pd.read_csv(ROOT / "integration/RNA_interaction_CN/WGS_chromosome_summary.csv").set_index("chrom")
shift["WGS_D0_median_segment_delta"] = wgs.median_segment_delta
shift["gained_FDR05"] = p[p.call == "gained"].chrom.value_counts()
shift["lost_FDR05"] = p[p.call == "lost"].chrom.value_counts()
shift = shift.fillna({"gained_FDR05": 0, "lost_FDR05": 0}).rename_axis("chrom").reset_index()
shift.to_csv(OUT / "chromosome_shift_by_passage.csv", index=False)
auto = shift[shift.chrom.str.fullmatch(r"chr\d+")]
rho_atac_wgs = stats.spearmanr(auto[["ATAC_p11", "ATAC_p12", "ATAC_p21"]].mean(axis=1), auto.WGS_D0_median_segment_delta)
rho_rna_wgs = stats.spearmanr(auto.RNA_median_logFC, auto.WGS_D0_median_segment_delta)
rho_atac_rna = stats.spearmanr(auto[["ATAC_p11", "ATAC_p12", "ATAC_p21"]].mean(axis=1), auto.RNA_median_logFC)

# 5. Motif sensitivity: matched CMH tests (as in the differential-heatmap supplement) with targets on
#    chr18 (gained) or chr13/chr18 (lost) removed; controls unchanged. statsmodels StratifiedTable: MH OR, 95% CI, CMH test (continuity-corrected).
REG = ROOT.parents[1] / "SEPT_11_2026_THESIS/figures/chapter3/source_data/ht1080_atac_accessibility/region_heatmap_regions.csv"
reg = pd.read_csv(REG)


def cmh(sub, grp, m):
    tabs = []
    for _, d in sub.groupby("stratum"):
        t = d[d.group == grp]; c = d[d.group == grp + "_ref"]
        if len(t) and len(c):
            tabs.append(np.array([[t[m].sum(), len(t) - t[m].sum()], [c[m].sum(), len(c) - c[m].sum()]], dtype=float))
    st = StratifiedTable(tabs); lo, hi = st.oddsratio_pooled_confint()
    T = sub[sub.group == grp]; C = sub[sub.group == grp + "_ref"]
    return dict(targets=len(T), pct_targets=100 * T[m].mean(), pct_controls=100 * C[m].mean(), MH_OR=st.oddsratio_pooled,
                CI_low=lo, CI_high=hi, P=st.test_null_odds(correction=True).pvalue)


mrows = []
for grp, motifs, drop in [("lost", ["TFAP2", "EBF", "ZIC", "TEAD", "RUNX", "GATA"], ["chr13", "chr18"]), ("gained", ["SOX6"], ["chr18"])]:
    sub = reg[reg.group.str.startswith(grp)]
    for m in motifs:
        for lab, s in [("all targets", sub), ("excluding " + "/".join(drop), sub[~((sub.group == grp) & sub.chrom.isin(drop))])]:
            mrows.append(dict(direction=grp, motif=m, target_set=lab, **cmh(s, grp, m)))
mot = pd.DataFrame(mrows)
mot["BH_FDR"] = stats.false_discovery_control(mot.P)
mot.to_csv(OUT / "motif_CMH_chromosome_exclusion.csv", index=False)

summary.update(
    chromosome_shift_spearman=dict(ATAC_vs_WGS=[rho_atac_wgs.statistic, rho_atac_wgs.pvalue],
                                   RNA_vs_WGS=[rho_rna_wgs.statistic, rho_rna_wgs.pvalue],
                                   ATAC_vs_RNA=[rho_atac_rna.statistic, rho_atac_rna.pvalue], n_autosomes=len(auto)),
    chr18=dict(gained_FDR05=int((p[p.call == "gained"].chrom == "chr18").sum()), gained_total=int((p.call == "gained").sum()),
               ATAC_shift_by_passage=shift.set_index("chrom").loc["chr18", ["ATAC_p11", "ATAC_p12", "ATAC_p21"]].tolist(),
               RNA_shift=float(shift.set_index("chrom").loc["chr18", "RNA_median_logFC"]),
               WGS_D0=float(shift.set_index("chrom").loc["chr18", "WGS_D0_median_segment_delta"])),
)
json.dump(summary, open(OUT / "summary.json", "w"), indent=2, default=float)
print(shift.round(3).to_string(index=False))
print(summary["chromosome_shift_spearman"])
print(mot.round(4).to_string(index=False))
