"""Freeze HT1080 ATAC-seq figure inputs from the SEPT_08_2026 reprocess.

Reads completed analysis outputs only; the differential models are not refitted.
The one new calculation is the genomic-feature Fisher test against the tested
non-significant background (panel d of the main figure).
Run from any working directory with the Windows Python environment.
"""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[3]
ATAC = ROOT.parent / "ATACseq/SEPT_08_2026_REPROCESS"
OUT = ROOT / "figures/chapter3/source_data/ht1080_atac_accessibility"
OUT.mkdir(parents=True, exist_ok=True)
SAMPLES = ["CTRL_11", "CTRL_12", "CTRL_21", "LOY_11", "LOY_12", "LOY_21"]

# Representative promoter loci selected by the author on 22 September 2026:
# top-FDR promoter peaks at named genes (two gained, three lost).
LOCI = [("chr18:52339854-52340355", "DCC"), ("chr10:72297247-72297748", "DDIT4-AS1"),
        ("chr19:7294159-7294660", "INSR"), ("chr3:2098564-2099065", "CNTN4"),
        ("chr16:68448385-68448886", "SMPD3")]
LOCUS_HALF_WIDTH = 5000

inputs = {}
def read(rel, **kw):
    path = ATAC / rel
    inputs[rel] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path, **kw)

# QC and sample structure
qc = read("qc/qc_metrics_master.csv")
assert list(qc["sample"]) == SAMPLES
qc.to_csv(OUT / "library_qc.csv", index=False)
pc = read("qc/pca_scores.csv")
pcs = [c for c in pc.columns if c.startswith("PC")]
variance = pc[pcs].var() / pc[pcs].var().sum()
pc[["sample", "condition", "passage"] + pcs].to_csv(OUT / "pca_scores.csv", index=False)
variance.rename("variance_fraction").rename_axis("component").reset_index().to_csv(OUT / "pca_variance.csv", index=False)
read("qc/spearman_correlations.csv", index_col=0).loc[SAMPLES, SAMPLES].to_csv(OUT / "spearman_correlations.csv")
read("qc/matched_depth_peak_counts.tsv", sep="\t").to_csv(OUT / "matched_depth_peak_counts.csv", index=False)
frag, tss = [], []
for s in SAMPLES:
    f = read(f"qc/{s}_fragment_lengths.tsv", sep="\t"); f.insert(0, "sample", s); frag.append(f)
    t = read(f"qc/{s}_tss_profile.tsv", sep="\t"); t.insert(0, "sample", s); tss.append(t)
pd.concat(frag).to_csv(OUT / "fragment_lengths.csv", index=False)
pd.concat(tss).to_csv(OUT / "tss_profiles.csv", index=False)

# Differential accessibility (positive logFC = gain in Y-depleted)
ann = read("annotation/annotated_differential_peaks.csv")
csaw = read("differential/primary_TMM_csaw_support.csv")[["peak_id", "csaw_significant_same_direction"]]
bg = read("differential/background_TMM_all.csv")[["peak_id", "logFC", "FDR", "direction"]]
ps = read("differential/passage_TMM_all.csv")[["peak_id", "logFC", "FDR", "direction"]]
da = (ann[["peak_id", "chrom", "start", "end", "logFC", "logCPM", "PValue", "FDR",
           "TREAT_FDR_lfc0.5", "direction", "annotation", "nearest_symbol", "nearest_tss_distance_signed"]]
      .merge(csaw, on="peak_id", validate="1:1")
      .merge(bg.add_prefix("background_").rename(columns={"background_peak_id": "peak_id"}), on="peak_id", validate="1:1")
      .merge(ps.add_prefix("passage_").rename(columns={"passage_peak_id": "peak_id"}), on="peak_id", validate="1:1"))
assert len(da) == 84510
counts = da.direction.value_counts()
assert counts["gained"] == 95 and counts["lost"] == 663
da.to_csv(OUT / "differential_peaks.csv", index=False)
read("differential/analysis_summary.csv").to_csv(OUT / "analysis_summary.csv", index=False)
read("qc/peak_qc.tsv", sep="\t").to_csv(OUT / "peak_qc.csv", index=False)
csaw_regions = read("differential/csaw_regions.csv")
(csaw_regions.groupby("direction").size()
 .rename("regions").reset_index().to_csv(OUT / "csaw_region_counts.csv", index=False))

# Genomic-feature composition and Fisher tests: each direction vs tested ns peaks
rows = []
ns = da[da.direction == "ns"]
for direction in ["gained", "lost"]:
    fg = da[da.direction == direction]
    for feature in ["promoter", "gene_body", "intergenic"]:
        a, b = int((fg.annotation == feature).sum()), int((fg.annotation != feature).sum())
        c, d = int((ns.annotation == feature).sum()), int((ns.annotation != feature).sum())
        odds, p = fisher_exact([[a, b], [c, d]])
        rows.append(dict(direction=direction, feature=feature, foreground_in=a, foreground_total=a + b,
                         background_in=c, background_total=c + d, foreground_pct=100 * a / (a + b),
                         background_pct=100 * c / (c + d), odds_ratio=odds, fisher_p=p))
feat = pd.DataFrame(rows)
order = np.argsort(feat.fisher_p.to_numpy())
ranked = feat.fisher_p.to_numpy()[order] * len(feat) / np.arange(1, len(feat) + 1)
bh = np.empty(len(feat)); bh[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
feat["fisher_BH"] = bh
feat.to_csv(OUT / "feature_enrichment_fisher.csv", index=False)

# Locus definitions (coverage extracted separately under WSL with pyBigWig)
genes = read("config/genes_knownGene.tsv", sep="\t", dtype={"gene_id": str})
links = read("annotation/promoter_peak_gene_links.csv", dtype={"gene_id": str})
loci = []
for peak_id, label in LOCI:
    row = da.set_index("peak_id").loc[peak_id]
    mid = (row.start + row.end) // 2
    loci.append(dict(peak_id=peak_id, label=label, chrom=row.chrom, peak_start=row.start, peak_end=row.end,
                     window_start=mid - LOCUS_HALF_WIDTH, window_end=mid + LOCUS_HALF_WIDTH,
                     logFC=row.logFC, FDR=row.FDR, direction=row.direction,
                     promoter_genes=";".join(sorted(set(links[links.peak_id == peak_id].symbol.dropna())))))
loci = pd.DataFrame(loci)
loci.to_csv(OUT / "locus_definitions.csv", index=False)
gene_rows = []
for r in loci.itertuples():
    g = genes[(genes.chrom == r.chrom) & (genes.end > r.window_start) & (genes.start < r.window_end)]
    gene_rows.append(g.assign(peak_id=r.peak_id))
gene_rows = pd.concat(gene_rows)
from_entrez = links.drop_duplicates("gene_id").set_index("gene_id").symbol
gene_rows["symbol"] = gene_rows.gene_id.map(from_entrez)
gene_rows.to_csv(OUT / "locus_genes_unlabelled.csv", index=False)

json.dump({"source_run": str(ATAC), "input_sha256": inputs,
           "notes": "Fisher tests are new: each direction vs tested non-significant peaks, two-sided, BH across 6 tests."},
          open(OUT / "provenance.json", "w"), indent=1)
print(feat.round(4).to_string()); print(loci.to_string())
