"""HOMER comparison arm: all FDR < 0.05 peaks, no fold-change filter (author request, 22 Sep 2026).

Identical to the primary TF-extension HOMER arm (21_prepare_tf_inputs.py, 25_run_homer.py, 26_summarize_homer.py)
except that targets are primary_TMM peaks with FDR < 0.05 regardless of |log2FC| (663 lost, 95 gained, minus
N-containing peaks), instead of FDR < 0.05 and |log2FC| > 1 (609 / 47). Background: tested consensus FDR >= 0.20,
no N bases, same 5-point GC bins x logCPM deciles, up to 10 controls per target without replacement (seed
20260922, so control draws differ from the primary arm). Same JASPAR2024 HOMER motif file and options.
Existing primary results are not modified. Output: motifs/TF_extension/HOMER/primary_TMM_allFDR05/.

Run in WSL with python3 after restoring the TF tools:
  python3 code/84_homer_all_significant.py setup     # unpack verified HOMER archive, restore verified genome cache
  python3 code/84_homer_all_significant.py run
  python3 code/84_homer_all_significant.py summarize (Windows or WSL; pandas required)
"""
from pathlib import Path
import sys, os, json, re, subprocess, hashlib, shutil, time

ROOT = Path(__file__).resolve().parents[1]
TFX = ROOT / "motifs/TF_extension"
# Optional second argument "promoter": restrict targets AND background to promoter peaks (knownGene TSS +/- 2 kb,
# annotation/annotated_differential_peaks.csv), so enrichment is promoter-vs-promoter (author request, 22 Sep 2026).
SCOPE = sys.argv[2] if len(sys.argv) > 2 else "all"
assert SCOPE in ("all", "promoter", "distal")  # distal = gene body + intergenic (not TSS +/- 2 kb)
# Optional third argument: target peak FDR threshold (default 0.05; 0.10 added at author request, 22 Sep 2026).
TARGET_FDR = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
ARM = "primary_TMM_allFDR05" + ("" if SCOPE == "all" else f"_{SCOPE}") + ("" if TARGET_FDR == 0.05 else f"_fdr{round(TARGET_FDR * 100):03d}")
OUT = TFX / "HOMER" / ARM
TOOLS = Path("/tmp/atacseq-sept08-2026-tf-tools")
BIN = TOOLS / "homer_package/bin"
stage = sys.argv[1] if len(sys.argv) > 1 else "run"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


if stage == "setup":
    ref = json.loads((ROOT / "config/TF_extension/reference_verified.json").read_text())
    cache = Path(ref["linux_cache"])
    if not (BIN / "findMotifsGenome.pl").exists():
        subprocess.run([sys.executable, str(ROOT / "code/20b_unpack_homer.py")], check=True)
    if not cache.exists() or sha256(cache) != ref["sha256"]:
        cache.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ref["source"], cache)
    assert sha256(cache) == ref["sha256"], "genome cache hash mismatch"
    print("setup ok:", BIN / "findMotifsGenome.pl", cache)

elif stage == "run":
    import numpy as np, pandas as pd
    rng = np.random.default_rng(20260922)
    cov = pd.read_csv(ROOT / "motifs/peak_covariates.csv")
    df = pd.read_csv(ROOT / "differential/primary_TMM_all.csv").merge(cov[["peak_id", "GC", "N_fraction"]], on="peak_id",
                                                                      validate="one_to_one")
    df = df[df.N_fraction == 0].copy()
    if SCOPE != "all":
        ann = pd.read_csv(ROOT / "annotation/annotated_differential_peaks.csv", usecols=["peak_id", "annotation"])
        df = df.merge(ann, on="peak_id", validate="one_to_one")
        df = df[(df.annotation == "promoter") if SCOPE == "promoter" else (df.annotation != "promoter")].copy()
    df["gc_bin"] = np.minimum((df.GC / .05).astype(int), 19)
    df["accessibility_bin"] = pd.qcut(df.logCPM, 10, labels=False, duplicates="drop")
    OUT.mkdir(parents=True, exist_ok=True)
    selected = df.FDR < TARGET_FDR
    accounts = []
    for direction in ["gained", "lost"]:
        q = df[selected & ((df.logFC > 0) if direction == "gained" else (df.logFC < 0))]
        b = df[df.FDR >= .2]
        ratios = []
        for key, g in q.groupby(["gc_bin", "accessibility_bin"]):
            available = ((b.gc_bin == key[0]) & (b.accessibility_bin == key[1])).sum()
            if available:
                ratios.append(max(1, int(available) // len(g)))
        control_ratio = min([10] + ratios)
        qparts, bparts = [], []
        for key, g in q.groupby(["gc_bin", "accessibility_bin"]):
            pool = b[(b.gc_bin == key[0]) & (b.accessibility_bin == key[1])]
            if pool.empty:
                continue
            keep = min(len(g), len(pool))
            qparts.append(g if keep == len(g) else g.iloc[rng.choice(len(g), keep, replace=False)])
            bparts.append(pool.iloc[rng.choice(len(pool), control_ratio * keep, replace=False)])
        qt, bt = pd.concat(qparts), pd.concat(bparts)
        assert set(qt.peak_id).isdisjoint(bt.peak_id)
        for data, label in [(qt, "target"), (bt, "background")]:
            data.to_csv(OUT / f"{direction}_{label}.csv", index=False)
            data.assign(score=0, strand="+")[["chrom", "start", "end", "peak_id", "score", "strand"]].to_csv(
                OUT / f"{direction}_{label}.bed", sep="\t", header=False, index=False)
        accounts.append(dict(arm=ARM, direction=direction, eligible_target=len(q), matched_target=len(qt),
                             unmatched_target=len(q) - len(qt), background=len(bt), controls_per_target=control_ratio))
    pd.DataFrame(accounts).to_csv(OUT / "matching_accounting.csv", index=False)
    print(pd.DataFrame(accounts))
    ref = json.loads((ROOT / "config/TF_extension/reference_verified.json").read_text())
    env = os.environ.copy(); env["PATH"] = str(BIN) + os.pathsep + env["PATH"]
    expected_seq = dict(l.rstrip("\n").split("\t") for l in (ROOT / "config/TF_extension/tested_peak_sequences.tsv").open()
                        if not l.startswith("peak_id\t"))
    for direction in ["gained", "lost"]:
        folder = OUT / direction; folder.mkdir(exist_ok=True)
        args = [str(BIN / "findMotifsGenome.pl"), str(OUT / f"{direction}_target.bed"), ref["linux_cache"], str(folder),
                "-size", "given", "-bg", str(OUT / f"{direction}_background.bed"),
                "-mknown", str(ROOT / "config/TF_extension/JASPAR2024_CORE_vertebrates.homer"),
                "-nomotif", "-useOldBg", "-nofacts", "-keepFiles", "-p", "6",
                "-preparsedDir", str(TOOLS / "homer_preparsed")]
        (folder / "command.json").write_text(json.dumps(args, indent=2) + "\n")
        start = time.time()
        with (ROOT / f"logs/84_HOMER_{ARM}_{direction}.log").open("w") as log:
            subprocess.run(args, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        header = (folder / "knownResults.txt").read_text().splitlines()[0]
        n_t = sum(1 for _ in (OUT / f"{direction}_target.bed").open())
        assert int(re.search(r"Target Sequences with Motif\(of (\d+)\)", header).group(1)) == n_t
        checked = 0
        for line in (folder / "seq.tsv").read_text().splitlines():
            sid, seq = line.split("\t")[:2]; pid = sid[2:] if sid.startswith("BG") else sid
            assert seq.upper() == expected_seq[pid], f"sequence mismatch {pid}"
            checked += 1
        assert checked == n_t + sum(1 for _ in (OUT / f"{direction}_background.bed").open())
        (folder / "COMPLETE.json").write_text(json.dumps(dict(elapsed_seconds=time.time() - start, version="5.1",
                                                             command=args), indent=2) + "\n")
        print("HOMER", direction, "done", flush=True)

elif stage == "summarize":
    import numpy as np, pandas as pd
    mapping = pd.read_csv(TFX / "chromVAR/motif_archetype_mapping.csv")
    def bh(p):
        p = np.asarray(p); o = np.argsort(p); v = p[o] * len(p) / np.arange(1, len(p) + 1)
        a = np.minimum.accumulate(v[::-1])[::-1].clip(0, 1); out = np.empty(len(p)); out[o] = a; return out
    pieces = []
    for direction in ["gained", "lost"]:
        folder = OUT / direction; assert (folder / "COMPLETE.json").exists()
        raw = pd.read_csv(folder / "knownResults.txt", sep="\t")
        n_t = int(re.search(r"\(of (\d+)\)", raw.columns[5]).group(1)); n_b = int(re.search(r"\(of (\d+)\)", raw.columns[7]).group(1))
        d = pd.DataFrame(dict(motif_id=raw.iloc[:, 0].str.extract(r"^(MA\d+\.\d+)", expand=False),
                              log_PValue=pd.to_numeric(raw.iloc[:, 3]), target_motif_count=pd.to_numeric(raw.iloc[:, 5]),
                              background_motif_count=pd.to_numeric(raw.iloc[:, 7])))
        d["PValue"] = np.exp(d.log_PValue)
        assert len(d) == 879 and d.motif_id.is_unique
        d["target_total"] = n_t; d["background_total"] = n_b
        d["target_fraction"] = d.target_motif_count / n_t; d["background_fraction"] = d.background_motif_count / n_b
        d["fold_enrichment"] = d.target_fraction.div(d.background_fraction.replace(0, np.nan))
        d["arm"] = ARM; d["direction"] = direction; pieces.append(d)
    res = pd.concat(pieces, ignore_index=True); res["BH_FDR_both_directions"] = bh(res.PValue)
    res = res.merge(mapping, on="motif_id", validate="many_to_one").sort_values(["BH_FDR_both_directions", "PValue"])
    res.to_csv(OUT / "known_motif_all_results.csv", index=False)
    if SCOPE == "all":
        prim = pd.read_csv(TFX / "HOMER/known_motif_all_results.csv"); prim = prim[prim.arm == "primary_TMM"]
    else:  # promoter / distal scopes are compared with the all-peak FDR < 0.05 arm
        prim = pd.read_csv(TFX / "HOMER/primary_TMM_allFDR05/known_motif_all_results.csv")
    cmp = res.merge(prim[["motif_id", "direction", "fold_enrichment", "BH_FDR_both_directions", "target_total"]],
                    on=["motif_id", "direction"], suffixes=("_all", "_FC1"))
    cmp.to_csv(OUT / ("comparison_with_FC1_arm.csv" if SCOPE == "all" else "comparison_with_all_peak_arm.csv"), index=False)
    summ = cmp.groupby("direction").apply(lambda x: pd.Series(dict(
        targets_all=int(x.target_total_all.iloc[0]), targets_FC1=int(x.target_total_FC1.iloc[0]),
        sig_all=int((x.BH_FDR_both_directions_all < .05).sum()), sig_FC1=int((x.BH_FDR_both_directions_FC1 < .05).sum()),
        sig_both=int(((x.BH_FDR_both_directions_all < .05) & (x.BH_FDR_both_directions_FC1 < .05)).sum()),
        spearman_log2FE=np.log2(x.fold_enrichment_all.clip(lower=.25)).corr(np.log2(x.fold_enrichment_FC1.clip(lower=.25)),
                                                                           method="spearman"))), include_groups=False)
    summ.to_csv(OUT / "comparison_summary.csv")
    print(summ)
    for direction in ["gained", "lost"]:
        s = cmp[(cmp.direction == direction) & ((cmp.BH_FDR_both_directions_all < .05) | (cmp.BH_FDR_both_directions_FC1 < .05))]
        s = s.sort_values("BH_FDR_both_directions_all").drop_duplicates("archetype").head(14)
        print(f"\n{direction}:")
        print(s[["motif_name", "archetype", "target_fraction", "background_fraction", "fold_enrichment_all",
                 "BH_FDR_both_directions_all", "fold_enrichment_FC1", "BH_FDR_both_directions_FC1"]].round(4).to_string(index=False))
