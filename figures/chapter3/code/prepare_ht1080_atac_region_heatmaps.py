"""HOMER statistics for the representative motifs in the differential-region figure.

Region sets, GC + accessibility matching, motif presence and CMH tests are made by
prepare_ht1080_atac_region_sets.R. The earlier logCPM-only matching version of this
script is archived in archive/atac_region_matching_logCPM_only_2026-09-22/.
HOMER targets are the FDR < 0.05 and |log2FC| > 1 subsets (609 lost, 47 gained).
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "figures/chapter3/source_data/ht1080_atac_accessibility"
reps = pd.DataFrame({"family": ["TFAP2", "EBF", "ZIC", "TEAD", "RUNX", "GATA", "SOX6"],
                     "motif_id": ["MA0524.3", "MA0154.5", "MA0696.1", "MA0090.4", "MA0002.3", "MA1104.3", "MA0515.1"]})
homer = pd.read_csv(ROOT.parent / "ATACseq/SEPT_08_2026_REPROCESS/motifs/TF_extension/HOMER/known_motif_all_results.csv")
homer = homer[homer.arm == "primary_TMM"].merge(reps, on="motif_id")
assert len(homer) == 14
homer[["family", "motif_id", "motif_name", "direction", "target_total", "background_total", "target_fraction",
       "background_fraction", "fold_enrichment", "PValue", "BH_FDR_both_directions", "archetype"]].to_csv(
    DATA / "region_heatmap_homer_stats.csv", index=False)
