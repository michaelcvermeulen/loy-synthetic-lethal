# `data/`

Raw sequencing files (`*.bam`, `*.bw`), flow cytometry `.fcs` files, and large serialized R objects (`*.rds`) are excluded via `.gitignore`. Derived summary tables and manifests are stored within [`../analysis/`](../analysis/).

```text
data/
├── depmap_prism/          # DepMap, Sanger Score, NextGen 3D, PRISM, GDSC1/2, CTD2
├── sgrna_t2t_y/           # T2T-CHM13v2.0 chrY reference & CenSat v2.1 annotations
├── ddpcr/                 # QX200 ddPCR droplet counts & concentration tables
├── lpwgs/                 # T2T-CHM13v2.0 (Y-PAR-masked) BAMs & 100-kb bin tracks
├── rnaseq/                # STAR v2 (Y-PAR-masked, GENCODE v49) count matrices & external sets
├── atacseq/               # HT1080 BAMs, MACS3 narrowPeak consensus counts & TOBIAS tracks
├── crispr_screens/        # HT1080 & A549 Brunello sgRNA count matrices + external screens
└── comp_assay/            # BD FACSDiscover A8 & CytoFLEX S .fcs files
```
