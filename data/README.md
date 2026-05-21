# data/

**This file is the contract.** To reproduce the analyses in this repository,
place input files matching the layout below into `data/`. Nothing in this
directory is tracked by git except this README — the data itself stays local
(or, post-publication, is pulled from the public accessions listed below).

```
data/
├── crispr/
├── rnaseq/
├── atacseq/
├── chipseq/
├── scrnaseq/
├── depmap/
└── encode/
```

---

## crispr/ — CRISPR screen counts

- **Experiment**: Brunello genome-wide library, HT1080 isogenic LOY vs sgCTRL.
- **Files expected**: raw read counts per sgRNA (one file per sample / time
  point), plus a sample metadata sheet.
- **Filename convention**: TODO.
- **Source**: in-house (Babak Lab). Accession: `TODO: deposition target +
  accession at publication`.

## rnaseq/ — bulk RNA-seq

- **Experiments**: HT1080 parental vs LOY, LNCaP, HepG2, plus add-back
  conditions (KDM5D / UTY rescue).
- **Files expected**: gene-level count matrices (e.g., `salmon`/`STAR` output)
  per experiment, plus per-experiment sample metadata.
- **Filename convention**: TODO.
- **Source**: in-house. Accession: `TODO: GEO accession (to be deposited at
  manuscript submission)`.

## atacseq/ — bulk ATAC-seq

- **Experiment**: HT1080 parental vs LOY (Active Motif service).
- **Files expected**: processed BAMs, peak calls (narrowPeak), bigWig signal
  tracks.
- **Filename convention**: TODO.
- **Source**: Active Motif processing of in-house libraries. Accession:
  `TODO: GEO accession`.

## chipseq/ — ChIP-seq

- **Experiments**:
  - In-house, Postovit Lab collaboration: TODO mark targets and conditions.
  - Public: **GSE261834**, **GSE135580**.
- **Files expected**: BAMs / bigWigs / peak calls for each mark and condition.
- **Filename convention**: TODO.
- **Source**: see above; public sets fetched from GEO.

## scrnaseq/ — single-cell RNA-seq

- **Experiments**: public LOY atlas datasets re-analyzed here.
- **Files expected**: filtered count matrices (`.h5ad` or `10X` MTX bundles)
  plus per-dataset cell metadata.
- **Filename convention**: TODO.
- **Source**: GEO. Accessions: `TODO: list GEO accessions`.

## depmap/ — DepMap pan-cancer dependency data

- **Files expected**: CRISPR gene-effect matrix, cell-line metadata,
  expression matrix.
- **Filename convention**: as published by DepMap.
- **Source**: DepMap portal. Release: `TODO: pin DepMap release (e.g., 24Q4)`.

## encode/ — ENCODE tissue ATAC-seq

- **Files expected**: per-tissue ATAC signal tracks (bigWig) and peak calls
  (narrowPeak) for the panel used in analysis 08.
- **Filename convention**: ENCODE portal IDs preserved (`ENCFF*`).
- **Source**: ENCODE portal. IDs: `TODO: list ENCFF/ENCSR portal IDs`.
