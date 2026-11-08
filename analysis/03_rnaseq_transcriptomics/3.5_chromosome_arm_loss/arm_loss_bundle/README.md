# Chromosome-arm loss: is there a common transcriptional response?

Bundle date: 2026-09-29. Two linked analyses:

1. **GSE336056** (single-cell, RPTEC, CRISPR): the chromosome Yq arm and the 3p arm are each cut away in primary kidney tubule cells from one donor, and cells that lost the arm are compared with cells from the same guide that did not.
2. **Arm-loss meta-analysis** (bulk RNA-seq, 15 engineered losses in 10 cell backgrounds): asks whether the p53-up / proliferation-down / interferon response seen in (1) is a general consequence of losing any chromosome arm.

## Main findings

**GSE336056, single cell (days 5 and 10 after editing)**
- Yq loss: 194 Yq-null cells vs 6,407 Yq-intact cells within the sgYq arm. The break falls between PRKY (retained) and USP9Y (lost), at 7.4–12.5 Mb, which removes essentially the whole Yq arm. Yp is kept.
- 3p loss: 649 3p-loss cells vs 8,934 intact cells within the sg3p arm. One copy of the whole arm is lost (chr3p genes at log2 ratio −0.91).
- The two losses share most of their pathway response: 16 of the 20 Hallmark pathways significant in Yq are also significant in 3p, all in the same direction, and none go the opposite way (Hallmark NES Spearman ρ = 0.71). Down: E2F, G2M, MYC and mitotic spindle. Up: p53, apoptosis, TNFα–NF-κB, inflammatory response and EMT.
- The main difference is type I interferon: up in 3p loss, down in Yq loss.

**Meta-analysis, bulk (clones expanded after editing)**

| Question | Answer |
|---|---|
| Is there a program shared by all arm losses? | **No.** No pathway passes a sign test across the 15 contrasts (best is 13/15). The NES profiles of different arm losses don't resemble each other either: median pairwise ρ = −0.00 over 105 pairs. |
| Is the arm-loss response just p53 activation? | **No.** Compared with the p53-induction profile at fixed karyotype, median ρ = 0.02 (maximum 0.35, n = 30 pairs). |
| Does the response grow with the number of genes lost? | **No.** At the segment level (n = 7 comparable segments), p53 ρ = +0.50, P = 0.27. The contrast-level ρ of +0.73 came from counting the same 1q segment five times. |
| What replicates? | **Ribosomal proteins go down**: 9 of 15 contrasts significantly down, **none** up. It is the only axis that leans one way across the cohort. |
| Early vs late clones | **The response changes with time.** In the RPTEC 3p clones, between day 18 and day 90 interferon-alpha goes from +2.08 to -1.94 and E2F targets from -2.72 to +2.02. |

### The chrY contrast is not on the same dosage axis

`ARPE19_chrY` (PRJNA930361) is our own realignment of a CRISPR chrY-elimination line, added
2026-09-29. It is kept apart from the other contrasts in one specific way.

A male carries **one** chrY, so losing it is a 1 → 0 event. Every other contrast here is
2 → 1 (or 3 → 2 in the aneuploid cancer lines) — a halving. The dosage check shows it
plainly: the lost-segment log2FC is **−6.26** for chrY against −0.35 to −0.84 everywhere
else, because the denominator goes to zero rather than halving.

So chrY is included in the pathway questions but **excluded from the size-scaling fit**,
where the x-axis assumes a comparable per-gene dosage step. It is drawn as a triangle in
panel c and marked `Y*`. Excluding it leaves the scaling answer unchanged (ρ = +0.50,
P = 0.27, n = 7).

What it shows on its own: everything modestly down, nothing up — ribosomal proteins −1.73
(FDR 1e-03), G2M −1.66, TNFα–NF-κB −1.62, MYC −1.58, E2F −1.53, p53 −1.28 (n.s.). The
p53/arrest program is absent, and its correlation with pure p53 induction is 0.00–0.03.

### How the two analyses fit together

The single-cell data catch cells a few days after the cut, before selection. The bulk data
come from clones that were picked and grown for weeks. Comparing them: the sharp p53-target
induction seen per cell at days 5-10 is largely gone in the surviving clones. CDKN1A, the
gene that enforces the arrest, is flat, and only RPS27L, ZMAT3 and MDM2 stay modestly up.

Growing a clone selects against arrested cells. So a dataset built from clones cannot tell
you what acutely happens after an arm is lost, and the GSE336056 p53/arrest program should
not be described as "the consequence of arm loss" - it is what freshly cut cells do before
they are selected away.

A separate result points the same way: in GSE150686, switching on p53 in **diploid** RPE1
raises the TP53 transcript (+0.87) but its targets do not respond (CDKN1A +0.17, not
significant; 1 differentially expressed gene in total). The same induction in **monosomic**
cells drives CDKN1A +1.89 (adjusted P = 5e-163). Monosomy is what licenses p53 to act.

## Contents

```
README.md               this file
REFERENCES.md           datasets and papers, with accessions and PMIDs
figures/                main figures (PDF + SVG + PNG) and the data behind each panel
code/
  gse336056/            single-cell analysis (Yq and 3p arms)
  arm_loss_meta/        bulk meta-analysis (download -> annotation -> DE -> meta -> figure)
tables/
  arm_loss_meta/        per-contrast DE and GSEA, QC, meta-analysis outputs
  prjna930361/          ARPE-19 chrY realignment: count matrix + chrY/PAR check
  gse336056/            single-cell DE and GSEA, arm/breakpoint tables
  sessionInfo.txt       R version and package versions used
```

## Pipeline

Bulk meta-analysis, run in order from `code/arm_loss_meta/`:

| Step | Script | What it does |
|---|---|---|
| 0 | `00_download.sh` | Fetches the three deposited count sets from GEO plus the hg38 cytoband table. Nothing is realigned. |
| 0 | `00_gene_annotation.R` | Builds gene -> symbol, chromosome, coordinates and **chromosome arm**. Arm boundaries come from the cytoband file's centromere bands, not typed-in numbers. |
| 1 | `01_build_datasets.R` | Reads the three datasets into one object and labels every sample. |
| 2 | `02_de_gsea.R` | DESeq2 for each contrast, then GSEA. Writes the cis-dosage check. |
| 3 | `03_meta.R` | The three questions above: shared program, p53 comparison, size scaling. |
| 4 | `04_figure.R` | The figure. |

Two design choices matter:

- **One clone counts once.** Replicates of the same clone are added together, so an arm
  represented by 3 clones contributes 3 data points, not 6. Control libraries are repeats
  of a single parent genotype and are kept as they are.
- **Genes on the lost segment are set aside before pathway scoring.** They are at half
  dosage by construction and are what defines the contrast, so leaving them in would
  guarantee a result. They are reported separately as a check (below) and excluded from
  the pathway ranks - the same rule used for chrY and 3p genes in the single-cell work.

### The check that the loss is real

Every loss contrast shows the expected drop for the genes on the lost segment, with the
rest of the genome flat:

| Contrast | Lost-segment log2FC | Rest of genome |
|---|---|---|
| RPTEC 3p, day 18 / day 90 | -0.82 / -0.84 | +0.10 / +0.08 |
| RPE1 monosomy 13 / 19p / 10+18 | -0.67 / -0.71 / -0.50 | +0.02 / -0.01 / +0.04 |
| Cancer lines (1q, 7p, 8q) | -0.35 to -0.55 | +0.01 to +0.07 |
| ARPE-19 chrY (hemizygous, 1 -> 0) | **-6.26** | -0.03 |

A cell with two copies that loses one would give -1. RPTEC, a near-diploid primary cell,
comes closest. The cancer lines sit higher because they carry more than two copies to begin
with, so losing one is a smaller proportional change - not because the loss is partial.
chrY is the exception and is bounded only by the pseudocount, because a male has one copy
and loses all of it (see the section above).

## Caveats

- **Deposited counts, not realigned.** Gene annotation differs between datasets, so
  comparisons are made within a dataset and only the resulting directions are compared
  across datasets.
- **Two contrasts are thin.** A2058 1q and HCT116 8q each rest on 2 clones and produced
  6,619 and 66 differentially expressed genes respectively. With 2 clones the variance
  estimate is unstable; treat both as weak.
- **Clones, not acute loss.** Every bulk contrast uses clones grown for weeks. See the
  early-vs-late result above.
- **Ribosomal protein finding is relative, not global.** Counts are normalised per library,
  which cannot detect a change in total output per cell. Showing that would need a
  spike-in or a per-cell RNA measurement.
- **The A2780 libraries** use a different annotation from the rest of GSE222379 (48,162 vs
  58,735 genes) and are handled as their own block.
- **The p53 comparison** uses the two monosomic induction contrasts. The diploid one is
  left out because the induction does not actually activate p53 there.
- **ARPE19_chrY is a weak contrast on two counts**, independent of the hemizygosity issue:
  it is one clone in three replicates, and the comparator is the untreated parental line,
  so Cas9 damage is not separated from Y dosage (the same gap as GSE336576). The chrY loss
  itself is solid - every male-specific Y gene reads 0 and PAR halves to 0.462 - but the
  downstream expression differences carry that caveat.
