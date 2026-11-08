# Adopted realignment atlases — 28 October 2026

The author approved all three realignment atlases for Chapter 3 draft 2 and requested sex-chromosome karyotype/dosage context and conditional PAR expectations in their captions.

| Model | Current draft 2 | Historical central number | MSY regional median retention | Combined PAR1 / PAR2 medians |
|---|---|---|---|---|
| HT1080 | Fig 3.4 | Fig 3.10 | Yp 42.6% (67 bins); euchromatic Yq <1% (101 bins) | 67.2% / 53.2% |
| A549 | Fig 3.5 | Fig 3.12 | Yp 2.1% (67 bins); euchromatic Yq 2.7% (100 bins) | 66.8% / 63.7% |
| DU145 | Fig 3.6 | Fig 3.14 | Distal Yp 54.4% (23 bins); proximal Yp 6.1% (44 bins); euchromatic Yq 2.2% (101 bins) | 72.2% / 72.0% |

Every model has 24 PAR1 and four PAR2 bins. Regional figures show medians of paired-bin ratios, not ratios of regional read totals. Therefore whole-region CPM ratios need not equal the labels.

## Source and reproduction

Canonical renderer: [make_regional_y_atlas.py](../../code/make_regional_y_atlas.py). From the thesis root, run `python figures/chapter3/code/make_regional_y_atlas.py HT1080` (or `A549` / `DU145`). Use the project's plotting environment and bundled Open Sans fonts. `--output-dir` supports a separate review destination. The script reads the frozen per-model input directories here and saves derived tables and validation in `source_data/regional_y_atlases_ypar_20261028_rendered/<model>` under the output directory. It does not realign reads or refit copy number.

The three source directories freeze the approved broad 100-kb previews. All displayed DNA bin ratios, bin selections, regional medians, RNA observations, group means and mapped gene/probe positions were checked against those previews before adoption. Blue explanatory notes beneath a and c are absent in all three exports. DNA axes share 0–150% ticks with a small plotting margin; all points fit. Previous approved exports, captions, Results and renderer are preserved in [the adoption archive](../../../../archive/regional_y_ypar_adoption_20261028).

Adoption scripts and records: [atlas_adoption_20261028](../../../../research/model_development/atlas_adoption_20261028). Historical preview reports retain their original pre-adoption validations and must not be rerun as adoption checks: some deliberately assert that the formerly selected figures remain unchanged.

## Alignment, filtering and normalisation

The six BAMs were recovered to paired FASTQs and realigned to CHM13v2.0 with HG002-Y, masking only the Y copies of PAR. Their original primary read identities, mate assignments, sequences and qualities were preserved. Original BAMs were not overwritten. `LNCaP-sgCTRL` and `LNCaP-LOY-2X` are the legacy filenames of the DU145 samples. Extra drug-treated DU145 samples and A549-BC14 were not included in this realignment run.

For MSY, the reference mask requires ≥60% non-N sequence and ≥30% mappability; generic difficulty overlap is retained as annotation rather than a blanket exclusion. Both fitted samples must retain a bin, have finite corrected/segmented values, and sgCTRL corrected copy number must exceed 0.001. Displayed bins additionally require ≥50 primary MAPQ ≥37 control reads and must lie wholly within the indicated MSY region. No minimum depleted-sample read count or direction filter is applied. Centromere, Yq12 and region-boundary overlaps are excluded. Percentages use `100 * max(corrected_CN_LOY, 0) / corrected_CN_control`; tiny negative floating-point round-off is clipped to zero, with no positive numerator floor. Near-zero estimates are not a detection limit or proof of complete absence.

QDNAseq used `isPaired=TRUE` and MAPQ ≥37. Its inherited read parameters exclude unmapped, QC-failed and marked duplicate reads, but do **not** explicitly remove secondary/supplementary alignments. The additional support counts and PAR depth explicitly use primary reads only. Do not describe all QDNAseq counts as primary-only.

PAR uses primary-read aligned-base depth from the X copies, normalised by the same 21,446 strict-mask autosomal bins in each sample and then divided by sgCTRL. These bins require ≥20 control reads and a density ≥50 reads per 100 kb. PAR is not fitted using the old duplicated-reference mappability scores. PAR values combine homologous X- and Y-derived sequence; they are neither Y-specific copy number nor percentages of Y-negative cells. Terminal partial bins are retained when supported.

Zero-based half-open reference intervals: X PAR1 0–2,394,410; X PAR2 153,925,834–154,259,566; Y PAR1 0–2,458,320; Y PAR2 62,122,809–62,460,029. MSY Yp ends at 10,565,750, the centromere at 10,883,085 and euchromatic Yq at 27,449,937. Yq12 extends to 62,122,809. X-PAR positions are projected proportionally into labelled Y-PAR spans for display, not lifted over base by base. Only Yq12 is compressed. DU145's 4.80-Mb Yp split is an existing summary boundary, not a newly inferred breakpoint.

RNA data and sampling are unchanged: HT1080 has three passage samples per condition in each of two periods, normalised within period; A549 has four samples per condition; DU145 has three passage samples per condition. DU145 NLGN4Y relative expression remains undefined because the control mean is zero.

## Sex-chromosome dosage and expected PAR

An idealised control with `nX` X copies and `nY` Y copies would retain `nX/(nX+nY)` of combined PAR after complete loss of Y-derived PAR, provided X dosage and autosomal scaling remain unchanged. Thus XY → X predicts 50%, and XXY → XX predicts 66.7%. This calculation assumes that the relevant PAR copies track the chromosome dosage; rearrangements can invalidate that assumption.

The cultures were not cytogenetically karyotyped in this analysis. Captions distinguish inferred local dosage from a confirmed structural karyotype. Local checks used a fixed original strict reference mask: 1,178 non-PAR X bins, seven MSY Yp bins and 51 euchromatic Yq bins, with no sample-dependent bin selection. Primary-read depth was independently measured and normalised to the atlas's autosomal baseline; existing corrected copy-number estimates were inspected over exactly the same bins. These diagnostic bins differ from the broader display mask.

| Model | Local control X:Yq dosage | Depleted/control X dosage | Interpretation |
|---|---:|---:|---|
| HT1080 | 1.07:1 | 0.99–1.02 | Approximately XY-equivalent; nominal 50% PAR baseline. The empirical full-Y-loss baseline is about 51–53%. Retained Yp may elevate PAR1. |
| A549 | 1.91:1 | 1.00–1.01 | Approximately XXY-equivalent; nominal 67% PAR baseline. Empirical predictions from X/Yp or X/Yq are about 64–66%. |
| DU145 | 1.82–1.83:1 | 1.13–1.14 | Approximate XXY-equivalent dosage is useful as an idealised 67% baseline, but X dosage differs between pools. Using X/Yq and the observed X increase predicts about 73–74% under complete Y-PAR removal. |

DU145's Yp-based diagnostic prediction differs (approximately 66–69% after allowing for the X increase), illustrating that no single whole-chromosome model captures its regional/structural dosage. The caption identifies Yq explicitly as the baseline for the 73–74% calculation. Agreement with observed PAR does not prove complete Y-PAR deletion or identify its chromosome of origin. Nominal copy counts in aneuploid, heterogeneous cultures are models, not directly measured integer karyotypes.

Primary source context: [ATCC HT1080](https://www.atcc.org/products/ccl-121) reports a modal chromosome number of 46 with rearrangements; [ATCC A549](https://www.atcc.org/products/ccl-185) describes hypotriploidy and variable sex chromosomes in the initial stock, not a fixed XXY karyotype in this experiment. [ATCC DU145](https://www.atcc.org/products/htb-81) reports one X and a Y rearranged by translocation. Catalogue cytogenetics are not substituted for local dosage measurements. Consulted 28 October 2026.

Diagnostic source: `WGS/reanalysis/Y_PAR_masked_20261027/atlas_sex_dosage_20261028.py` and its output directory, alongside the earlier A549-only check. Copies of the summary, conditional predictions and provenance are frozen under `dosage_check/` here.

## Scope of adoption

Only the three regional atlases were replaced. Their captions and directly associated numerical Results claims were synchronised. Current draft 2 numbering (3.4–3.6), author edits elsewhere, and RNA analyses were preserved. Central captions retain their historical numbers and now state the draft 2 mapping. DOCX files were not regenerated.

The separate-track supplement remains the earlier unmasked-reference analysis with the old bin selection. Its caption now identifies that provenance instead of claiming that its bins are identical to the new main atlases. Other genome-wide and supplementary figures were not recalculated by this adoption.
