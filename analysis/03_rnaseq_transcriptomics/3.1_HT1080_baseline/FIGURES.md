# Publication figure specification

Required for SEPT_11_2026_REANALYZE, following the user's request on 11 September 2026. Applies to main figures, individual plots and supplementary analyses. Nine main and six supplementary figures are generated in [outputs/figures](outputs/figures), with [final captions](outputs/captions/figure_captions.md) and a [source-data manifest](outputs/figures/figure_manifest.csv). The expanded scope includes all-collection enrichment/CAMERA, ZMYND8, chromatin, HDAC gene drivers and the direct public perturbation benchmark. Complex, addback and ATAC integrations remain unrun extensions.

The confirmed core design is repeated sampling of the same established CTRL and LOY populations across August and May, with one shared nucleofection origin. Plot six passage measurements per condition and identify the three within-period passage/day matches. Main figures emphasize effect sizes, trajectories and repeatability; ordinary gene-level significance stars, population-level confidence intervals and independent-validation labels are not part of the primary display.

## Visual standard

Use a restrained Nature-style layout: white background, black text, aligned panels, compact legends and consistent scales. Use 89 mm or 183 mm widths and a maximum 170 mm figure height. Ordinary lettering is 7 pt at final size, with 6 pt permitted for dense annotations and a 5 pt floor. Use Arial with embedded, editable text. [Nature figure preparation guide](https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/)

Panel labels are bold, upright, lowercase **a**, **b**, **c**, in 8 pt Arial, aligned at the upper left. Preserve vector text and line work. [Nature artwork guide](https://www.nature.com/documents/nature-final-artwork.pdf)

These are the project defaults. Set physical dimensions before composing a figure; split an overcrowded figure instead of reducing all lettering. Use approximately 0.5-0.7 pt axes and key lines, sentence-case axis labels, explicit units and short panel headings where necessary. Put the full descriptive title in the caption. Share axes and legends only when their variables, scales and mappings are identical.

Use `ggplot2` for panels, `patchwork` for composition, `ggrepel` for selected labels and `ComplexHeatmap` for matrix displays. Implement one shared theme/export helper and matching heatmap settings; do not independently restyle each notebook. The machine has ggplot2 4.0.2, patchwork 1.3.2, ggrepel 0.9.7, ComplexHeatmap 2.26.1, ragg 1.5.0, svglite 2.2.2 and colorspace 2.1.2. Arial and Cairo support were verified on R 4.5.2.

## Consistent encodings

| Biological group | Display label | Colour |
|---|---|---|
| CTRL | sgCTRL | `#0072B2` |
| LOY | LOY | `#D55E00` |
| PARENTAL | Parental | `#666666` |
| LOY_LXHA | LOY + empty vector | `#333333` |
| LOY_UTY_AB | LOY + UTY | `#009E73` |
| LOY_KDM5D_AB | LOY + KDM5D | `#CC79A7` |
| LOY_KDM5D_UTY_AB | LOY + UTY + KDM5D | `#E69F00` |

Use circles for August and triangles for May wherever batch needs a second encoding. Keep explicit group labels; colour alone must not carry essential information. Facet or add shapes when similar colours are hard to distinguish. Inspect the composed figures under protanopia/deuteranopia simulations and in grayscale.

Effect direction has a separate, labeled key: increased in the stated numerator = vermilion, decreased = blue, below the declared descriptive magnitude threshold = light gray. Do not label gray primary-analysis points "not significant". Missing/ineligible values are distinct from zero and from small effects. Use cividis for sequential abundance/density and a zero-centred blue-white-vermilion scale for signed effects, with symmetric limits. Preserve shared limits across directly compared panels. State and mark any display clipping; preserve the unclipped values in source data.

## Plot choices for this study

- Sample expression: show all observations with a clearly defined summary. Include passage-trajectory plots for the confirmed repeated CTRL and LOY population measurements, with passage order on the axis and population/batch identities visible. CTRL/LOY observations may also be connected within their confirmed batch-specific passage/day matches; captions must identify these as timing matches rather than independent biological pairs. Use inferential intervals only from a justified model and label within-population variability separately. Only connect parental/addback groups if their matching is established. Avoid fitted violin densities for groups of three.
- DE and batch agreement: MA plots, August-versus-May effect plots, and passage trajectories. Add volcano plots and effect/95% CI panels only where a justified inferential model supplies their statistics. State numerator and denominator and distinguish biological replication from repeated sampling of established populations. Label genes using a recorded selection rule. If points are too dense, rasterize only the point layer; keep axes and lettering vector.
- Heatmaps: show sample metadata, the gene-selection rule and the exact transformation. Label row z-scores separately from log2 fold changes. Display underlying sample-level abundance for central claims.
- Pathways and TFs: signed score/effect dot plots with gene/target coverage and statistical support; selected enrichment curves and leading-edge displays. Keep competitive tests, self-contained tests and TF target-pattern statistics distinguishable.
- Dosage/CN: align chromosome/segment coordinates and use matched effect scales where appropriate. Show the available DNA specimen count. Label count-derived paralog totals as relative abundance proxies.
- Rescue: plot effects with uncertainty and the actual comparison groups. Display factorial interactions explicitly. A nonsignificant residual difference is not labeled complete rescue.

## Planned figure set

| ID | Subject | Panels to assemble from saved outputs |
|---|---|---|
| Fig01 | Experimental design and QC | Sample/design summary; PCA; assignment/library QC; sample relationships |
| Fig02 | LOY expression differences and temporal repeatability | MA/effect plots; per-period agreement; passage trajectories; observed variability and influence |
| Fig03 | Pathways and regulatory responses | Pathway summary; leading-edge expression; TF activity; supporting targets |
| Fig04 | Y-linked expression, paralog dosage and CN context | Yp/Yq/PAR profiles; paralog abundance; pair totals; matched CN/RNA view when available |
| Fig05 | All MSigDB collections and CAMERA | Resource/eligibility counts; exploratory BH counts; method agreement; voom precision trend |
| Fig06 | ZMYND8 follow-up | ZMYND8/KDM5C/KDM5D RNA; original regulon scores; all eight target effects; MMP1/MMP3 edge-omission sensitivity |
| Fig07 | Chromatin-associated signatures | Source-defined enrichment; paired UP/DN response arms; reference-program overlap; gene-exclusion sensitivity |
| SuppFig04 | Chromatin regulator expression | All 34 mechanism-selected regulators across six passage matches; separate scale for UTY/KDM5D |

Panel membership and main-versus-supplementary placement depend on completed analyses and scope. Number the final set once the contents are fixed. Supplementary figures receive the same style and full captions, including filter/method sensitivities and additional QC. Do not imply that optional analyses have run.

## Exports and captions

For every completed figure, produce a vector PDF using Cairo, an editable SVG, and a 600 dpi PNG at the same physical dimensions. The PNG is a high-resolution project/report copy. Keep text editable in vector files and verify font embedding rather than assuming the export succeeded. Add journal-specific submission conversions only when preparing a submission; Nature distinguishes main and Extended Data formats. [Nature export guidance](https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/)

Each figure must have:

- A matching `<figure_id>.caption.md` and a combined `figure_captions.md` / manuscript-ready `figure_captions.docx`.
- Source data per panel in CSV/TSV, including plotted points, summaries, intervals and complete supporting statistics. Gene IDs remain available even when labels use symbols.
- A manifest entry linking panel IDs, output paths, source data, saved model/contrast, statistical family, selection rule and export settings.

Captions begin with a short descriptive title, then describe panels in letter order. Target fewer than 250 words, with exact sample numbers, replicate definitions, statistical tests and tails, multiplicity adjustment and definitions of error bars where applicable. Essential definitions stay in the caption; extensive numerical tables accompany it as source data. [Nature caption and statistical guidance](https://www.nature.com/nature/for-authors/initial-submission)

Specify the denominator for every n: libraries, independent cultures, matched pairs, DNA specimens, genes or targets. Report exact P or adjusted P values where readable; use source tables for high-dimensional results. Distinguish the zero-effect DE test from TREAT and distinguish TF target-pattern tests from biological-replicate tests. Explain expression transformations, colour/shape keys, dashed thresholds, abbreviations and any exclusions. Do not use unexplained stars or report a rounded value as P = 0.

Draft captions are in [FIGURE_CAPTIONS.md](FIGURE_CAPTIONS.md). Fill their fields from the sample manifest and saved results; never invent sample sizes or conclusions. Render captions beside figures in the report and export the same text separately for manuscript use. Keep captions outside the artwork itself.

## Final figure review

The HDAC step-2 addition contains Figure 8 (183 x 180 mm; four panels) and Supplementary Figure 5 (183 x 215 mm; three panels), using the same Arial, PDF/SVG and 600 dpi PNG standard. The first separates catalogue memberships, complete-source timing scores, post hoc readouts and gene influence. The supplement shows all libraries, separate response arms and fixed process annotations. The final captions are generated beside the figures; archived planning captions remain separate.

The step-3 addition contains Figure 9 (183 x 230 mm; four panels) and Supplementary Figure 6 (183 x 240 mm; three panels). These use full-page heights as documented exceptions to the 170 mm default, preserving legible 6-7 pt text for the complete contrast catalogue instead of omitting negative controls. Panel labels remain 8 pt. The 15-figure set has 54 panels. Figure 9 shows all primary contrasts and all active-transgene controls, all 54 enrichment tests and the same eight readouts. The supplement retains all modeled source libraries, published compatibility checks and all 27 contrasts under four fixed masks.

Before calling an output publication ready, inspect the composed figure at its intended print size, check labels/legends for overlap or clipping, verify colour accessibility and font embedding, and check PDF/SVG/PNG agreement. Confirm every panel matches its source data, the caption describes the actual plot, and all caption fields are resolved. Check that uncertainty and replication are represented accurately. Record the review status in the figure manifest; creating a high-DPI file alone does not complete this review.

The PRC2 extension adds Figure 10 (183 x 200 mm) and Figure 11 (183 x 240 mm), each with four panels. These documented full-page layouts retain 5.5-7 pt Arial text and 8 pt panel labels. Figure 10 shows complete PRC2 response arms, six timing matches, exhaustive EED/mark partitions and both directions of post hoc readouts. Figure 11 distinguishes all available promoter states, size-ineligible and unavailable combinations, comparator RNA arms, descriptive exclusions and regulator RNA. Signed numeric labels preserve heatmap direction in grayscale; source-arm and period symbols supplement color. Final project totals are 17 figures and 62 source-data panels.

The subsequent GSEA plot request adds Figure 12 and Figure 13 (both 183 x 240 mm, six panels each), plus all twelve panels individually at 89 x 78 mm. Each panel retains an exact running ES curve, complete membership barcode with colored leading-edge hits and the original ranked log2-effect track. All panels share ES limits of -0.70 to +0.50 and the identical ranked-effect scale. Original NES, represented count and global BH are shown separately from the running ES. Source context, response direction and the distinction between source chromatin marks and source RNA responses are explicit. The current composite totals are 19 figures and 74 panels. A dedicated export manifest tracks the 36 individual-panel PDF/SVG/PNG files.
