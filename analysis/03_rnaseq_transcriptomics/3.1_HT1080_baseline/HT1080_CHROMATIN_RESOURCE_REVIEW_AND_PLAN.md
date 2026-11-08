# HT1080 chromatin resources and integration plan

Reviewed 14 September 2026. Status: resource discovery and provenance audit completed; new chromatin-state assignments and RNA/ATAC integration have not been run.

We have enough directly measured HT1080 reference data to build useful cell-line-specific chromatin gene sets now. Start with the local GSE261834 NEED-seq panel, validate key marks against independent HT1080 ChIP-seq, and overlay the current RNA effects and our own ATAC effects. Keep the origin of each layer visible: reference chromatin describes another laboratory's HT1080 cells; our ATAC measures accessibility differences between our populations; RNA measures transcriptional differences. Neither reference chromatin nor ATAC alone measures a LOY-associated change in H3K27me3 or HDAC activity.

## Local resources recovered

The requested results folder contains **581 files (1.31 GiB)**, including 37 HTML reports, 131 CSVs, 19 RDS objects and publication figures. Its main ChIP-related analyses cover promoter co-marking, TSS signal profiles and permutations. A substantial ATAC workflow is nested under `ACTIVE_MOTIF_ATAC`.

The actual public source files are under `A:/multiome/data/ChIPseq`; the main local scripts are under `A:/multiome/scripts/ChIPseq`. Some rendered reports were produced on HPC and reference `/global/home/hpc5173/projects/LOY/ChIPseq/`. Their existence does not establish that the original tracks are available on this Windows machine.

| Resource | What is available locally | Reuse assessment |
|---|---|---|
| GSE261834 | 80 files including downloaded archive, extracted peaks/bigWigs, lamin BEDs, Hi-C valid pairs and two GEO metadata files | Best starting point for measured reference annotations |
| GSE135580 | Two GEO metadata files and rendered signal analyses | Recover/download signal tracks; not present in the inspected local dataset folder |
| GSE162690 | Metadata, old HT1080/multicell ATAC report | Recover tracks and establish their assembly |
| GSE181247 | Metadata and old confined-migration ATAC report | Perturbation context; recover tracks if used |
| GSE147672 | Brain metadata/report | Cross-tissue context only |
| GSE324952 | Mouse Ranbp2 ATAC report using mouse orthologs | Cross-species context only |
| ENCODE multitissue ATAC report | Survey of 103 tissues/cell types | Original accessions/builds are not recoverable from the displayed configuration alone |
| Own ATAC | Legacy Active Motif/nf-core results plus newer `projects/ATACseq/SEPT_08_2026_REPROCESS` | Reuse measurements/QC, revise interpretation for successive passages |

Machine-readable inventories and a 13-resource catalogue are in [provenance/chipseq_resources](provenance/chipseq_resources/). These are scoped inventories, not an exhaustive scan of every drive or HPC location.

## The measured HT1080 panel is broader than the old analysis used

The old configuration used only **GSM8152736 H3K4me3** and **GSM8152724 H3K27me3** for promoter classification. The local GEO metadata confirms both are HT1080. The same study provides the following standard panel, with local peaks and bigWig signal for every row:

| Mark | GEO samples labelled as replicates | Main proposed use |
|---|---|---|
| H3K4me3 | GSM8152736–GSM8152737 | Active/poised promoter context |
| H3K27ac | GSM8152734–GSM8152735 | Acetylated promoters and candidate active enhancers |
| H3K4me1 | GSM8152732–GSM8152733 | Candidate enhancer context |
| H3K27me3 | GSM8152724–GSM8152725 | Polycomb-associated repression |
| H2AK119ub1 | GSM8152720–GSM8152723 | PRC1-associated context alongside H3K27me3 |
| H3K9me3 | GSM8152726–GSM8152727 | Constitutive heterochromatin context |
| H4K20me3 | GSM8152728–GSM8152729 | Additional repressive chromatin context |
| H4K20me1 | GSM8152730–GSM8152731 | Additional chromatin context, interpreted with other marks |

These are measured NEED-seq profiles, not conventional ChIP-seq. Metadata specifies hg38 and RPKM-normalized bigWigs. Processed tracks are associated with the first accession for each mark; the paper describes replicate merging, so the number of local files must not be interpreted as the number of independent experiments. Confirm exactly which replicates feed each deposited track. Two metadata entries have a title/condition discrepancy (`H3K4me4` in the condition field); retain that flag until checked against antibodies/supplements. [GEO study record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261834), [NEED-seq paper](https://academic.oup.com/nar/article/53/2/gkae1317/7958221).

Additional local measurements include lamin B1/B2 profiles and domains, H3K9me2 under a size-selection protocol, c-FOS/c-JUN, NPAT, p53 and RB1. Hi-C valid pairs are present; EM-seq has three HT1080 sample records but methylation calls were not found locally. Separate standard NEED-seq, size-selection and MNase comparison protocols. The lamin/Hi-C analyses operate at domain scales and should not be interpreted as precise promoter or enhancer contacts.

H3K36me3 and H3K14ac are not in this standard local panel. In particular, H3K4me1 alone cannot establish the dual H3K4me1–H3K14ac recognition mechanism discussed for ZMYND8.

## Audit of the old promoter lists

The old workflow:

- Used EnsDb.Hsapiens.v86 gene-level TSSs, promoters extending 2.5 kb in each direction, and one retained entry per gene symbol.
- Converted zero-based peak starts to one-based coordinates correctly (`start + 1`).
- Required peak −log10(q) ≥5 for H3K4me3 and ≥2 for H3K27me3, with threshold sensitivity plots.
- Called a promoter “Bivalent” when its window overlapped both marks. The two peaks need not overlap each other within that 5-kb window.
- Exported **432 co-marked genes**, of which **209** were in the old 13,676-row RNA table.

Reuse the source peaks and signal. Rebuild the gene assignments using the current GENCODE v49 identifiers and fixed RNA universe. The current analysis has 25,461 representative ranked genes and a different statistical interpretation; the old DE/Fisher results must not be carried forward.

Call the reconstructed category **H3K4me3/H3K27me3 co-marked promoters** or **candidate bivalent promoters**. Separate bulk profiles cannot establish co-occupancy in the same cell or nucleosome. Replace “Unmarked” with “no peak detected under this definition”; low sensitivity, alternate TSS usage and missing coverage remain possible.

The old table labels CHRNB2, CHRNA7, MMP1, MMP3 and ZMYND8 as “Unmarked”. This is historical annotation, not evidence that those loci lack regulation. In particular, membership of CHRNB2 in a Mikkelsen GSEA leading edge does not establish H3K27me3 at its HT1080 promoter. That is a useful hypothesis for the rebuilt annotation and locus views to test.

## Public resources to add

**Independent measured chromatin.** GSE135580 includes two HT1080 ChIP replicates for H3K27me3, H3K9me2, H3K9me3, HP1α and CTCF, with inputs; it also offers Protect-seq and replication timing. Use it to distinguish Polycomb-associated genes from other inaccessible domains. Some processed tracks/domains are coarse: verify resolution before promoter scoring. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135580), [Protect-seq paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7026637/).

**A second H3K27me3/H3K27ac reference.** GSE153869 includes WT HT1080 as well as DDIT3 and fusion-expressing derivatives. WT H3K27me3 samples are GSM4656285/6291, H3K27ac samples GSM4656288/6290, and CTCF samples GSM4656284/6289, with WT inputs. Processed BEDs and raw data are public; metadata specifies hg38. Isolate WT and verify controls/replicates before integration. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153869), [WT H3K27me3 sample](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4656291).

**Active transcription and accessibility.** GSE86504 provides HT1080 H3K4me3, MYC, MEN1 and Pol II ChIP with associated knockdown RNA. Its H3K4me3 coordinates are hg19 and generic shRNA metadata needs resolution. GSE162690 contains real HT1080 bulk ATAC samples GSM4957247/7248, despite also containing other cell lines. Its GEO build is hg19, whereas the old rendered analysis uses hg38 EnsDb v86 TSSs. This is an unresolved build risk: the omitted HPC configuration may have supplied converted/reprocessed tracks, so an actual mismatch has not been established. [GSE86504](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86504), [H3K4me3 sample](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2304703), [ArchR study](https://www.nature.com/articles/s41588-021-00790-6).

**A useful migration comparison.** GSE181247 measures HT1080 accessibility under different collagen/confinement conditions. This could contextualize the ECM and migration programs, but represents a different perturbation. [Primary paper](https://www.sciencedirect.com/science/article/pii/S2589004222012500).

**Ready-made states, with a provenance qualification.** ENCODE confirms that ENCFF461DZD / ENCSR295VRM is a released GRCh38 HT1080 18-state ChromHMM annotation, from EpiMap biosample BSS00709. EpiMap uses a mixture of observed and imputed inputs. The exact observed/imputed status of its HT1080 marks remains to be traced; the file has no `derived_from` list. Use this as a model-derived secondary annotation and compare it with measured marks. [ENCODE annotation](https://www.encodeproject.org/annotations/ENCSR295VRM/), [EpiMap primary paper](https://www.nature.com/articles/s41586-020-03145-z).

A 2025 HT1080 paper cites this segmentation and six ENCODE signal files. Our direct API audit found that **all six signal files are Avocado imputations**, not new HT1080 measurements. Also, **ENCFF408PEC is H3K9ac**, although the paper describes it as H3K9me3. Preserve the correct accession/target pairing and do not count these predictions as independent experimental validation. Details are in `ENCODE_provenance_audit.csv` and cached JSON. [Paper](https://academic.oup.com/nar/article/53/14/gkaf720/8219897), [H3K9ac file](https://www.encodeproject.org/files/ENCFF408PEC/), [imputed H3K27me3 file](https://www.encodeproject.org/files/ENCFF633FBI/).

**A particularly relevant HDAC/CoREST lead.** An April 2026 preprint studies HT1080 treated with romidepsin, the LSD1 inhibitor GSK-2879552, or both, using NicE-seq, NEED-seq and RNA analyses. Reported features include H3/H3K4me1/H3K4me2/H3K9ac, RCOR1, RUNX and JunB. This could supply an HT1080 perturbation benchmark for our HDAC and neural-program findings. However, the inspected full text states that data are in GEO without providing an accession, and gives differing LSD1 exposure durations in Results and Methods. Treat it as an acquisition lead, not an analysis-ready resource. Drug responses would not isolate HDAC1 versus HDAC2 catalytic activity. [Preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC13081975/), [publication status](https://pubmed.ncbi.nlm.nih.gov/41993470/).

For ZMYND8/KDM5D, GSE82260 remains a useful binding/perturbation reference from **DU145**, not HT1080. This audit did not verify a ready-to-use, genome-wide HT1080 HDAC1/2, EZH2/SUZ12 or ZMYND8/KDM5D occupancy panel. That is a search result, not proof that none exists. Targeted HT1080 hTERT ChIP literature supports mechanistic possibilities but cannot supply a genome-wide target list. [GSE82260](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82260), [HT1080 HDAC/REST/LSD1 promoter study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5602382/).

## Our ATAC is useful, but its design must follow the new clarification

**User clarification, 14 September 2026:** CTRL/LOY_11, _12 and _21 are successive passages; their direct RNA collection partners are not completely known.

The newer September 8 ATAC workflow already contains fresh peak calling, QC, differential accessibility, csaw, motif analyses and footprint predictions. It uses the same six libraries as older processing variants; these are not additional biological replicates. Its documentation assumed independent replicates and integrates an older RNA model. The current clarification supersedes that assumption for this integration.

Reuse verified counts, intervals and technical QC. Inspect sample effects and passage sensitivity, then report descriptive condition differences with the same population-level limitation as RNA. Do not carry the older ATAC or RNA significance labels into a new joint causal claim. Do not construct sample-wise RNA/ATAC correlations, paired multimodal tests or invented collection dates. Group-effect integration is possible without exact pairs, provided the temporal mismatch remains explicit. A passage covariate does not create independent LOY-generation experiments.

The existing ATAC pipeline's strong sensitivity of significance to passage adjustment reinforces the need for this revision. Copy-number context can be added where build and specimen relevance are documented; D0 WGS is not automatically a sample-matched correction for later passages. Existing source ATAC reports are historical and were not changed during this audit.

## Recommended work order

1. **Build a measured HT1080 promoter/region atlas.** Freeze the eight standard NEED-seq inputs, provenance, signal scales and controls; record checksums and replicate pooling. Use the current hg38/GENCODE v49 annotation and stable IDs. Keep TSS-level assignments before summarizing to genes. Use a prespecified primary promoter window and report ±1 kb, ±2 kb and the old ±2.5 kb as sensitivity, rather than choosing whichever produces enrichment. Use broad-domain coverage alongside peak overlap for H3K27me3/H3K9me3. Do not threshold RPKM signal identically across antibodies. Export all eligible genes, including explicit unknown/no-peak flags, before looking at RNA direction.

2. **Define a compact, interpretable set collection.** Include H3K27me3-positive, H2AK119ub1-positive and jointly positive promoters; H3K4me3/H3K27me3 co-marked promoters; H3K4me3/H3K27ac active promoters; H3K9me3/H4K20me3 repressive domains; and lamin-associated genes. Keep individual mark features alongside derived categories. Distal H3K4me1/H3K27ac regions are candidate enhancers, with gene links labelled by evidence. Baseline Hi-C compartments are not validated enhancer–gene links. Export a feature matrix, gene sets in GMT format and full membership/provenance tables.

3. **Connect to the current RNA findings.** Examine full effect distributions, August/May direction and passage consistency. Test the prespecified custom sets with the existing fgsea ranking and CAMERA framework, retaining fixed- and estimated-correlation results as sensitivity analyses. Use a separately labelled BH family for the new custom collection; never replace the original global GSEA q-values. All statistics remain conditional on these established populations. Report gene-set overlap, leading-edge dependence and abundance/CpG/gene-length/CN sensitivity. Do not generate the sets from selected RNA leading edges and then test them against the same ranking.

4. **Validate reference annotations and integrate own ATAC.** Prioritize independent H3K27me3 support from GSE135580 and WT GSE153869. Compare EpiMap as a secondary model-derived layer. Join own ATAC promoter effects to current RNA at gene/group level; retain discordant genes. A reference H3K27me3-positive promoter that opens and rises in RNA is compatible with derepression, but does not measure H3K27me3 loss. A gene with reduced accessibility and RNA in a reference acetylated region suggests loss of regulatory activity without identifying the responsible enzyme.

5. **Resolve the mechanistic branches.** Ask whether Mikkelsen leading-edge genes actually have measured HT1080 Polycomb marks; whether the downregulated EED-response genes lie in the same or different states; whether HDAC-associated RNA changes localize to acetylated/CoREST-associated regions; and whether IFN/ECM changes remain after chromosomal-domain sensitivity. Review CHRNB2/CHRNA7, MMP1/MMP3 and representative IFN/ECM loci in context. Pursue the HT1080 HDAC/LSD1 dataset if its accession becomes available; keep drug-response and direct-binding gene sets distinct.

Use `GenomicRanges`/`rtracklayer` for assembly-aware overlap and bigWig extraction, the pinned GTF for promoter definitions, and `deepTools` or equivalent for signal profiles. Reuse existing fgsea/limma machinery. A de novo ChromHMM model is a later option if the measured panel, controls and replicate structure support it; it is not required to make useful sets, and the absent H3K36me3 track limits standard transcription-state resolution.

## Planned publication figures and caption requirements

- **Reference atlas:** measured mark heatmaps at prespecified promoters, with RNA effects and source-study support as annotations. Caption identifies source cultures, assays, windows, units, ordering and no-peak handling.
- **Chromatin-state/RNA relationships:** enrichment curves plus effect-distribution panels for the frozen custom sets, with exact eligible set sizes, multiple-testing family, CAMERA variant and repeated-population scope.
- **RNA–ATAC integration:** condition-level effect scatter and selected passage profiles, with no implied sample matching. Caption states the six ATAC libraries are successive passages and the RNA partners are uncertain.
- **Locus views:** measured baseline histone tracks, own CTRL/LOY ATAC tracks and GENCODE transcripts for selected genes. Distinguish measured and imputed tracks visually; label normalization, pooling and copy-number context. No accessibility-to-histone-change inference in the caption.

Export vector PDF/SVG and 600-dpi PNG with readable type, consistent axes and restrained colours, accompanied by standalone captions and plotted-source tables. These figures are planned, not generated by this audit.

**First concrete implementation:** the measured GSE261834 feature matrix and compact gene-set collection. This can proceed without RNA–ATAC sample matching and directly tests whether the current epigenetic GSEA interpretation transfers to measured HT1080 chromatin.

## Audit reproducibility

`audit_chromatin_resources.py` inventories sources and extracts legacy report text; `--upstream-only` refreshes the upstream inventory/metadata without extracting HTML again. `fetch_chromatin_metadata.py` caches small public ENCODE/article records. `catalogue_chromatin_resources.py` builds the curated tables. All outputs from these scripts are confined to `provenance/chipseq_resources/`. Counts, source peaks, source scripts, older ATAC analyses and completed RNA/GSEA results were not modified. Inventory creation validates discovery and metadata consistency; it does not substitute for alignment/peak quality review before quantitative reuse.
