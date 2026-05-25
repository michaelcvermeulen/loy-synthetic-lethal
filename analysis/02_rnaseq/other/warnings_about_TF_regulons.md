# Critical Assessment of decoupleR ULM + CollecTRI for TF Activity Inference on Bulk RNA-seq, with Specific Reference to Small-N Isogenic (LOY vs Y+) Designs

## TL;DR
- **decoupleR-ULM with CollecTRI is a reasonable hypothesis-generating tool on bulk RNA-seq, but for a 3 vs 3 isogenic LOY comparison it is structurally a re-weighted differential-expression signal and is empirically noisy at small N — it should not be treated as orthogonal evidence of TF protein activity without independent validation (ChIP-seq, ATAC-seq, motif enrichment, or perturbation).**
- **The dominant failure modes that matter for the user's design are: (i) extreme correlation/non-independence among TF activity scores driven by overlapping regulons; (ii) sensitivity to library size, normalization, and batch through the underlying gene-level statistic; (iii) confounding by cell-cycle, stress, interferon, and senescence programs (which is especially concerning for LOY because Y-loss is associated with a DNA-damage/injury response signature, per Wilson et al. 2024 *Genome Biology*); and (iv) one well-conducted independent benchmark (Trescher & Leser, 2019, *Scientific Reports*, DOI: 10.1038/s41598-019-46053-7) showed that across 54 TF knockdown experiments, the knocked-down TF ranked in the top 5% in only 15/54 cases, leading the authors to conclude verbatim that "the investigated methods cannot yield robust TF activity estimates in knockdown schemes."**
- **Best-practice safeguards: feed the t- or Wald-statistic (not normalized counts) from a properly modeled DE analysis (DESeq2/limma with batch covariates) into ULM; use the consensus score across methods (ULM + MLM + wsum) as recommended by Badia-i-Mompel et al. 2022; require minsize ≥ 5 targets; correct for multiple testing across TFs; report and interrogate correlation structure among "significant" TFs; and treat any TF call as a hypothesis to be validated by ATAC-seq motif enrichment on Y-loss–matched samples, or by knockdown/rescue.**

---

## Key Findings

1. **CollecTRI + ULM is the current saezlab-recommended default**, replacing DoRothEA, based on the benchmark in Müller-Dott et al., 2023, *Nucleic Acids Research* (DOI: 10.1093/nar/gkad841), which reported improved AUROC/AUPRC against four other meta-resources (DoRothEA, ChEA3, RegNetwork, Pathway Commons) using the KnockTF perturbation database as ground truth. CollecTRI contains "43,175 signed TF–gene interactions for 1186 TFs" (Müller-Dott et al. 2023, PMC10639077).

2. **ULM is mathematically a weighted regression of gene-level statistics on the TF–target sign vector.** The decoupleR Bioconductor vignette is explicit: "the obtained t-value of the slope is the score." This means the TF score is, by construction, a linear summary of the input gene statistics weighted by regulon membership and sign — it is **not** an independent measurement of protein activity and inherits all noise, batch, and library-size properties of the upstream DE pipeline.

3. **Sample-size guidance is weak in the official literature.** The decoupleR bulk vignette demonstrates the workflow on the standard 3 WT vs 3 KO FOXA2 pancreatic-cancer dataset, implicitly endorsing 3 vs 3 as workable, but no formal power analysis or false-positive calibration is provided. The Holland et al. 2020 benchmark (*Genome Biology*, DOI: 10.1186/s13059-020-1949-z) reported DoRothEA AUROC of 0.752 on the Perturb-seq (7d) dataset at full gene coverage, with performance degrading specifically below ~2000 covered genes per profile (SCENIC AUROC 0.631 and D-AUCell AUROC 0.629 trended toward near-random in low-coverage regimes); variance in TF-activity ROC curves shrinks as both library depth and sample number increase — i.e., small-N bulk comparisons sit on the noisy end of the operating curve.

4. **Regulon overlap is a real, formally proven problem for univariate methods like ULM.** The EPEE paper (Madar et al., *Bioinformatics* 2019, DOI: 10.1093/bioinformatics/btz398) showed that in context-specific FANTOM5 networks, "even the least complex context has overlaps among 96.4% of all the TF-TF pairs, whereas the median was 98%," and proved by causal-inference argument that univariate inference cannot disentangle co-regulated TFs — when one TF is genuinely perturbed, its co-regulating peers appear correlated and significant. CollecTRI partially mitigates this with sparser regulons (mean Jaccard 0.01 vs other GRNs) but does not eliminate intra-CollecTRI overlap.

5. **Independent benchmarking is critical of these methods.** Trescher & Leser (2019, *Scientific Reports*, DOI: 10.1038/s41598-019-46053-7) tested major regulon-based TF inference tools on 54 TF knockdown experiments and found the actually-knocked-down TF ranked in the top 5% in only **15 of 54 cases**, with poor agreement across methods. The newer TFActProfiler benchmark (Sept 2025 bioRxiv, DOI: 10.1101/2025.10.05.680506) shows CollecTRI/DoRothEA outperform broad unfiltered priors (ChIP-Atlas, CellOracle base) but are themselves outperformed by uncertainty-aware models, and explicitly flags "small case-control bulk RNA-seq studies without matched perturbation datasets" as a regime where existing methods underperform. Hu et al. (2025) *Advanced Science* (DOI: 10.1002/advs.202410745) reported that across 109 TF perturbations, average AUROC was 0.69 — i.e., ~30% of TFs at or near random performance — and DoRothEA/SCENIC were outperformed in 80/109 and 66/109 TFs respectively.

6. **For the LOY-specific context, two known confounders should worry the user.** (a) Wilson et al. 2024 *Genome Biology* (DOI: 10.1186/s13059-024-03173-2) studied LOY in kidney proximal tubule cells and concluded: "LOY is an indicator of increased DNA damage and potential marker of cellular senescence that can be applied to single-cell datasets in other tissues," and that "LOY was associated with decreased expression of DNA repair and cell cycle genes across multiple cell types—suggesting that LOY is associated with an injury response in multiple cell states." This will drive E2F, TP53, RB-related, and stress-response TF signals that look biological but are downstream of genomic instability rather than a specific Y-chromosome perturbation. (b) Bladder-cancer LOY work (Abdel-Hafiz et al. 2023 *Nature*) showed altered TGF-β and EMT-related programs; SMAD2/3/4 and SNAI/ZEB are heavily represented in CollecTRI regulons.

---

## Details

### CATEGORY 1 — Application papers using decoupleR-ULM (or DoRothEA) with CollecTRI on bulk RNA-seq

The CollecTRI/ULM combination is recent (CollecTRI itself was published November 2023), so the corpus of high-impact application papers using **exactly** this configuration on bulk RNA-seq is small and still emerging. The strongest examples identified:

- **Müller-Dott et al., 2023, *Nucleic Acids Research*** (DOI: 10.1093/nar/gkad841) — the CollecTRI paper itself includes a case study applying CollecTRI + ULM to bulk RNA-seq from three CPTAC cancer types (UCEC, LUAD, CCRCC) and benchmarks against KnockTF perturbation data. The most authoritative demonstration but methodological, not biological discovery. **Caveat:** KnockTF overlaps with the literature curations that built CollecTRI — the TFActProfiler authors explicitly excluded KnockTF for this reason, calling the overlap "data leakage." CollecTRI's reported AUROC is therefore likely optimistic.

- **Baudre et al., 2025, *Cancer Research*** (DOI: 10.1158/0008-5472.CAN-25-0995) — TNBC drug-tolerant persisters across multiple PDX models. Methods state: "Regulon activity was inferred using DecoupleR (v2.10.0)… a consensus score of TF enrichment scores was obtained by combining 3 methods (ulm, wsum and mlm) on the T-value from the differential expression analysis obtained with limma." Critically, they used the **consensus** score, not ULM alone, and the **t-value**, not normalized counts. TF activity is load-bearing; orthogonal validation comes from paired bulk epigenomics on the same persister populations.

- **Afenteva et al., 2025, *Cell Reports Medicine*** (DOI: 10.1016/j.xcrm.2025.102316) — chemo-refractory HGSC, n=58 patients (38 sensitive vs 20 refractory) bulk RNA-seq deconvolved by PRISM. Methods: "univariate linear model method executed via decoupler.run_ulm()… the stat-value of the DESeq2 output as input." Identified attenuation of IRF1/2/9 as driver of refractoriness. **Load-bearing finding**, validated by paired scRNA-seq and functional IFN-I experiments. Note: cohort scale (n=58) is far larger than a typical isogenic 3v3.

- **Cambridge CRUK PDAC study (Lin et al.), 2024–2025** (bioRxiv DOI: 10.1101/2024.04.28.591518; GEO GSE274684) — KPC mouse PDAC, Smad4 WT vs KO bulk RNA-seq, used ULM with CollecTRI on the DESeq2 *stat* column; identified JAK/STAT TF family as activated in Smad4-KO. Load-bearing; validated by in vivo JAK/STAT inhibitor sensitivity and paired scRNA-seq.

- **Foundational DoRothEA + VIPER applications (pre-CollecTRI era):**
  - **Garcia-Alonso et al., 2019, *Genome Research*** (DOI: 10.1101/gr.240663.118) — the DoRothEA benchmark paper itself, applied DoRothEA+VIPER to TCGA bulk RNA-seq and CCLE/GDSC, validated against ChIP-seq and TF essentiality (DEMETER/ATARiS) screens. The standard for "good practice" in this field.
  - **Holland, Szalai, Saez-Rodriguez, 2020, *BBA Gene Regulatory Mechanisms*** (DOI: 10.1016/j.bbagrm.2019.194431) — mouse DoRothEA, perturbation-benchmarked.

**Pattern across well-conducted applications:** they (i) use a DE *statistic* (t- or Wald-stat) as input, not raw counts; (ii) use the **consensus** of multiple decoupleR methods, not ULM alone; (iii) treat TF activity as a hypothesis-generation step paired with orthogonal data (ATAC, ChIP, scRNA-seq, or functional inhibitor experiments); and (iv) have substantially larger sample sizes (n ≥ ~10/group, often cohort-scale) than a 3v3 isogenic design.

### CATEGORY 2 — Benchmarking, validation, and critique literature

#### 2A. The foundational benchmarks (saezlab, friendly)

- **Garcia-Alonso et al., 2019, *Genome Research*** (DOI: 10.1101/gr.240663.118): DoRothEA outperforms other resources on TF perturbation datasets when coupled with VIPER. AUROC for confidence-A interactions is moderately good (~0.7) on cancer cell perturbation benchmarks; confidence levels D/E perform near random. **Implication for user: only the highest-confidence subset of any regulon is reliable.**

- **Holland et al., 2020, *Genome Biology*** (DOI: 10.1186/s13059-020-1949-z): Established that performance is **more sensitive to the choice of gene set/regulon than to the statistical method**. DoRothEA AUROC was 0.752 on Perturb-seq (7d) at full coverage, degrading below ~2000 covered genes. Coupling DoRothEA with simpler scores (AUCell, ULM) gave broadly similar performance to VIPER. *Critical implication:* the ULM vs VIPER vs wsum vs MLM choice is not the dominant axis of variation — the regulon definition and the input gene-level statistic are.

- **Badia-i-Mompel et al., 2022, *Bioinformatics Advances*** (DOI: 10.1093/bioadv/vbac016): The decoupleR software paper benchmarks ULM/MLM/wsum/wmean/GSVA/AUCell/VIPER on transcriptomic and phosphoproteomic perturbation data. Headline conclusion (verbatim): "simple linear models and the consensus score across top methods perform better than other methods at predicting perturbed regulators." ULM and MLM are both endorsed; the **consensus** is the recommended default — not ULM in isolation.

- **Müller-Dott et al., 2023, *NAR*** (DOI: 10.1093/nar/gkad841): CollecTRI outperformed DoRothEA, ChEA3, RegNetwork and Pathway Commons on KnockTF; explicitly discusses TF overlap and notes that calculated activity scores can correlate with **regulon size** — a confound the user must check (very small regulons produce unstable scores; very large ones risk picking up generic transcriptional state).

#### 2B. Independent / critical benchmarks (this is where the bad news lives)

- **Trescher & Leser, 2019, *Scientific Reports*** (DOI: 10.1038/s41598-019-46053-7; PMC6610105): the most damaging independent benchmark. Tested major TF-activity inference tools on 54 TF knockdown experiments. **Only 15/54 KD TFs ranked in the top 5%.** Tools disagreed substantially. Pathway-/network-neighbors of the KD TF were not ranked above random. Verbatim conclusion: "the investigated methods cannot yield robust TF activity estimates in knockdown schemes." This applies directly to the user's regime (a defined perturbation in a clean isogenic system).

- **Kang et al. ("Inferring TF activities… with constraints from TF perturbation data"), 2021, *Bioinformatics*** (DOI: 10.1093/bioinformatics/btab019; PMC8189679) — argues that almost all validation of TFA methods has been "small-scale and often qualitative, highlighting successes with just a few TFs," and that reasonable inference requires perturbation-derived constraints, which off-the-shelf ULM does not have.

- **Madar/EPEE et al., 2019, *Bioinformatics*** (DOI: 10.1093/bioinformatics/btz398) — formal causal-inference argument that overlapping regulons render univariate TF estimates **dependent on each other** by construction (>96% TF-pair overlap median in FANTOM5 contexts), and that activity must be inferred jointly (multivariate). This is a direct theoretical critique of ULM: when several CollecTRI TFs share targets (extremely common), their ULM scores are mathematically coupled and a single underlying transcriptional change will light up multiple correlated TFs. EPEE benchmarks show ISMARA loses performance with TF overlap, while EPEE (multivariate) does not. ULM is closer to ISMARA than to EPEE in this respect.

- **TFActProfiler, October 2025 (bioRxiv DOI: 10.1101/2025.10.05.680506)** — most recent independent benchmark; CollecTRI/DoRothEA outperform broad unfiltered priors but are themselves outperformed by uncertainty-aware models. Explicitly notes that "regression-based models often require substantial task-specific training data and can sometimes fail to outperform simple baselines, limiting their utility in common scenarios, such as small case-control bulk RNA-seq studies." This is *exactly* the user's scenario being flagged.

- **MetaTF, Hu et al., 2025, *Advanced Science*** (DOI: 10.1002/advs.202410745) — across 109 TF perturbations, average AUROC was 0.69, and DoRothEA/SCENIC were outperformed in 80/109 and 66/109 TFs respectively. Reinforces that even the best regulon-based methods have substantial per-TF noise.

#### 2C. Confounder-specific evidence relevant to LOY

- **Cell-cycle / proliferation confound:** TF activity profiles split robustly into "stress response" and "proliferation/homeostasis" modules (e.g., intronic TF-activity work in *PLoS Computational Biology*, PMC8782462). E2F1–E2F8, MYC, FOXM1, MYBL2 are major CollecTRI TFs and respond to proliferation rate per se. Any difference in growth rate between Y+ and Y- cells will produce E2F/MYC/FOXM1 "activity" differences with no specific causal role. **Recommendation:** measure proliferation explicitly (Ki-67, growth curve, EdU) and either regress it out or restrict interpretation.

- **DNA damage / senescence / injury confound (most important for LOY):** Wilson et al. 2024 *Genome Biology* (DOI: 10.1186/s13059-024-03173-2) studied LOY in kidney proximal tubule cells and reported that "LOY was associated with decreased expression of DNA repair and cell cycle genes across multiple cell types—suggesting that LOY is associated with an injury response in multiple cell states," concluding that "LOY is an indicator of increased DNA damage and potential marker of cellular senescence that can be applied to single-cell datasets in other tissues." Note this finding is from proximal tubule cells; the cross–cell-type generality is the authors' inference and should not be assumed for the user's specific cell line. Nevertheless, this directly drives TP53, NFKB, ATM, E2F, and senescence-associated secretory phenotype TF signals in CollecTRI.

- **Interferon / EMT confounds:** Afenteva et al. 2025 above showed IRF1/2/9 activity is highly informative but also highly sensitive to deconvolution and immune contamination. Bladder LOY studies (Abdel-Hafiz et al., *Nature* 2023) show altered TGF-β/EMT — SMAD2/3/4 and SNAI/ZEB are CollecTRI TFs that will move.

#### 2D. Library size, normalization, and batch sensitivity

- The decoupleR documentation requires "normalized log-transformed counts" or a DE statistic as input. The score *is* a linear function of those inputs, so **any unmodeled batch or library-size effect in the input propagates directly into the TF score.** Holland et al. 2020 demonstrate clearly that AUROC variance shrinks with deeper libraries and more cells/samples; the inverse is also true.

- **Strong recommendation from saezlab's own vignettes:** the bulk vignettes recommend using DESeq2/limma t-stat or Wald-stat as input rather than per-sample normalized counts when doing a contrast — because the DE step explicitly models library size, sample-level variance, and (if specified) batch covariates. Feeding `vst`/`rlog` counts per sample into ULM and then comparing sample-level scores by t-test is more sensitive to nuisance variance than running DE first and feeding the DE statistic in once.

#### 2E. Multiple testing and significance interpretation

- ULM produces a p-value per TF per sample (or per contrast). With ~1,186 TFs in CollecTRI, naive p < 0.05 will return ~59 false positives per analysis. The Müller-Dott paper and the saezlab tutorials all use multiple-testing correction; the user must do the same (BH-FDR across TFs). Additionally, because regulons overlap, the BH assumption of weak dependence is **violated** — empirical FDR is unknown and likely worse than nominal. Permutation-based empirical FDR using shuffled condition labels is the defensible answer.

#### 2F. Is ULM-CollecTRI just a re-weighted DEG signal?

- Mechanistically: substantially yes. ULM regresses gene-level statistics on the signed regulon vector; the score is a weighted sum of DE statistics restricted to regulon members. This is acknowledged implicitly by every saezlab tutorial showing TF activity computed *from* the limma/DESeq2 output. The added value is (i) directionality (signed mode-of-regulation weights) and (ii) collapsing thousands of gene tests into ~1,186 interpretable TF tests with statistical power gain from aggregation. But it does **not** measure protein activity, nuclear localization, post-translational modification, cofactor availability, or DNA-binding — all of which are real determinants of TF activity that are invisible to mRNA-target footprinting. Dugourd & Saez-Rodriguez (2019, *Curr Opin Syst Biol*, DOI: 10.1016/j.coisb.2019.04.002) acknowledge this is a *footprint* method, i.e., it observes the consequences of TF activity rather than the activity itself, and is only as reliable as the regulon's correspondence to context-relevant targets.

---

## Recommendations

**Stage 1 — Before running decoupleR at all:**
1. Run rigorous QC. Confirm Y+ vs Y- samples cluster primarily on Y status (not on batch, library size, or growth state) in PCA after vst normalization. If they don't, decoupleR cannot fix it.
2. Quantify proliferation/cell-cycle explicitly (Ki-67, EdU, cell-cycle gene-signature score, doubling time). Include cell-cycle scores or proliferation rate as a covariate in the DE model if there is any group difference.
3. Include batch as a covariate in DESeq2/limma. Apply ComBat-seq if multiple sequencing batches exist.
4. Screen for senescence/DNA-damage signatures (e.g., CDKN1A/p21, CDKN2A/p16, γH2AX target panel, SASP signature, p53-target signature) — given the Wilson et al. 2024 finding, expect these to be elevated in Y- cells; this is itself biology but it will confound downstream TF inference.

**Stage 2 — Run decoupleR correctly:**
5. Use the **DESeq2/limma t- or Wald-statistic from the contrast as the single input matrix** (one column = the contrast), not per-sample normalized counts. This is the saezlab-recommended path for differential TF activity and lets the DE model handle library size/batch.
6. Run **the consensus** across at least ULM + MLM + wsum (per Badia-i-Mompel et al. 2022), not ULM alone. Disagreement among methods is a red flag for that TF.
7. Use `minsize = 5` (default) or higher. Discard TFs with <10 expressed targets in your data, and report regulon size alongside scores (Müller-Dott et al. 2023 explicitly tested regulon-size confounding).
8. Apply BH-FDR across TFs; treat FDR < 0.05 with extreme caution because of overlap-induced dependence; consider permutation-based empirical FDR using shuffled condition labels.

**Stage 3 — Interpret defensively:**
9. Inspect the **target-level evidence** for each significant TF (decoupleR provides `plot_targets`) — does a small handful of strongly DE genes drive the entire signal? Is it just MYC targets that are also E2F targets that are also FOXM1 targets?
10. Cluster significant TFs by target-set Jaccard similarity. If many TFs cluster, treat them as a single regulatory module, not as independent hits.
11. Cross-check whether the "winners" line up with known cell-cycle/IFN/senescence/EMT/DNA-damage modules. If so, downweight them and look for the residual — what remains after stripping the genomic-instability/injury signature?

**Stage 4 — Validate orthogonally before claiming:**
12. Run ATAC-seq on the same isogenic Y+/Y- lines if possible — chromVAR motif accessibility scores against JASPAR/HOCOMOCO motifs are an orthogonal readout that does not share the regulon assumption. Concordance between an ULM-significant TF and its motif accessibility is much stronger than either alone.
13. Confirm at least one top TF call via siRNA/CRISPR knockdown and rescue of the Y-loss phenotype, or via ChIP-qPCR at known targets.
14. Replicate in an independent isogenic clone (different Y-deletion clone, or a complementary Y-knock-back) and look for TFs that reproduce — most published applications use multiple models or cohort-scale data precisely to filter noise.

**Benchmarks/thresholds that would change the recommendation:**
- If sample size could be expanded to n ≥ 6/group with explicit batch design, the noise concerns shrink substantially and ULM becomes more defensible.
- If matched ATAC-seq is generated, decoupleR-ULM becomes a corroborative layer rather than primary evidence, and the standard becomes much more defensible.
- If a TF survives the consensus call, the cell-cycle/IFN/senescence filter, and one orthogonal validation (motif, ChIP-qPCR, knockdown), it is publishable as a load-bearing finding.

---

## Caveats

1. **The CollecTRI + ULM application literature on bulk RNA-seq is still young** (CollecTRI was published November 2023). Several of the application papers cited above (Baudre 2025, Afenteva 2025, the Cambridge PDAC study) should have their final publication details verified directly before manuscript citation.
2. **There is no published benchmark of ULM-CollecTRI specifically on 3 vs 3 bulk isogenic perturbations.** The closest is Trescher & Leser 2019, which used somewhat older methods/networks and found serious failure rates (15/54 TFs in top 5%); the saezlab benchmarks are on larger cohorts or simulations. There is therefore a real evidence gap on the user's exact use case.
3. **The KnockTF ground truth used to benchmark CollecTRI overlaps with the literature curations that built CollecTRI** — the TFActProfiler 2025 preprint explicitly excluded KnockTF for this reason, calling it data leakage. CollecTRI's reported AUROC is therefore likely optimistic relative to truly held-out perturbations.
4. **No method here measures TF protein activity.** All produce mRNA footprints. For a mechanistic claim about a TF in LOY, biochemical validation (Western for nuclear protein, ChIP-qPCR, or knockdown/rescue) is required.
5. **The user's design (isogenic 3v3) is *favorable* in the sense that genetic background is controlled, but *unfavorable* because the perturbation (LOY) is pleiotropic** — it removes the male-specific Y-linked coding genes (the literature variously gives "more than 70 protein-coding genes," "around 100," or "70 to 200" depending on the inclusion of non-coding elements; Forsberg et al. 2020, *J Intern Med*, PMC6996334 state "the human Y chromosome harbors more than 70 protein-coding genes, none of them are essential for the survival of somatic cells"), including epigenetic regulators (KDM6A/UTY, KDM5D), and creates secondary genomic-instability signatures. This is the worst case for footprint-based TF inference: many simultaneous primary perturbations producing convergent downstream signatures that mimic specific TF activity changes.
6. **Bottom line for a committee meeting or manuscript review:** present ULM-CollecTRI results as "transcriptional signatures consistent with altered activity of TF X" — never as "TF X is activated/inactivated" — unless paired with at least one orthogonal piece of evidence (motif accessibility, knockdown, ChIP, or independent isogenic replication). This phrasing is both honest about the method's limits and defensible to a skeptical reviewer who has read Trescher & Leser 2019 or the EPEE paper.