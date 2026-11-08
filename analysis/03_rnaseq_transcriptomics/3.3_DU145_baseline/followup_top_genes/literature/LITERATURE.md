# Literature check: DU145 LOY RNA-seq top genes (GRHL2 axis, Y genes, interferon, NE hypothesis)

Compiled 2026-09-23. Every citation below was retrieved from PubMed/PMC or the publisher page during
this search; verification status is stated per item. "Verified" = I read the abstract (and, where
noted, relevant full-text passages). "Unverified" = found only as a secondary mention; not used to
support a claim.

**The result being checked.** In Y-depleted DU145 (partial Y loss: RPS4Y1/ZFY retained at ~50%;
USP9Y, DDX3Y, UTY, KDM5D, EIF1AY ~95% depleted), the largest non-Y change is **GRHL2 down ~15x**,
with coordinate loss of epithelial-identity/junction genes (TACSTD2, CLDN7, AP1M2, EPS8L1, SCEL,
MAL2, TMC4, ST14, OVOL1, S100A14), loss of interferon/inflammatory genes (BST2, CXCL1, LTB; Hallmark
IFN-γ / IFN-α / inflammatory response negatively enriched), TF-target inference showing GRHL2 targets
down (top hit), STAT1 and CEBPD down, **ZEB1 targets up**, and up-genes skewed to matrix/mesenchymal
(COL4A2, CREB3L1, F13A1, RSPO3, ADAM12).

---

## Q1. GRHL2–ZEB1 reciprocal loop, and which DU145 down-genes are established direct GRHL2 targets

### 1.1 The double-negative loop — well established, but in breast/ovarian, not prostate

| Finding | Citation | Status |
|---|---|---|
| GRHL2 is directly repressed by ZEB1, and ZEB1 is a direct target of repression by GRHL2 → the two form a double-negative feedback loop; GRHL2 knockdown in breast lines induces EMT-like morphology; GRHL2 protein is lost at the invasive front of primary breast tumours | Werner S et al. *J Biol Chem* 2013;288(32):22993–23008. "Dual roles of the transcription factor grainyhead-like 2 (GRHL2) in breast cancer." DOI 10.1074/jbc.M113.456293; PMID 23814079 | Verified (abstract) |
| Reciprocal GRHL2/ZEB1 loop controls epithelial vs mesenchymal phenotype; GRHL2 blocks Six1-driven ZEB1 transactivation; TGF-β + Wnt repress GRHL2 via direct ZEB1 binding at the GRHL2 promoter; GRHL2 is down-regulated in recurrent, EMT-like tumours | Cieply B, Farris J, Denvir J, Ford HL, Frisch SM. *Cancer Res* 2013;73(20):6299–309. DOI 10.1158/0008-5472.CAN-12-4082; PMID 23943797 | Verified (abstract) |
| GRHL2 as "gatekeeper" of the epithelial state in ovarian cancer via direct transcriptional activation of miR-200b/a → ZEB1 → CDH1; GRHL2 knockdown → EMT, migration, invasion; loss of GRHL2 raises H3K27me3 at GRHL2 sites | Chung VY et al. *Sci Rep* 2016;6:19943. DOI 10.1038/srep19943; PMID 26887977 | Verified (abstract + full text) |
| GRHL2 knockdown causes CpG hypermethylation and loss of H3K4me3/H3K27ac at epithelial gene promoters, mirroring progressive EMT states — i.e. GRHL2 loss is an *epigenetic* collapse of the epithelial program, not just a TF stoichiometry change | Chung VY et al. *Commun Biol* 2019;2:272. DOI 10.1038/s42003-019-0506-3; PMID 31372511 | Verified (abstract + full text) |
| **Important caveat**: in luminal breast cancer, GRHL2 ChIP-seq + nascent (Bru-seq) transcription after conditional GRHL2 knockout confirm CLDN4 and OVOL2 as direct targets, but explicitly **rule out CDH1 and ZEB1 as direct GRHL2 targets in that context**. The loop is therefore context-dependent, and "GRHL2 directly represses ZEB1" should not be stated as universal | Wang Z et al. *Cell Commun Signal* 2023;21(1):15. DOI 10.1186/s12964-022-01029-5; PMID 36691073 | Verified (abstract + full text) |

Chung 2016 likewise reports **no GRHL2 binding at the ZEB1 promoter** in ovarian cells (only a weak
intron-1 peak), with ChIP-qPCR failing to confirm ZEB1 as a direct target. So: the loop is real at the
phenotypic level and at the ZEB1→GRHL2 direction (Werner 2013, Cieply 2013), but the
GRHL2→ZEB1 direction is at least partly indirect (via miR-200b/a, Chung 2016) and is not reproduced in
every lineage.

### 1.2 Direct-target table vs the DU145 down-gene list

"Direct" = GRHL2 ChIP-seq/ChIP-qPCR binding at the locus **plus** loss of expression on GRHL2
knockdown/knockout. "KD-regulated" = expression change on GRHL2 perturbation, binding not at a
promoter or not validated. "Co-expression only" = correlated with GRHL2 in epithelial signatures.

| DU145 down-gene | Evidence class | Key citation(s) | Status |
|---|---|---|---|
| **TACSTD2 (TROP2)** | **Direct** — GRHL2 binding at promoter + marked down-regulation after GRHL2 knockdown (RT-qPCR, OVCA429) | Chung 2016 (PMID 26887977). Also GRHL2/TACSTD2 co-expression noted in Wang 2023 (PMID 36691073) | Verified (full text) |
| **ST14** (matriptase) | **KD-regulated + bound** — down after GRHL2 KD, with GRHL2 binding at intronic/distal (not promoter) sites; GRHL2 site at ST14 shifts active→poised during EMT | Chung 2016; Chung 2019 (PMID 31372511) | Verified (full text) |
| **MAL2** | **KD-regulated + bound** (intronic/distal GRHL2 site) | Chung 2016 | Verified (full text) |
| **S100A14** | **KD-regulated, promoter CpG at a GRHL2 site** — down after GRHL2 KD, promoter gains methylation at GRHL2-bound CpGs; restored by GRHL2 re-expression | Chung 2019 | Verified (full text) |
| **ESRP1** | **Direct/near-direct** — GRHL2-bound, down after KD; up on GRHL2 re-expression; also a GRHL2-dependent NNE gene in mouse neural tube | Chung 2016; Chung 2019; Ray HJ & Niswander LA, *Development* 2016;143(7):1192–204, DOI 10.1242/dev.129825, PMID 26903501 | Verified |
| **CLDN4** | **Direct** — GRHL2 binds the Cldn4 core promoter (mouse) and the CLDN4 promoter in human cells; required for Cldn4 expression; nascent-transcription drop after GRHL2 loss | Werth M et al. *Development* 2010;137(22):3835–45, DOI 10.1242/dev.055483, PMID 20978075; Senga K et al. *Mol Biol Cell* 2012;23(15):2845–55, DOI 10.1091/mbc.E12-02-0097, PMID 22696678; Wang 2023 | Verified |
| **RAB25** | **Direct** — GRHL2 target validated by ChIP-qPCR; part of the Grhl2/Ovol2→Cldn4/Rab25 network | Aue A et al. *J Am Soc Nephrol* 2015;26(11):2704–15, DOI 10.1681/ASN.2014080759, PMID 25788534; Senga 2012; Chung 2016 | Verified |
| **CDH1** | **Direct in most systems, NOT in luminal breast** — Grhl2 binds an intron-2 cis-element that loops to the promoter; but ruled out as direct in MCF7 KO/Bru-seq | Werth 2010; Chung 2016; **contra** Wang 2023 | Verified (both directions) |
| **OVOL2** | **Direct** — GRHL2 promoter binding, conserved across luminal breast lines; Grhl2→Ovol2→Cldn4/Rab25 in kidney | Aue 2015; Wang 2023 | Verified |
| **OVOL1** | Co-expression / same epithelial module; no GRHL2 ChIP evidence found | Kohn KW et al. *PLoS One* 2014;9(6):e99269, DOI 10.1371/journal.pone.0099269, PMID 24940735 | Verified (co-expression only) |
| **CLDN7** | Co-expression / epithelial module; falls with GRHL2 in the SSEA-4+ DU145 subpopulation (see Q2) | Kohn 2014; Sivasubramaniyan 2015 (PMID 25978997) | Verified (correlative) |
| **AP1M2** | Co-expression only (tight-junction/polarity module) | Kohn 2014 | Verified (correlative) |
| **TMC4** | Co-expression only | Kohn 2014 | Verified (correlative) |
| **SCEL** | **No GRHL2 evidence found** | — | Searched, nothing found |
| **EPS8L1** | **No GRHL2 evidence found.** (EPS8 — the paralog, not EPS8L1 — is GRHL2-bound in Chung 2016 but was only modestly regulated) | Chung 2016 | Searched; paralog only |

Related GRHL2-dependent epithelial genes that are *not* in the DU145 list but would be worth checking
in the data: SPINT1 (Walentin K et al. *Development* 2015;142(6):1125–36, DOI 10.1242/dev.113829,
PMID 25758223; Chung 2019), ERBB3 (Werner 2013; Chung 2016), ELF3, EPCAM, ARHGEF19, PROM2, CGN, PVRL4
(Chung 2016/2019), CLDN3 (Senga 2012), NKX2-1 (Varma S et al. *J Biol Chem* 2012;287(44):37282–95,
DOI 10.1074/jbc.M112.408401, PMID 22955271 — direct Grhl2 target, positive feedback loop).

**Bottom line for Q1:** roughly half the DU145 epithelial down-list (TACSTD2, ST14, MAL2, S100A14,
ESRP1, plus CLDN4/RAB25/OVOL2/CDH1 if present) are documented GRHL2 targets by ChIP + knockdown. The
rest (CLDN7, AP1M2, TMC4, SCEL, EPS8L1, OVOL1) belong to the same co-regulated epithelial/tight-junction
module (Kohn 2014) but have no direct GRHL2 binding evidence. That is enough to call the DU145 down-set
"a GRHL2-target-enriched epithelial program", not "all direct GRHL2 targets".

### 1.3 The NKX2-1 link (SFTPA1, NAPSA)

**Varma 2012 (PMID 22955271) is directly relevant and easy to miss.** Nkx2-1 is a *direct
transcriptional target of Grhl2*, and the two form a positive feedback loop in lung epithelium; Grhl2
knockdown lowers Nkx2-1. So SFTPA1/NAPSA falling alongside GRHL2 in DU145 is **internally coherent**
with a GRHL2→NKX2-1→surfactant/napsin cascade rather than being an unrelated oddity. Verified
(abstract). Caveat: this is lung epithelium, and the DU145 NKX2-1 expression itself is atypical.

---

## Q2. GRHL2 in prostate cancer, and in DU145 specifically

| Finding | Citation | Status |
|---|---|---|
| GRHL2 is an **AR coregulator**: colocalises with AR on chromatin, is amplified/upregulated in prostate cancer, maintains AR expression, enhances AR transcriptional activity, and is itself AR-regulated (positive feedback). Also binds and regulates AR splice variants. **These oncogenic functions are counterbalanced by GRHL2 suppressing EMT and invasion**; mechanistic evidence that AR *assists* GRHL2 in maintaining the epithelial phenotype | Paltoglou S et al. *Cancer Res* 2017;77(13):3417–3430. DOI 10.1158/0008-5472.CAN-16-1616; PMID 28473532 | Verified (abstract + full text scan) |
| Editorial framing of the above ("dichotomous functions") | Thoma C. *Nat Rev Urol* 2017;14(7):391. DOI 10.1038/nrurol.2017.78; PMID 28534537 | Verified (title/citation) |
| CAF-conditioned medium induces stable EMT in **PC3 and DU145 but not LNCaP/C4-2B**; DNMT3A-dependent promoter methylation **silences CDH1 and GRHL2**, explicitly naming GRHL2 as "the direct repressor of ZEB1"; DNMT3A knockdown prevents EMT entry. TCGA-PRAD clustering on this signature defines a **CAF+/AR− worst-prognosis** cluster | Pistore C et al. *Oncogene* 2017;36(40):5551–5566. DOI 10.1038/onc.2017.159; PMID 28581528 | Verified (abstract) |
| **SSEA-4+ subpopulations of DU145** spontaneously lose the epithelial phenotype: Claudin-7, E-cadherin, **ESRP1 and GRHL2 are all down-regulated** in the SSEA-4+ fraction of DU145 (and HCT-116), which forms fibroblast-like colonies, is enriched for pluripotency markers, gains EMT/mesenchymal markers, and is more tumorigenic and invasive | Sivasubramaniyan K et al. *Glycobiology* 2015;25(8):902–17. DOI 10.1093/glycob/cwv032; PMID 25978997 | Verified (abstract) |
| GRHL2 in whole blood is a prognostic biomarker in mHSPC/mCRPC (correlative, CTC-based) | Kwan EM et al. *Transl Androl Urol* 2021;10(4). DOI 10.21037/tau-20-1444; PMID 33968657 | Verified (citation/title only) |
| **DU145 in the CRPC taxonomy**: DU145 and PC3 define the AR-negative **"stem-cell-like" (CRPC-SCL)** subtype, driven by AP-1 (FOSL1) with YAP/TAZ–TEAD; CRPC-SCL is the second most common CRPC subtype (~28% of patients). DU145 was one of the two experimental models used | Tang F et al. *Science* 2022;376(6596):eabe1505. DOI 10.1126/science.abe1505; PMID 35617398 | Verified (abstract + full text scan) |
| ZEB1 silencing in DU145 reduces motility, invasion and colony formation; ZEB1 overexpression drives EMT (VIM up, CDH1 down) in prostate lines | Orellana-Serradell O et al. *Asian J Androl* 2018;20(3):294–299. DOI 10.4103/aja.aja_61_17; PMID 29271397 | Verified (abstract) |
| ZEB1 is a driver of EMT and docetaxel resistance in DU145 and PC-3 docetaxel-resistant sublines, via repression of E-cadherin; ZEB1 rises in prostate tumours after docetaxel | Hanrahan K et al. *Mol Oncol* 2017;11(3):251–265. DOI 10.1002/1878-0261.12030; PMID 28133913 | Verified (abstract) |

**No study was found that knocked down or overexpressed GRHL2 in DU145 or PC3 directly.** I searched
PubMed for GRHL2/grainyhead × DU145/PC3/PC-3 and found only Pistore 2017 (GRHL2 silenced indirectly by
CAF-driven methylation) and Sivasubramaniyan 2015 (GRHL2 down in the SSEA-4+ DU145 fraction). This is a
genuine gap — state it as such rather than implying a DU145 GRHL2 loss-of-function experiment exists.

**Why the AR-negativity of DU145 matters.** Paltoglou's oncogenic arm of GRHL2 runs through AR. DU145
is AR-negative, so the AR-dependent, growth-promoting half of GRHL2 biology does not apply; what remains
is the AR-independent, EMT-suppressing half. That is a clean argument for why GRHL2 loss in DU145 should
read as *epithelial de-differentiation* rather than as loss of an AR-axis oncogene — and it is a
prediction the thesis can state.

---

## Q3. Y genes → GRHL2 / epithelial state / EMT

### 3.1 Direct Y-gene → epithelial/EMT evidence

| Finding | Direction | Citation | Status |
|---|---|---|---|
| **KDM5D (JARID1D) is an anti-invasion factor in prostate cancer.** Knockdown in **DU145** (the exact line) increases invasion in vitro and in vivo; KDM5D represses MMP1/2/3/7 and **Slug (SNAI2)** by demethylating H3K4me3 at their promoters; KDM5D is down-regulated and frequently deleted in metastatic prostate tumours | KDM5D loss → *more* invasive, Slug up | Li N et al. *Cancer Res* 2016;76(4):831–43. DOI 10.1158/0008-5472.CAN-15-0906; PMID 26747897 | Verified (abstract + full-text methods confirming DU145 shRNA, H3K4me3 ChIP-seq in DU145) |
| KDM5D loss confers docetaxel resistance in prostate cancer via AR; KDM5D interacts with AR and demethylates H3K4me3; KDM5D low in CRPC and associated with worse prognosis | KDM5D loss → AR dysregulation, docetaxel resistance | Komura K et al. *PNAS* 2016;113(22):6259–64. DOI 10.1073/pnas.1600420113; PMID 27185910 | Verified (abstract). Note: mechanism is AR-dependent (LNCaP), so it does **not** transfer to AR-negative DU145 |
| **KDM5D upregulation represses regulators of the epithelial tight junction** and MHC-I components in KRAS-mutant colon cancer; *Kdm5d* deletion **increases tight junction integrity** and decreases invasiveness; forced Kdm5d increases invasion | KDM5D **high** → tight junctions down. **Opposite direction to DU145** | Li J et al. *Nature* 2023;619(7970):632–639. DOI 10.1038/s41586-023-06254-7; PMID 37344599 | Verified (abstract) |
| KDM5D suppresses EMT in gastric cancer by demethylating the CUL4A promoter; KDM5D knockdown raises **ZEB1** | KDM5D loss → ZEB1 up | Shen X et al. *J Cell Biochem* 2019;120(8):12247–12258. DOI 10.1002/jcb.27308; PMID 30864186 | Verified (abstract). Indirect (gastric) |
| **LOY drives EMT and lineage plasticity in lung adenocarcinoma** using isogenic A549 single-cell clones ± Y; mechanism = haploinsufficiency of dosage-sensitive regulators → widespread DNA **hypo**methylation at EMT gene promoters (THY1, LOX); LOY increases epigenetic heterogeneity, destabilises chromatin, and facilitates engraftment/metastasis | LOY → EMT/plasticity | Schlüter K, Chen M-JM, Altun G, … Plass C, Llamazares-Prada M. bioRxiv, posted 5 June 2026. "Loss of the Y chromosome drives epigenetic and transcriptomic plasticity in lung adenocarcinoma." DOI 10.64898/2026.06.02.729627 | Verified (title, author list, date, DOI, abstract retrieved from bioRxiv). **Preprint, not peer-reviewed.** Abstract names no Y gene, no GRHL2, no ZEB1, no interferon |
| LOY bladder cancer cells: DDR2 deletion in LOY cells reduces glycolysis, **reduces EMT and stemness features** — i.e. LOY cells carry EMT/stemness features to begin with | LOY → EMT/stemness (correlative) | Subramaniyan B et al. bioRxiv 2025.08.26.672471, posted 31 Aug 2025. DOI 10.1101/2025.08.26.672471; PMID 40909597 | Verified (abstract). Preprint |
| ZEB1 knockdown in male GBM lines **reduces KDM5D and UTY** expression; ZEB1 prognostic effect is sex-dependent | ZEB1 → Y genes (reverse direction!) | Whittaker BE et al. *Cells* 2026;15(15):1428. DOI 10.3390/cells15151428; PMID 42587836 | Verified (abstract). Correlative, different tissue |

### 3.2 UTY / KDM6A and epithelial identity

| Finding | Citation | Status |
|---|---|---|
| UTY(KDM6C) is a genuine, if weak, H3K27 demethylase (not catalytically dead) | Walport LJ et al. *J Biol Chem* 2014;289(26):18302–13. DOI 10.1074/jbc.M114.555052; PMID 24798337 | Verified |
| UTX/UTY suppress leukemogenesis through **noncatalytic** functions; loss changes H3K27ac, chromatin accessibility and (predominantly) **H3K4me1**, alters ETS/GATA factor binding, and acts via **coordination of the COMPASS complex** | Gozdecka M et al. *Nat Genet* 2018;50(6):883–894. DOI 10.1038/s41588-018-0114-z; PMID 29736013 | Verified |
| KDM6A loss induces **squamous-like, metastatic** pancreatic cancer selectively in females via COMPASS deregulation and super-enhancer activation of ΔNp63/MYC/RUNX3; in males the same subtype required **concomitant UTY and KDM6A loss** → overlapping, largely demethylase-independent tumour-suppressor roles. This is the cleanest published demonstration that **UTY buffers KDM6A for lineage identity** | Andricovich J et al. *Cancer Cell* 2018;33(3):512–526.e8. DOI 10.1016/j.ccell.2018.02.003; PMID 29533787 | Verified |
| UTX(KDM6A) promotes epithelial differentiation **noncatalytically** in skin, directly regulating Cdh1, **Grhl3**, Ctnnb1 and retinoid genes — and the phenotype appears **only in homozygous females, because UTY compensates in males** | Pacella GN et al. *PNAS* 2025;122(20):e2422971122. DOI 10.1073/pnas.2422971122; PMID 40372430 | Verified. **Closest published link from a Y gene to a Grainyhead-family factor and epithelial identity** — but it is *Grhl3*, not GRHL2, and it is mouse skin |
| KDM6A suppresses proliferation and EMT in prostate cancer cells; KDM6A is down in CRPC vs primary; depletion activates EMT via HMGB2 | Ni Y et al. *Urol Oncol* 2026;44(10):652–661. DOI 10.1016/j.urolonc.2026.07.029; PMID 42636791 | Verified. Prostate, but X-linked paralog |
| **Uty** disruption in haematopoietic cells phenocopies LOY's profibrotic heart-failure attributes; Uty loss causes epigenetic alterations and TGF-β-dependent pathology | Horitani K et al. *Nat Cardiovasc Res* 2024;3(3):343–355. DOI 10.1038/s44161-024-00441-z; PMID 39183958 | Verified. Already in the project's context notes |
| **GRHL2 physically interacts with KMT2C and KMT2D (MLL3/MLL4)**; these interactions drive MET, ICAM-1 expression and NK sensitivity | MacFawn I et al. *Mol Immunol* 2019;105:137–149. DOI 10.1016/j.molimm.2018.11.006; PMID 30508726 | Verified |
| KMT2C/D are required to **amplify** GRHL2-driven enhancer activation: GRHL2 binds its sites independently of KMT2C/D, but without them H3K4me1/2, P300 recruitment, H3K27ac and transcriptional output collapse | Boileau RM, Chen KX, Blelloch R. bioRxiv 2026.03.10.710694, posted 16 Mar 2026. DOI 10.64898/2026.03.10.710694; PMID 41890106 | Verified (abstract). Preprint, mouse ESCs |
| An MLL3/GRHL2 complex regulates the squamous-cancer epigenome and an MLL3/GRHL2–**IRF1** axis promoting Th1 chemokines and anti-tumour immunity | Nam C et al. *J Exp Med* 2025;222(4):e20240758. DOI 10.1084/jem.20240758; PMID 39964485 | Verified |

**The mechanistic chain that the literature *permits* but does not demonstrate** (flag as a hypothesis,
not a finding): UTY/KDM6A act largely non-catalytically through COMPASS/KMT2C-D coordination
(Gozdecka 2018; Andricovich 2018); GRHL2's enhancer output depends on KMT2C/D (MacFawn 2019;
Boileau 2026); UTX loss in mouse skin dysregulates Grhl3 and epithelial differentiation and is
normally buffered by UTY (Pacella 2025). **No paper tests UTY → GRHL2.** Every link in the chain is
published; the chain itself is not.

### 3.3 What is *not* supported

- **No published evidence that any Y-linked gene regulates GRHL2.** I searched GRHL2 × (KDM5D, UTY,
  KDM6C, DDX3Y, USP9Y, EIF1AY, Y chromosome) and found nothing.
- **Schlüter et al. do not report GRHL2, ZEB1 or interferon in the abstract**, and name THY1 and LOX as
  the mechanistic examples. Their methylation direction is **hypo**methylation at EMT promoters (i.e.
  mesenchymal genes de-repressed), which is mechanistically *different* from — though phenotypically
  compatible with — the Chung 2019 model in which epithelial promoters gain methylation when GRHL2 is
  lost. Do not conflate the two.
- Schlüter's model is **selected single-cell clones**, matching the project's existing
  acute-vs-selected distinction; DU145 Y-depleted populations are closer to acute/pooled.

---

## Q4. LOY and tumour-intrinsic interferon signalling

**This is where the DU145 direction needs the most care: the published LOY literature points the
opposite way, or is silent on the cell-intrinsic question.**

| Study | System | What was reported | Direction relative to "LOY → IFN down (cell-intrinsic)" |
|---|---|---|---|
| Abdel-Hafiz HA et al. *Nature* 2023;619(7970):624–631. "Y chromosome loss in cancer drives growth by evasion of adaptive immunity." DOI 10.1038/s41586-023-06234-x; PMID 37344596 (**Verified**, abstract + full text) | MB49 bladder, Y+/Y− (pooled clones and CRISPR whole-Y deletion), in vivo | Y− tumours are more aggressive **only in immunocompetent hosts**, T-cell dependent; CD8 exhaustion; better anti-PD-1 response. Bulk RNA-seq of **whole tumours grown in wild-type mice** gave 958 DEGs; GSEA showed **Y+ tumours enriched for immune responses including "response to IFNγ" (NES 1.87)** — i.e. **Y− tumours lower on the IFN-γ response axis**. Uty and Kdm5d each partly rescue. Y+/Y− grow identically in vitro | **Same direction as DU145** for IFN-γ — but that signature came from a whole tumour containing immune cells, so it is *not* a clean cell-intrinsic claim |
| Li J et al. *Nature* 2023;619(7970):632–639 (**Verified**) | iKAP mouse colon, KRAS-mutant | KRAS*→STAT4 drives **KDM5D up**; high KDM5D **represses MHC-I components** and tight-junction regulators; *Kdm5d* deletion improves CD8 killing | **Opposite gene direction**: here it is Y-gene *gain*, not loss, that suppresses antigen presentation. Cited correctly, this is the "KDM5D-directionality paradox", not support |
| Chen X et al. *Nature* 2025;642(8069):1041–1050. "Concurrent loss of the Y chromosome in cancer and T cells impacts outcome." DOI 10.1038/s41586-025-09071-2; PMID 40468066 (**Verified**) | Pan-cancer, 29 tumour types, scRNA-seq + mouse | LOY in CD4+/CD8+ T cells carries **immunosuppression transcriptomic signatures**; epithelial LOY predicts LOY in benign cells; concurrent epithelial + T-cell LOY = worst survival | Supports LOY↔immune suppression, but the transcriptomic effect is described in **T cells**, not in tumour cells |
| MacFawn I et al. *Mol Immunol* 2024;170:156–169. "Grainyhead-like-2, an epithelial master programmer, promotes interferon induction and suppresses breast cancer recurrence." DOI 10.1016/j.molimm.2024.04.012; PMID 38692097 (**Verified**) | Breast cancer cells | **GRHL2 promotes type-I and type-III IFN responses to dsRNA**, enhancing IRF3 and relA/NF-κB activation and **IRF1** expression; a functional GRHL2 binding site sits in the **IFNL1 promoter**; multiple interferon response genes **anti-correlate with EMT** | **This is the single most useful citation for the DU145 result.** It supplies a *cell-intrinsic, GRHL2-mediated* mechanism by which loss of the epithelial program lowers ISG expression — no immune cells required |
| MacFawn 2019 (PMID 30508726, **Verified**) | Breast | GRHL2 raises ICAM-1 and NK sensitivity via KMT2C/D; IFN response genes anti-correlate with EMT | Same direction |
| Nam C et al. *J Exp Med* 2025 (PMID 39964485, **Verified**) | Squamous (UASCC) | MLL3/GRHL2–**IRF1** axis drives Th1 chemokines and T-cell infiltration | Same direction; consistent with STAT1/IRF1-target loss |

**Conclusion for Q4.** There is **no published paper showing that Y-deficient tumour cells have reduced
tumour-intrinsic ISG expression in a pure cell-autonomous setting**. Abdel-Hafiz 2023 is the closest and
is directionally concordant, but its IFN-γ GSEA is from bulk tumours in immunocompetent mice. The
mechanistic bridge that *does* exist is **GRHL2 → IFN** (MacFawn 2024): if LOY collapses GRHL2, reduced
ISG/STAT1-target expression is the predicted downstream consequence. That makes the DU145 IFN result a
**coherent consequence of the GRHL2 result rather than an independent second finding** — which is a
stronger and safer way to write it.

---

## Q5. Neuroendocrine hypothesis

| Finding | Citation | Status |
|---|---|---|
| **Direct LOY→NE claim**: an AUA 2025 abstract (IP14-09) titled "Y chromosome loss promotes neuroendocrine differentiation in bladder cancer", authors listed as Yiling Shen, Xingyu Chen, Suhyeon Choi, Simon R.V. Knott, Dan Theodorescu (Los Angeles, CA), *J Urol* 2025 AUA supplement | I could **not** retrieve the abstract text: the AUA PDF returned HTTP 403 and it is not indexed in PubMed. The title, author list and session code come from a search-result snippet only | **UNVERIFIED — do not cite as read.** Cite as a conference abstract with the caveat, or obtain the PDF via institutional access before use |
| **ZEB1-marked intermediate state is the cell of origin for NEPC**: prostate cells pass through a Zeb1+ state with EMT, stemness and NE features; lineage tracing shows Zeb1+ epithelioid cells are putative NEPC cells of origin; Zeb1 drives glycolysis → histone lactylation → chromatin opening → neuro-gene induction | Wang D et al. *Cell Death Differ* 2024;31(6):779–791. DOI 10.1038/s41418-024-01295-5; PMID 38654072 | Verified (abstract). **Strong indirect support**: ZEB1-target-up is exactly the state this paper places upstream of NE transdifferentiation |
| Zeb1 and the SK3 channel are up in NE-CRPC and in enzalutamide-neurodifferentiated LNCaP; both are drivers of NE differentiation | Bery F et al. *Cancers* 2021;13(12):2947. DOI 10.3390/cancers13122947; PMID 34204608 | Verified (abstract) |
| Review: EMT and stemness are widely activated in, and critical to, NE differentiation of prostate cancer | Li YF et al. *Oncol Rep* 2025;54(3). DOI 10.3892/or.2025.8942; PMID 40641143 | Verified (citation/title). Review |
| **DU145 is not a pre-NE line**: chromatin-based CRPC taxonomy places DU145 and PC3 in the AR-negative **stem-cell-like (CRPC-SCL)** subtype driven by AP-1/FOSL1 + YAP/TAZ–TEAD, which is **distinct from the CRPC-NE subtype** (SYP/CHGA/DLL3-high). CRPC-SCL is high in CD44 and **TACSTD2** | Tang F et al. *Science* 2022 (PMID 35617398) | Verified (abstract + full text) |

**No paper was found linking GRHL2 loss to neuroendocrine transdifferentiation in prostate cancer.**
Searches for GRHL2 × (neuroendocrine, NEPC, lineage plasticity, castration) returned nothing on point.

**The honest read for the thesis:** the DU145 data support an **epithelial→mesenchymal/plastic** shift
(GRHL2 down, ZEB1 targets up), and the ZEB1+ intermediate state is documented as *upstream of, and
permissive for,* NE transdifferentiation (Wang 2024; Bery 2021). But DU145 sits in the AP-1/YAP-driven
SCL subtype, not the NE subtype (Tang 2022), and nothing in the reported gene list is an NE marker. So:
**the result is consistent with a step toward lineage plasticity, and is not evidence of an NE program.**
A direct test — scoring the DEGs against an NE/neuronal signature and against Tang's SCL signature — is
the right follow-up, and the SCL comparison is probably the more informative one for DU145.

---

## What the literature supports for the thesis

1. **GRHL2 is a bona fide master regulator of the epithelial state, and the DU145 down-list is a
   GRHL2-target-enriched epithelial program.** TACSTD2, ST14, MAL2, S100A14 and ESRP1 are documented
   GRHL2 targets by ChIP + knockdown (Chung 2016 PMID 26887977; Chung 2019 PMID 31372511); CLDN4,
   RAB25, OVOL2 and (context-dependently) CDH1 are canonical GRHL2 targets (Werth 2010 PMID 20978075;
   Senga 2012 PMID 22696678; Aue 2015 PMID 25788534; Wang 2023 PMID 36691073). A 15-fold GRHL2 drop
   with this exact downstream set is internally coherent and is the expected signature.
2. **GRHL2-down / ZEB1-targets-up is the textbook configuration of the double-negative feedback loop**
   (Werner 2013 PMID 23814079; Cieply 2013 PMID 23943797; Chung 2016). The direction is right.
3. **The interferon arm is explained by the GRHL2 arm, cell-autonomously.** GRHL2 directly promotes
   type-I/III IFN induction and IRF1 expression, with a GRHL2 site in the IFNL1 promoter, and IFN
   response genes anti-correlate with EMT (MacFawn 2024 PMID 38692097; MacFawn 2019 PMID 30508726;
   Nam 2025 PMID 39964485). STAT1-target loss and the negative Hallmark IFN enrichments follow from
   epithelial-program collapse. This is the thesis's best mechanistic sentence.
4. **GRHL2 loss in an AR-negative prostate line is the "unopposed EMT-suppressor loss" scenario.**
   Paltoglou 2017 (PMID 28473532) splits GRHL2 into an AR-dependent oncogenic arm and an AR-independent
   EMT-suppressing arm; DU145 has no AR, so only the second arm is in play.
5. **DU145 is already known to lose GRHL2 when it de-differentiates.** Sivasubramaniyan 2015
   (PMID 25978997) shows GRHL2, CLDN7, ESRP1 and E-cadherin all fall together in the SSEA-4+
   DU145 subpopulation; Pistore 2017 (PMID 28581528) shows CAF-driven DNMT3A methylation silences GRHL2
   and CDH1 specifically in DU145/PC3 and not LNCaP. So the axis is *operative and dynamic in this very
   line* — the LOY result is landing on a known, mobile switch in DU145, not an inert one.
6. **A Y-gene→epithelial-identity link exists in principle.** KDM5D knockdown in DU145 itself increases
   invasion and de-represses Slug (Li N 2016 PMID 26747897); KDM6A/UTY redundancy governs squamous
   lineage identity (Andricovich 2018 PMID 29533787); UTX regulates Grhl3/Cdh1 in skin with UTY
   compensating in males (Pacella 2025 PMID 40372430); UTY/UTX act non-catalytically through COMPASS
   (Gozdecka 2018 PMID 29736013), and GRHL2's enhancer output depends on KMT2C/D (MacFawn 2019;
   Boileau 2026 PMID 41890106).
7. **LOY→EMT/plasticity is now reported independently in another line.** Schlüter et al. (bioRxiv 2026,
   DOI 10.64898/2026.06.02.729627) report LOY→EMT and lineage plasticity in isogenic A549 clones. Cite
   as a preprint.
8. **The ZEB1-up state is documented as the gateway to NEPC** (Wang 2024 PMID 38654072; Bery 2021
   PMID 34204608), so the plasticity framing is defensible even though the NE framing is not.

## What the literature does NOT support

1. **No evidence that any Y gene regulates GRHL2.** The UTY→COMPASS→GRHL2-enhancer chain is a
   hypothesis assembled from three separate literatures. It must be written as a proposed mechanism.
2. **No GRHL2 knockdown or overexpression experiment in DU145 or PC3 exists.** Do not imply one does.
3. **CLDN7, AP1M2, TMC4, SCEL, EPS8L1 and OVOL1 are not established direct GRHL2 targets.** They are
   co-regulated members of the same epithelial module (Kohn 2014 PMID 24940735). Write "GRHL2-associated
   epithelial program", not "GRHL2 targets".
4. **"GRHL2 directly represses ZEB1" is not universally true.** Wang 2023 (PMID 36691073) and Chung 2016
   both fail to find GRHL2 binding at the ZEB1 promoter; Chung's route is via miR-200b/a. Similarly CDH1
   is a direct target in some systems and explicitly not in luminal breast. State the loop at the
   phenotypic level and be specific about which direction has ChIP support.
5. **No published cell-intrinsic LOY→ISG-down result.** Abdel-Hafiz 2023's IFN-γ GSEA (Y+ > Y−) comes
   from bulk tumours grown in immunocompetent mice, and the same paper reports Y+/Y− grow identically in
   vitro. Li J 2023 runs the *opposite* gene direction (KDM5D **up** represses MHC-I). Present DU145's
   ISG loss as a **new, cell-autonomous observation** whose mechanism is supplied by GRHL2, and name the
   KDM5D-directionality paradox rather than glossing it.
6. **No support for a neuroendocrine interpretation of DU145.** No GRHL2–NE link in prostate exists in
   the literature; DU145 belongs to the AP-1/YAP-driven CRPC-SCL subtype, not CRPC-NE (Tang 2022
   PMID 35617398); and the Shen AUA 2025 abstract could not be retrieved and is **unverified**. The NE
   hypothesis should be retired for DU145 in favour of "epithelial de-differentiation / lineage
   plasticity", with the ZEB1-state→NEPC literature cited only as what plasticity could *later* permit.
7. **Directional mismatch with Schlüter worth noting rather than hiding:** their mechanism is promoter
   **hypo**methylation de-repressing mesenchymal genes (THY1, LOX), whereas the GRHL2-loss literature
   describes epithelial promoters **gaining** methylation (Chung 2019). Same phenotype, different
   mechanistic entry point.

## Suggested analyses this reading implies

- Score the DU145 DEGs against (a) a published GRHL2-knockdown signature (Chung 2016 / Wang 2023
  supplementary target lists), (b) Tang 2022's CRPC-SCL vs CRPC-NE signatures, (c) an NE/neuronal
  signature — the SCL comparison being the most informative for DU145.
- Check whether SPINT1, ERBB3, ELF3, EPCAM, CLDN3, CGN, PROM2, PVRL4 and ARHGEF19 also fall: these are
  GRHL2 targets the current top-gene list does not mention, and their behaviour would strongly
  discriminate "GRHL2 program collapse" from "generic epithelial drift".
- Check GRHL1 and GRHL3 (the paralogs) — GRHL3 is the factor UTX regulates in skin (Pacella 2025).
- Check IFNL1, IRF1, ICAM1 and NT5E specifically, given the GRHL2→IFN mechanism (MacFawn 2024/2019;
  Coban B et al. *iScience* 2024;27:109738, DOI 10.1016/j.isci.2024.109738, PMID 38706844 —
  GRHL2 suppression of NT5E/CD73; verified by citation/title only).
- The partial-Y-loss structure of the DU145 population (distal Yp retained, ~95% loss of USP9Y, DDX3Y,
  UTY, KDM5D, EIF1AY) is itself a natural experiment separating RPS4Y1/ZFY from the male-specific-region
  genes — worth stating explicitly, since it narrows the candidate causal genes to the depleted set and
  matches the UTY/KDM5D genes implicated in Abdel-Hafiz 2023 and Horitani 2024.
