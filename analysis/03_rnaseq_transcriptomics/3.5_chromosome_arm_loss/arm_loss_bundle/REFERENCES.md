# Datasets and references

## Datasets analysed

All are public. Author-deposited raw count matrices were used throughout; nothing was
realigned from sequencing reads. Accessions were checked against NCBI before download.

### GSE336056 — single cell, RPTEC, CRISPR Yq and 3p
10x 3' v3.1 with cell hashing, 12 samples, 25,123 cells, one male donor line of primary
kidney tubule cells. Cas9 with a non-targeting guide, a 3p guide (two cuts flanking VHL,
PBRM1, BAP1, SETD2) and a Yq guide, collected at day 5 and day 10.
Dahiya R, van Belzen IAEM, et al. "A single chromosome 3p break initiates ccRCC evolution."
bioRxiv 2026. PMID 42465403. doi 10.64898/2026.07.10.737056

### GSE336576 — bulk, same cell system and paper as above
21 libraries: three 3p-loss clones at day 18 and day 90, three replicates each, plus three
unedited controls. Same citation as GSE336056.

### GSE222379 — bulk, ReDACT engineered arm loss
62 libraries across six isogenic backgrounds: A2058 melanoma (1q, 7p, 8q), A2780 ovarian
(1q), AGS gastric (1q), HCT116 colorectal (8q), MCF10A and MCF10A-HRAS (1q). Two to five
independent clones per arm. A2780 and AGS also include "TK" clones that went through the
same engineering without being selected for arm loss.
Girish V, Lakhani AA, et al. "Oncogene-like addiction to aneuploidy in human cancers."
Science 2023;381:eadg4521. PMID 37410869. doi 10.1126/science.adg4521

### PRJNA930361 — bulk, ARPE-19 CRISPR chrY elimination (realigned here, not a deposit)
6 libraries: parental ARPE-19 x3 vs a chrY-eliminated clone (A7123) x3. Male retinal
pigment epithelium. chrY removed with multiplexed Cas9 — guides against TSPY1 (x2),
RBMY1A1 and the DYZ3 centromere, roughly 60 simultaneous cuts, plus 16 h of the ATM
inhibitor KU-55933 to block repair.
**No GEO record and no deposited count matrix exist for this study.** The counts here were
produced by realigning SRR23563058-63 with STAR + featureCounts against the 1000 Genomes
GRCh38 analysis set (chrY pseudoautosomal regions hard-masked) and GENCODE v49 — the same
reference used for every human batch in this project. That masking is what makes the PAR
dosage readout possible: with chrY's PAR present, those reads map equally well to chrX and
chrY, become multi-mappers and are discarded.
Celli L, Gasparini P, et al. "CRISPR/Cas9-mediated loss of the Y chromosome…"
Cell & Bioscience 2024;14:19. PMID 38291538. doi 10.1186/s13578-024-01198-5
Caveats: one clone in three replicates, and the comparator is the untreated parental line —
no non-targeting or cutting-but-neutral guide, so Cas9 damage is not separated from Y
dosage.

### GSE150686 — bulk, RPE1 monosomies and inducible p53
33 libraries: monosomy 10 and 18, monosomy 13, monosomy 19p, each against its matched
diploid parent, in p53-knockout and p53-knockdown backgrounds; plus an inducible-p53 set
(plus or minus doxycycline) in diploid, monosomy 10+18 and monosomy 13 cells.
Chunduri NK, Menges P, et al. "Systems approaches identify the consequences of monosomy in
somatic human cells." Nature Communications 2021;12:5576. PMID 34552071.
doi 10.1038/s41467-021-25288-x

Note on sample naming: the GEO file named `RM_10_ip53` is annotated in the series metadata
as "Monosomy Chr. 10 and 18". The file name alone would mis-state which chromosomes are
affected; the metadata karyotype was used.

### Rejected after checking: GSE165208 (mouse)
"Mosaic loss of chromosome Y leads to clonal hematopoiesis and promotes leukemogenesis in
mice." Mouse bulk RNA-seq with sgSsty1 / sgSsty2 / sgKdm5d arms, each against its own
sgScramble — on paper the cutting-but-neutral control missing from all the human data.
Downloaded and checked; **not used**, for reasons that are worth recording:
- The Ssty arms show no chrY loss at all (chrY total 0.97 and 0.94 of scramble).
- Only 4 of 12 mouse Y genes are expressed in blood (Kdm5d, Uty, Ddx3y, Eif2s3y). Ssty1
  and Ssty2 are themselves at zero in every sample including the controls, so the guides
  cannot be verified from RNA.
- The in-vivo arm is sex-mixed (Xist 3.2-127.9 CPM within arms) and chrY rises in the
  knockout arm, which is incoherent.
- The in-vitro Kdm5d arm is a clean single-gene knockout (Kdm5d 0.68, other Y genes flat),
  not chromosome loss.
The general point: bulk RNA cannot detect Y loss unless the Y-null fraction is near-clonal.

## Reference data

- Gene models: GENCODE v49 (GRCh38).
- Chromosome arm boundaries: UCSC hg38 `cytoBand.txt.gz`, centromere ("acen") bands.
- Gene sets: MSigDB Hallmark and GO Biological Process, via msigdbr 26.1.0. The
  ribosomal-protein set is defined here as protein-coding genes whose symbol starts with
  RPL or RPS followed by a digit (97 genes), which excludes the many pseudogenes.

## Background reading

Grouped by the claim each supports. Directions that disagree between papers are marked.

**Chromosome loss triggers a damage response, not only a dosage response**
- Janssen A, et al. "Chromosome segregation errors as a cause of DNA damage and structural
  chromosome aberrations." Science 2011;333:1895. PMID 21960636. Missegregated chromosomes
  are often cut during cytokinesis, so missegregation-based systems carry DNA breaks by
  construction.
- Tsuchida CA, Brandes N, et al. "Mitigation of chromosome loss in clinical CRISPR-Cas9-
  engineered T cells." Cell 2023;186:4567. PMID 37794590. In edited primary T cells the p53
  damage and apoptosis modules are the most raised in cells that lost the targeted
  chromosome, and p53 level is inversely related to how often loss occurs.
- Haapaniemi E, et al. Nature Medicine 2018;24:927. PMID 29892067; and Ihry RJ, et al.
  Nature Medicine 2018;24:939. PMID 29892062. Cas9 cutting alone activates p53.
- Soto M, et al. Cell Reports 2017;19:2423. PMID 28636931. p53 responds in only a subset of
  cells with altered chromosome content; it tracks structural damage more than dosage.

**Monosomy specifically**
- Chunduri NK, et al. 2021 (above). Ribosome biogenesis and translation fall in every
  monosomy; the proliferation defect tracks the number of genes on the lost chromosome;
  monosomy and p53 expression are hard to combine.
- Hintzen DC, et al. PLoS One 2022;17:e0268579. PMID 35776704. Clones carrying only
  monosomies do not show the instability that trisomy-bearing clones do.

**Interferon after chromosome loss — reports disagree**
- Up, acute: Mackenzie KJ, et al. Nature 2017;548:461. PMID 28738408; Harding SM, et al.
  Nature 2017;548:466. PMID 28759889; Santaguida S, et al. Developmental Cell 2017;41:638.
  PMID 28633018.
- Up but not through cGAS-STING: Wang RW, et al. EMBO Reports 2021;22:e52032. PMID 34105235.
  Reports no significant cGAS-STING activation; NF-kB carries the phenotype.
- Down, or silenced by selection: Tripathi R, et al. Oncoimmunology 2019;8:1657374.
  PMID 31646083; Li J, et al. Nature 2023;620:1080. PMID 37612508; Hintzen DC, et al.
  EMBO Reports 2024;25:5169. PMID 39294502 — the same karyotype loses its inflammatory
  signalling after prolonged culture.
- Down by gene dosage, shown directly: Barriga FM, Tsanov KM, et al. Nature Cancer
  2022;3:1367. PMID 36344707. Deleting the interferon gene cluster along with a neighbouring
  region lowers interferon output and helps tumours evade immunity.
- Up by gene dosage, in a gain rather than a loss: Sullivan KD, et al. eLife 2016;5:e16220.
  PMID 27472900. Trisomy 21 raises interferon signalling through the extra copies of the
  four interferon receptors on chromosome 21, and lowers ribosomal protein transcripts.

**Dosage reaching RNA and protein**
- Schukken KM, Sheltzer JM. Genome Research 2022;32:1254. doi 10.1101/gr.276378.121. On
  loss, transcripts fall about 15% and proteins about 8%; the protein change is smaller
  than the gene-copy change would predict.
- Di Y, Li W, et al. Molecular Cell 2025;85:4268. PMID 41270725. Gains and losses are
  buffered by different mechanisms, so a single "amount of imbalance" axis is not right.
- Hwang S, et al. PNAS 2021;118:e2014723118. PMID 33526671. No evidence for widespread
  compensation across human autosomes.

**Ribosomal proteins and translation fall after loss — the most reproducible signal**
- Chunduri 2021 (above), every monosomy.
- Isogenic 45,X versus 46,XY cells: American Journal of Human Genetics 2025, PMID 39922196.
  Transcripts carrying a 5' terminal oligopyrimidine motif, which includes ribosomal protein
  and mitochondrial genes, are reduced.
- Sarel-Gallily R, et al. Stem Cell Reports 2025;20:102471. PMID 40185088. Across more than
  2,650 stem-cell samples, those that lost chromosome Y show most ribosomal protein genes
  reduced, and the same is seen in Turner syndrome and Diamond-Blackfan anaemia samples.
- Sullivan 2016 (above) reports the same direction in a gain, which is a caution that this
  may be a general stress response rather than specific to loss.

**Does the effect scale with size**
- For: Chunduri 2021, where the proliferation defect tracks the number of genes lost;
  Williams BR, et al. Science 2008;322:703. PMID 18974345.
- Against: Durrbaum M, et al. BMC Genomics 2014;15:139. PMID 24548329, where complex
  aneuploidy gives the same transcriptional changes as a single trisomy; Sack LM, et al.
  Cell 2018;173:499. PMID 29576454, where the relevant quantity is which genes are affected,
  not how many.
- No published work fits a quantitative law of the form "effect scales with amount lost".
  The analysis in this bundle finds no such scaling either.

**Related, not used here**
- Bosco N, et al. Cell 2023;186:1985. PMID 37075754. KaryoCreate, which moves a chosen
  chromosome into the wrong daughter cell using a catalytically dead Cas9, so no programmed
  DNA break is made. Its transcriptome deposit (GSE217325) is p53-null and has one clone per
  group, so it could not serve as a break-free comparison here.
