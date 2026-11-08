# Run status

2026-09-08: ATAC processing, downstream analysis, report, and verification complete.

Open [the standalone report](report/ATACseq_report.html). Final checks are recorded in [verification.txt](report/verification.txt).

TF computational extension complete and verified on 2026-09-08: JASPAR2024 chromVAR, mapping-audited motif groups, TF RNA comparison, requested DA effect-size plots, four HOMER comparisons, and TOBIAS ATACorrect/ScoreBigwig/BINDetect. Primary HOMER: 1 gained and 31 lost enriched motifs at FDR <0.05. Independent-model chromVAR: 7 significant individual motifs. Six models pass both tests and all six also have negative TOBIAS differential scores. They are related motif models, not six independently established TF regulators.

TOBIAS compared 879 motifs using verified pools of 81,576,998 CTRL reads and 81,569,618 LOY reads. The CTCF QC uses 20,585 shared motif sites. Signal/reference checks, strict multiple-testing checks, TF RNA source comparisons, and final report checks passed. All nine TF figures were reviewed and their current PNGs are embedded in the standalone HTML report; PDF exports are retained. Full checks are in motifs/TF_extension/verification.txt and report/verification.txt. Pooled-site TOBIAS statistics are not biological-replicate tests, and footprints are occupancy predictions.

Technical-batch assignments remain unconfirmed. The primary ATAC design retains the user's independent replicates, with passage adjustment labeled as sensitivity. RNA continues to use the selected HT1080 combined-batch DESeq2 batch + treatment result. ZFY RNA is decreased but has no named motif in this JASPAR set; sex-chromosome TF coverage and unavailable RNA statistics are explicitly accounted for. See TF_EXTENSION_PLAN.md and motifs/TF_extension/METHODS.md.

- Six source BAMs audited: matching hg38 sequence dictionaries and processing history, readable indexes, valid BGZF EOF markers, SHA256 checksums recorded.
- R dependencies available. Isolated Linux environment installed at `/home/micha/.local/share/atacseq-sept08-2026-venv` in Ubuntu-24.04 WSL. MACS3 3.0.4 and pysam 0.23.3; package lock saved in config/.
- Three preparation edge-case checks passed.
- CTRL_11 pilot completed: 39,019,300 source records; 19,377,567 retained fragments; no orphaned primary records; TSS enrichment 13.6948; 85,553 MACS3 narrow peaks. Preparation elapsed 752.2 seconds.
- All six samples completed preparation and fresh MACS3 peak calling. No orphaned primary records were found.
- Consensus: 84,811 regions; 84,510 tested after chrY/PAR and abundance exclusions. Own-peak FRiP 0.473-0.506; TSS enrichment 11.69-15.30.
- Primary independent-replicate edgeR analysis: 95 gained and 663 lost regions at FDR <0.05; 32 pass the minimum-effect TREAT test at |log2FC| >0.5 and its own FDR <0.05.
- Background normalization: 65 gained, 692 lost; 656 primary calls remain significant with the same direction.
- Passage adjustment: 8,788 gained, 10,241 lost. Significance is sensitive to passage treatment despite similar effects; this is a secondary analysis.
- Annotation, promoter gene summaries, GO-BP, and motifs complete. All 758 primary significant peaks received GC/accessibility-matched controls.
- Matched-depth peak discovery complete at 13,595,988 fragments per sample. LOY_11 full-depth BAMPE and matched-depth BEDPE peak files are byte-identical.
- csaw analyzed 2,251,773 enriched windows and produced 155,533 merged regions; 4 gained and 115 lost regions have FDR <0.05. It corroborates 117 primary significant peaks in the same direction (4 gained, 113 lost).
- Original motifmatchr/JASPAR2022 screen: 17 significant lost-region motifs, no significant gained-region motifs at BH FDR <0.05 across both directions and 692 human motif matrices. The later HOMER/JASPAR2024 extension uses a different peak subset and test.
- Standalone HTML report assembled, figures reviewed, and all final output checks passed. Code/configuration hashes are recorded in logs/code_config_checksums.csv.
- Recovery: the first differential run stopped during plotting because edgeR no longer exports decideTestsDGE. Updated to decideTests and reran the full differential stage successfully; original error logs are retained.
- RNA selection confirmed by user: HT1080 combined-batch result used in Ben's thesis, `~ batch + treatment`. The source notebook identifies DESeq2 with the LOY/CTRL contrast; the thesis manuscript itself was not located.
- RNA integration completed from `HT1080_BATCHES/DE/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv`; source snapshots and SHA256 hashes are in config/RNAseq/ and integration/RNAseq/.
- Common eligible promoter-gene universe: 11,898. Joint RNA/ATAC-promoter significance: 45 genes, comprising 1 concordant gain and 44 concordant losses. Applying RNA |log2FC| >0.75 retains 43. ATAC significance remains peak-level.
- Seven ATAC gene IDs have ambiguous RNA Entrez mappings, 6,191 are absent from the selected RNA table, eight lack finite RNA statistics, and one overlaps RNA PAR coordinates; all are accounted for.
- Integration compares group effects. Direct ATAC/RNA sample matching and lpWGS dataset selection/assay matching remain pending; no copy-number adjustment has been applied.
- Promoter 0.10/0.10 gene review complete: 14 concordant increases, 76 concordant decreases, four opposing directions. Separate RNA-up/down lists, paired effect plots, chromosome summaries, and exploratory GO results are in integration/promoter_FDR010_review/. No GO term passes BH FDR 0.05 across both directions; 19 decreased genes map to chr19. Gene-set and direction checks passed.

Original source BAMs are unchanged. The full independent-replicate primary results remain the headline analysis; sensitivity and corroborated subsets are retained separately. The report includes QC, normalization and passage sensitivity, annotations, gene summaries, GO, motifs, locus views, and matched-depth peak counts.
