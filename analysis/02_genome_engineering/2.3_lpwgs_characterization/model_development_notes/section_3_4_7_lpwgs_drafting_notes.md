# 3.4.7 - lpWGS drafting notes

**Current structure, approved 17 September 2026:** Development and lpWGS characterisation now form Section 3.4.6, organised as HT1080 (3.4.6.1), A549 (3.4.6.2), DU145 (3.4.6.3), H2452 (3.4.6.4), HepG2 (3.4.6.5) and efficiency-improvement attempts (3.4.6.6). The earlier section assignments below are historical. Lead each genomic-characterisation account with lpWGS before integrating matched ddPCR and expression. See the [adopted structure](integrated_model_development_plan.md); complete HT1080 lpWGS next.

## Author direction (17 September 2026)

**A549 drafting update:** the author subsequently requested work on A549 and confirmed that its sequenced edited population received two successive CP + Q53 nucleofections. Section 3.4.6.2 now includes verified ddPCR summaries and a source-checked lpWGS account. See [the A549 drafting notes](a549_results_drafting_notes.md) for exact values, normalization distinctions and outstanding specimen links. This confirmation does not establish recovery timing or the identities of later molecular samples.

Lead the genomic-characterisation Results with low-pass whole-genome sequencing. Start with HT1080, then A549 and DU145, before integrating locus-specific ddPCR and expression. Keep the short pharmacological optimisation account in 3.4.6: KU-55933 (ATM inhibitor) before nucleofection and panobinostat afterwards, per the author clarification.

## Located lpWGS sources

- `A:/Dropbox/LOY_paralogs/projects/thesis/Methods/METHODS_lpWGS.md`: written methods and sample names, including HT1080-sgCTRL-D0/HT1080-LOY-D0, A549-CTRL2/A549-LOY, and historical LNCaP-sgCTRL/LNCaP-LOY-2X plus treatment samples. This is a source route, not a verified population lineage or a final Methods authority.
- `A:/Dropbox/LOY_paralogs/projects/WGS/figures/HT1080_YposYneg/`: saved figures, `figure_source_stats.tsv`, `figure_chrY_region_stats.tsv` and figure captions.
- `A:/Dropbox/LOY_paralogs/projects/WGS/figures/A549_YposYneg/`: corresponding A549 outputs.
- `A:/Dropbox/LOY_paralogs/projects/WGS/figures/lpWGS_forBEN/`: another exported bundle; compare versions before selecting a canonical output.
- `.../lpWGS_forBEN/relative_copy_number.md`: describes QDNAseq 100-kb bins, GC/mappability correction and genome-median relative depth. It discusses approximate HT1080 2-copy and A549 3-copy baselines; verify these against the retained analysis rather than generalising to every sample.

## Checks before numerical prose

1. Match each lpWGS specimen and passage to the edited/control populations and ddPCR measurements. Do not automatically alias HT1080 D0_LOY to K5/MIX-MIX.
2. Verify which figure pipeline is used. The written Methods describes mosdepth windows and autosomal-median normalisation; the exported figure notes describe QDNAseq correction and genome-median normalisation. Keep the measurements and denominators explicit rather than combining them as one scale.
3. Verify the actual reference FASTA/assembly and Y sequence used. Do not copy an imprecise CHM13/Y description from an older Methods draft.
4. Confirm coverage QC, mappability filters and regional boundaries before assigning Yp/Yq retention. Sparse or repetitive-region coverage is not by itself a structural breakpoint call.
5. Recheck the previously reported HT1080 residual Yp and 13q gain. Use numerical source estimates with the correct comparator and uncertainty appropriate to the data.
6. For DU145, verify specimen-specific reference-locus dosage and distinguish baseline, second-nucleofection and treatment samples. Public DepMap WGS is contextual evidence only.

## Draft figure sequence

Main: HT1080 edited/control genome-wide profiles followed by a regional Y view with ddPCR target positions. Extend to A549 and DU145 with consistent scales where justified. Supplement: complete per-sample profiles, QC, coverage-bin distributions, normalisation checks and expanded region summaries. Do not set a fixed figure limit.

No new lpWGS analysis or numerical finding was generated in this planning update.

## HT1080 lpWGS and expression figure integration - 17 September 2026

Three Ben_Thesis-derived figures now support 3.4.6.1: genome-wide lpWGS with chromosome 13q detail, regional Y copy-number profiles, and chromosome-ordered Y/PAR RNA expression. Source statistics, floor handling and the unchanged expression values are recorded in [the figure notes](../../figures/chapter3/ch03_main_ht1080_lpwgs_confirmation_notes.md). The author assigns HT1080-LOY-D0 to the CP + Q53 2x line. Exact D0/K5 DNA-specimen identity and WGS/RNA passage matching remain unresolved; they are not necessary for describing the separate assays.
