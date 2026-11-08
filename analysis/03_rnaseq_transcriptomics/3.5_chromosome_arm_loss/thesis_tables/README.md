# Chromosome-loss comparison: data and analysis tables

Four public datasets; two single-cell contrasts and 14 bulk contrasts. These are additional to the 18 external LOY contrasts in Section 3.3.4.4.

| Accession | Included comparison | Sampling |
|---|---|---|
| [GSE336056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336056) | Two within-guide, RNA-defined Yq/3p contrasts | Days 5 and 10; four HTO samples per guide; one RPTEC background |
| [GSE336576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336576) | Two 3p-loss-versus-unedited contrasts | Three clones at days 18 and 90; 21 libraries |
| [GSE222379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222379) | Nine ReDACT extra-arm-removal contrasts | Six backgrounds; 61 included libraries |
| [GSE150686](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150686) | Three RPE1 monosomy contrasts | p53-deficient loss/control lines; 15 included libraries |

| Bulk contrast | Case units | Control units | Case unit definition | Control unit definition |
|---|---:|---:|---|---|
| A2058_1q | 2 | 4 | Loss clones; repeat libraries summed within clone | Separate control libraries retained |
| A2780_1q_vsWT | 3 | 2 | Loss clones; repeat libraries summed within clone | Separate control libraries retained |
| A2780_1q_vsTK | 3 | 2 | Loss clones; repeat libraries summed within clone | Separate control libraries retained |
| AGS_1q_vsTK | 5 | 2 | Loss clones; repeat libraries summed within clone | Two TK control clones |
| MCF10A_1q | 4 | 2 | Loss clones; repeat libraries summed within clone | Separate control libraries retained |
| MCF10Ahr_1q | 4 | 2 | Loss clones; repeat libraries summed within clone | Separate control libraries retained |
| RPTEC_3p_early | 3 | 3 | Three clones; repeat libraries summed within clone/time | Three libraries from unedited control |
| RPTEC_3p_late | 3 | 3 | Three clones; repeat libraries summed within clone/time | Three libraries from unedited control |
| A2058_7p | 3 | 4 | Loss clones; repeat libraries summed within clone | Separate control libraries retained |
| A2058_8q | 2 | 4 | Loss clones; repeat libraries summed within clone | Separate control libraries retained |
| HCT116_8q | 2 | 2 | Loss clones; repeat libraries summed within clone | Separate control libraries retained |
| RPE1_mono1018 | 3 | 3 | Three samples from one established monosomic line | Three samples from matched p53-deficient parental line |
| RPE1_mono13 | 3 | 3 | Three samples from one established monosomic line | Three samples from matched p53-deficient parental line |
| RPE1_mono19p | 3 | 3 | Three samples from one established monosomic line | Three samples from matched p53-deficient parental line |

Shared samples: A2058 contrasts reuse WT controls; A2780 contrasts reuse loss clones; RPTEC contrasts reuse controls and sample the same clones at two times; RPE1 10+18 and 13 reuse p53-KO controls.

[Dataset provenance and Methods checks](../../../research/molecular_characterisation/arm_loss_integration_notes_2026-09-28.md).

[Sample-to-contrast membership](sample_contrast_membership.tsv) records every included library and its aggregation unit. [Design and regional RNA QC](contrast_design_qc.tsv) gives library counts, model-unit counts and regional expression summaries. [Single-cell classification and regional-expression tables](../../../figures/chapter3/source_data/arm_loss_comparison/README.md) retain the Yq/3p calls and checks formerly shown in the main figure; the standalone RNA-check supplement was retired from Results on 8 October 2026.

Complete result tables (unchanged copies from the supplied bundle):

- [A2058_1q_DE.tsv.gz](A2058_1q_DE.tsv.gz)
- [A2058_1q_GSEA.tsv.gz](A2058_1q_GSEA.tsv.gz)
- [A2058_7p_DE.tsv.gz](A2058_7p_DE.tsv.gz)
- [A2058_7p_GSEA.tsv.gz](A2058_7p_GSEA.tsv.gz)
- [A2058_8q_DE.tsv.gz](A2058_8q_DE.tsv.gz)
- [A2058_8q_GSEA.tsv.gz](A2058_8q_GSEA.tsv.gz)
- [A2780_1q_vsTK_DE.tsv.gz](A2780_1q_vsTK_DE.tsv.gz)
- [A2780_1q_vsTK_GSEA.tsv.gz](A2780_1q_vsTK_GSEA.tsv.gz)
- [A2780_1q_vsWT_DE.tsv.gz](A2780_1q_vsWT_DE.tsv.gz)
- [A2780_1q_vsWT_GSEA.tsv.gz](A2780_1q_vsWT_GSEA.tsv.gz)
- [AGS_1q_vsTK_DE.tsv.gz](AGS_1q_vsTK_DE.tsv.gz)
- [AGS_1q_vsTK_GSEA.tsv.gz](AGS_1q_vsTK_GSEA.tsv.gz)
- [GSE336056_sg3p_3pLoss_vs_intact_GSEA.tsv.gz](GSE336056_sg3p_3pLoss_vs_intact_GSEA.tsv.gz)
- [GSE336056_sg3p_3pLoss_vs_intact_MAST.tsv.gz](GSE336056_sg3p_3pLoss_vs_intact_MAST.tsv.gz)
- [GSE336056_sgYq_YqNull_vs_intact_GSEA.tsv.gz](GSE336056_sgYq_YqNull_vs_intact_GSEA.tsv.gz)
- [GSE336056_sgYq_YqNull_vs_intact_MAST.tsv.gz](GSE336056_sgYq_YqNull_vs_intact_MAST.tsv.gz)
- [HCT116_8q_DE.tsv.gz](HCT116_8q_DE.tsv.gz)
- [HCT116_8q_GSEA.tsv.gz](HCT116_8q_GSEA.tsv.gz)
- [MCF10A_1q_DE.tsv.gz](MCF10A_1q_DE.tsv.gz)
- [MCF10A_1q_GSEA.tsv.gz](MCF10A_1q_GSEA.tsv.gz)
- [MCF10Ahr_1q_DE.tsv.gz](MCF10Ahr_1q_DE.tsv.gz)
- [MCF10Ahr_1q_GSEA.tsv.gz](MCF10Ahr_1q_GSEA.tsv.gz)
- [RPE1_mono1018_DE.tsv.gz](RPE1_mono1018_DE.tsv.gz)
- [RPE1_mono1018_GSEA.tsv.gz](RPE1_mono1018_GSEA.tsv.gz)
- [RPE1_mono13_DE.tsv.gz](RPE1_mono13_DE.tsv.gz)
- [RPE1_mono13_GSEA.tsv.gz](RPE1_mono13_GSEA.tsv.gz)
- [RPE1_mono19p_DE.tsv.gz](RPE1_mono19p_DE.tsv.gz)
- [RPE1_mono19p_GSEA.tsv.gz](RPE1_mono19p_GSEA.tsv.gz)
- [RPTEC_3p_early_DE.tsv.gz](RPTEC_3p_early_DE.tsv.gz)
- [RPTEC_3p_early_GSEA.tsv.gz](RPTEC_3p_early_GSEA.tsv.gz)
- [RPTEC_3p_late_DE.tsv.gz](RPTEC_3p_late_DE.tsv.gz)
- [RPTEC_3p_late_GSEA.tsv.gz](RPTEC_3p_late_GSEA.tsv.gz)
