# A549 primary RINe-adjusted nuclear-gene analysis

Author-selected on 2 October 2026. The primary model is unpaired limma-voom with condition and continuous RINe, robust empirical Bayes moderation, four sgCTRL and four Y-depleted libraries, and Benjamini–Hochberg correction across all 22,843 expressed nuclear genes. TMM effective library sizes remain those estimated from the 22,792 autosomal/non-PAR X genes. Positive logFC denotes higher expression in Y-depleted cells.

| File | Contents |
|---|---|
| [Full primary results](DE_nuclear_RINe_adjusted.csv) | All 22,843 genes; 27 with FDR <0.05 |
| [Significant genes](significant_genes.csv) | 17 non-PAR Y genes and ten non-Y/non-PAR genes (eight higher, two lower) |
| [Condition-only sensitivity](DE_nuclear_unadjusted.csv) | Same 22,843 genes; 16 significant, all non-PAR Y |
| [PAR follow-up](PAR_RINe_adjusted.csv) | Full-family FDR and exploratory BH_PAR_family across 31 PAR features |
| [Sample metadata](sample_metadata.csv) | Same eight samples and RINe values |
| [Design matrix](design_RINe.csv) | Condition plus RINe |

The model and gene-family choice supersede the earlier primary 22,792-gene condition-only reporting. Earlier fits remain in the [historical sensitivity tables](../a549_RINe_sensitivity/README.md). No statistical model was rerun for this adoption; saved model coefficients, P values and FDR were checked against the fitted R objects. Shared-gene effects differ negligibly from the earlier RINe-adjusted fit; adding strong Y-linked P values changes BH ranks and accounts for the ten shared-gene discoveries. See the [source audit](../../../figures/chapter3/source_data/a549_primary_nuclear_RINe/README.md).

PCA/QC still use the 22,792 autosomal/non-PAR X genes. Regional CPM percentages and X-paralog descriptive effects are unchanged. Existing pathway rankings and CAMERA tests are not RINe-adjusted; CAMERA retains the condition-only design and default fixed inter-gene correlation of 0.01. Cross-line descriptive comparisons are unchanged.
