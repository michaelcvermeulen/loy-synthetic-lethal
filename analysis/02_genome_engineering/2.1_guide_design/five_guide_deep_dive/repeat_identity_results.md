# Repeat identities of the five selected Y-targeting guides

Suggested Results wording:

The five selected sgRNAs targeted three distinct repetitive sequence classes distributed across the Y chromosome. CP and CN targeted the centromeric DYZ3 alpha-satellite higher-order repeat array, Q53 targeted DYZ1 (HSat3) satellite repeats in Yq12, and P47 and P51 targeted TSPY-associated composite repeats on Yp. In the HG002 T2T-Y reference, these guides had 54, 53, 532, 47 and 51 exact NGG-adjacent matches, respectively, all within their corresponding repeat annotations.

| Guide | Sequence (5' to 3') | Repeat identity | Exact Y matches | Exact matches overlapping the stated annotation |
|---|---|---|---:|---:|
| CP | ACTATCTTCGCATGAGAACC | DYZ3 centromeric alpha-satellite higher-order repeats; S4CYH1L | 54 | 54/54 |
| CN | ATGATAGGTTGAACTCCCGG | DYZ3 centromeric alpha-satellite higher-order repeats; S4CYH1L | 53 | 53/53 |
| Q53 | TCGCCAGGAATTGATGTGAA | DYZ1 / HSat3, A6 subfamily, Yq12 | 532 | 532/532 |
| P47 | TTCGCGCGCAGTCCCTTAGG | TSPY-associated composite repeats, Yp | 47 | 47/47 |
| P51 | ATTAAGGCCCACCTGAGTAA | TSPY-associated composite repeats, Yp | 51 | 51/51 |

P51 includes historical aliases P53 and p51. The TSPY tandem repeat array is historically termed DYZ5, but “TSPY-associated composite repeats” directly matches the annotation used here. Membership in a TSPY composite repeat does not by itself establish overlap with a TSPY coding exon. “censat_Y” is an annotation label covering several repeat types, not a specific satellite family; avoid calling the P47/P51 targets centromeric alpha-satellite repeats.

## Evidence and interpretation

Exact 20-nt matches on both strands were independently rescanned in the local HG002 T2T-Y FASTA (CP086569.2; 62,460,029 bp) and reproduced every saved exact-match coordinate. The strand-appropriate adjacent PAM was checked at every site: all 737 matches have canonical NGG PAMs. Each 20-nt interval was intersected with the published T2T-CHM13v2.0 Cen/Sat v2.1 and composite-repeat 2022DEC tracks, converting the source 1-based positions to zero-based half-open BED coordinates. CP and CN overlap `hor_Y_1(S4CYH1L)`; Q53 overlaps 21 `hsat3_Y_*(A6/DYZ1)` intervals; all P47 and P51 sites overlap `TSPY` composite-repeat intervals.

These counts describe reference sequence matches, not measured cleavage events, repeat-unit copy numbers, or cell-line-specific target counts. This focused analysis does not use the legacy CFD scores and does not assess mismatched sites.

The site-level evidence is in [exact_Y_repeat_sites.csv](exact_Y_repeat_sites.csv), grouped overlaps in [repeat_annotation_summary.csv](repeat_annotation_summary.csv), and source hashes in [repeat_annotation_sources.csv](repeat_annotation_sources.csv). Reproduce with `python sgRNA_guide_development/code/annotate_five_guide_repeats.py` from the projects directory.

## Sources

- Rhie et al., 2023. [The complete sequence of a human Y chromosome](https://doi.org/10.1038/s41586-023-06457-y). Reference for DYZ3/S4CYH1L, DYZ1/HSat3, and TSPY composite-repeat nomenclature and architecture.
- [T2T reference annotation resource](https://github.com/marbl/CHM13).
- [Cen/Sat v2.1 BED](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed), retrieved 2026-09-16.
- [Composite repeats 2022DEC BED](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_composite-repeats_2022DEC.bed), retrieved 2026-09-16.
- Manz et al., 1993. [TSPY-related sequences represent a microheterogeneous gene family organized as constitutive elements in DYZ5 tandem repeat units on the human Y chromosome](https://pubmed.ncbi.nlm.nih.gov/8244388/).
