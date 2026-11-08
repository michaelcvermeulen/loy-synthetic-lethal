# Guide filtering from saved reference-match tables

## Procedure for this reanalysis

Candidate 20-nt guide sequences were filtered using saved reference-match counts, combining matches in both orientations. Guides were required to have at least 10 exact Y-chromosome matches and no exact matches to non-Y chromosomes. Among guides with available mismatch annotations, we retained those with no non-Y sites with one or two mismatches. Guides lacking a required mismatch count were classified as unassessed unless an observed criterion established rejection. Three- and four-mismatch counts and saved CFD summaries were retained as annotations. Retained guides were ranked by decreasing exact Y-match count, with alphabetical spacer order used to resolve ties. The same filters were applied to all candidates, including the named historical guides.

## Results and coverage of this run

The input contained **14,723 guides**. Exact-match specificity and Y-abundance filters retained **6,901**. Required mismatch counts were available for **215** of these guides and missing for **6,686**. The final decisions were **215 retained**, **6,686 unassessed** and **7,822 rejected** across the input table.

Run status: **partial_mismatch_coverage**. The number of guides with available annotations is a data-coverage count, not an additional biological filter. Historical shortlist overlap was 215 guides. Available annotations originate from that historical shortlist; overlap therefore does not independently reproduce the original selection from the full candidate universe.

## Scope for thesis writing

This text describes the new table-based reanalysis. It should not be presented as a recovered record of the original experimental guide-selection procedure. No genome search or new activity prediction was performed. CFD values are legacy annotations, not recalculated predictions. They are not filters or ranking criteria. Non-Y CFD already excludes chrY; do not subtract Y CFD.

Reference context: Saved exact-match tables from the T2T-labelled workflow; HG002 T2T-Y (CP086569.2) and CHM13 non-Y reference sequences.
