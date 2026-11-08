"""One-time synchronisation of the reviewed dot-plot/circular-STRING draft."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
def edit(rel,pairs):
    p=ROOT/rel;s=p.read_text(encoding='utf-8')
    for old,new in pairs:
        assert s.count(old)==1,(rel,old,s.count(old))
        s=s.replace(old,new)
    p.write_text(s,encoding='utf-8')

p=ROOT/'Chapter4_figure_captions.md';s=p.read_text(encoding='utf-8')
heading='## main — HT1080 screening recovers DDX3X and nominates additional dependencies'
start=s.index('**a,**',s.index(heading));end=s.index('\n\n[PDF]',start)
caption='''**a,** Genome-wide drugZ normZ scores for 19,112 genes, ordered by depletion rank on a logarithmic axis. Negative values indicate relative guide depletion in LOY endpoints; positive values indicate enrichment. At FDR <0.25, 48 genes are depleted and nine enriched, including 28 and two, respectively, at FDR <0.05. The pooled non-targeting-control category is omitted; original scores, ranks and FDRs are retained. **b,** Scores and depletion FDRs for the four established X-linked gametolog dependencies. **c,** Ranked dot plots with horizontal labels for the 28 depleted and two enriched hits at directional FDR <0.05, plus TAF5L (depletion FDR 0.171), retaining all six SAGA candidates. Depleted genes are ordered across the two columns by depletion rank, followed by the two enriched hits. Bold names identify SAGA components; both columns share the score axis. **d,** STRING v12.5 network queried with all 48 depleted candidates and no added proteins. At functional score ≥0.4, 26 proteins have 66 connections; 22 candidates without a retained connection within this set are listed alongside. Nodes are circles of equal size within the connected network, coloured by drugZ depletion FDR class. Dark solid edges have physical score ≥0.7; lighter solid edges have physical score 0.4–<0.7; dashed edges are other functional associations at score ≥0.4. Shading identifies the six SAGA components. Positions are illustrative; no enrichment test is shown. These database associations do not establish direct binding or LOY-specific synthetic lethality. **e,** Endpoint abundance of all four Brunello guides for each follow-up gene, shown as log₂(CPM + 1) after sample-wise normalisation across all 77,441 guides. Blue, CTRL; red, LOY. Circles and triangles denote outgrowth cultures R1 and R2 from one established population per genotype. †, DDX3X guide 4 has zero reads in all six samples. normZ is a standardised statistic, not a fold change; the primary endpoint contrast does not adjust for baseline guide differences.'''
s=s[:start]+caption+s[end:]
s=s.replace('Purple solid borders/circles indicate drugZ FDR <0.05; ochre dashed borders/squares indicate 0.05–<0.25.',
    'Purple circles indicate drugZ FDR <0.05; ochre circles indicate 0.05–<0.25. Connected nodes have equal area, with gene labels positioned beside them.')
p.write_text(s,encoding='utf-8')

edit('draftv1/chapter4/chapter4_results.md',[(
'a, Genome-wide drugZ ranking, with labels restricted to directional FDR <0.05 hits plus the SAGA candidates; b, the saved Metascape PPI network from the MAGeCK candidate list; c, the four established gametolog dependencies; d, selected SAGA and SREBF candidates; e, endpoint guide abundance for the four follow-up targets.',
'a, Genome-wide drugZ ranking; b, the four established gametolog dependencies; c, ranked dot plots for directional FDR <0.05 hits plus the SAGA candidates; d, the complete drugZ-candidate STRING network with circular nodes; e, endpoint guide abundance for the four follow-up targets.')])

edit('figures/chapter4/ch04_ht1080_candidate_discovery_notes.md',[
('a, genome-wide drugZ ranking with directional FDR <0.05 hits and the six SAGA candidates labelled; b, the saved supporting PPI network; c, four established gametologs; d, selected SAGA/SREBF candidates; e, endpoint guides for the four follow-up targets. Physical size 9.6 × 7.6 inches, reduced from 9.6 × 8.6 with a and b sharing the first row.',
 'a, compact genome-wide drugZ ranking; b, four established gametologs; c, ranked dot plots with horizontal gene labels and FDRs for 29 depleted and two enriched candidates; d, the full 48-candidate STRING network with circular nodes; e, endpoint guides for the four follow-up targets. Physical size 9.6 × 10.3 inches.'),
('Main c–d and candidate-support supplementary d show depletion FDR; candidate-support supplementary e shows enrichment FDR.',
 'Main b and candidate-support supplementary d show depletion FDR; main c uses the directional FDR for each gene; candidate-support supplementary e shows enrichment FDR.'),
('## Saved PPI network provenance',
 '## Current full-candidate STRING network\n\nMain panel d now incorporates all 48 primary drugZ depleted candidates using the frozen STRING v12.5 query. There are 26 connected proteins and 66 associations at functional score ≥0.4, with 22 isolated candidates displayed alongside. Of the retained edges, 45 also meet physical score ≥0.4, including 24 at ≥0.7. All markers are circles; connected-node areas are equal, and gene names are offset from the nodes. The [network source notes](ch04_ht1080_all_sl_string_notes.md) describe exact requests, mapping, edge evidence and the full-size standalone view. This display replaces the older seven-node panel described below.\n\n## Historical saved PPI network provenance (superseded display)'),
('No new STRING query, database-version claim, PPI enrichment calculation or candidate selection was performed. The canvas/caption identify the MAGeCK input so this cannot be mistaken for a complete network of all 48 drugZ candidates.',
 'That historical rendering used no new STRING query or enrichment calculation, and its caption identified the MAGeCK input. The active main figure now uses the full drugZ network documented above.'),
('the network input membership and all seven nodes/11 unique links,',
 'the historical seven-node input audit and the active STRING network containing all 48 candidates and 66 unique links,'),
('Checks require the selected 29 depleted and two enriched labels specifically in main panel a, with no overlapping gene-name text boxes,',
 'Checks require the selected 29 depleted and two enriched names in main panel c, all 48 network names in d, equal circular connected-node markers, and no overlapping PDF text spans,'),
('## Main panel a: FDR <0.05 plus SAGA labels, 14 October 2026',
 '## Ranked hit dot plots and circular STRING nodes, 14 October 2026')])
p=ROOT/'figures/chapter4/ch04_ht1080_candidate_discovery_notes.md';s=p.read_text(encoding='utf-8')
start=s.index("At the author's request, panel a labels only")
s=s[:start]+'''The author selected ranked dot plots as the replacement for dense leader-line labels. Panel c retains the same 28 depleted and two enriched genes with directional FDR <0.05, plus TAF5L (FDR 0.171), giving 29 depleted names and all six SAGA candidates. Gene names are horizontal, SAGA names are semibold, and the two columns use the same normZ limits with the original FDR beside each row. The original selected SAGA/SREBF score panel is consolidated into this display. Panel a retains the complete screen distribution without individual gene callouts; the four prespecified gametologs remain in b. Panel d incorporates the full STRING network with circular nodes. Guide data and supplements are unchanged.

Reproduce the main with `make_ht1080_candidates.py --main-only`, which delegates to [make_ht1080_candidate_layouts.py](code/make_ht1080_candidate_layouts.py). Reproduce the standalone network with [make_ht1080_all_sl_string.py](code/make_ht1080_all_sl_string.py). Active named-hit data are in [main_named_hits.csv](source_data/ht1080_candidate_discovery/main_named_hits.csv); `main2a_label_positions.csv` describes the superseded callout layout. The prior main and rectangular-node network are preserved in `archive/ht1080_main2_new_layout_circles_20261014`.
'''
p.write_text(s,encoding='utf-8')

edit('figures/chapter4/ch04_ht1080_all_sl_string_notes.md',[
('This is a standalone draft for review before integration into the discovery multipanel figure; the existing main panel b remains the historical seven-node Metascape display.',
 'The full network is now integrated into discovery main panel d, with this full-size standalone view retained. Both displays use circular nodes and offset gene labels; the historical seven-node Metascape display is superseded.'),
('used in the existing discovery figure.', 'used in the previous discovery figure.'),
("Node border colour reproduces the candidate figure's purple/ochre drugZ FDR classes. Solid versus dashed borders, and circle versus square symbols for isolated candidates, preserve class identification without colour. Labels identify all 48 genes once.",
 "Circular node fill reproduces the candidate figure's purple/ochre drugZ FDR classes. Connected nodes have equal area; the separate isolated-gene list uses smaller circular bullets. Gene names are offset from connected nodes, with borderless white label backgrounds for readability and editable vector text. Labels identify all 48 genes once. Colours and line types were inspected in greyscale."),
])

p=ROOT/'figures/chapter4/ch04_ht1080_candidate_figures_plan.md';s=p.read_text(encoding='utf-8')
start=s.index('| Panel | Proposed display |');end=s.index('## supp — Candidate-level support',start)
s=s[:start]+'''| Panel | Current display | Purpose and scope |
|---|---|---|
| a | Compact all-gene drugZ rank plot | Retain 19,112 genes and the original directional FDR colours. Named hits are presented separately in c, avoiding crowded callouts. |
| b | Four X-linked gametologs: DDX3X, EIF1AX, ZFX and RPS4X | Preserve the positive benchmark and nonsignificant prespecified genes, with saved normZ and depletion FDR. |
| c | Two columns of ranked dot plots with horizontal names and FDRs | Show 28 depleted and two enriched FDR <0.05 hits plus TAF5L, retaining all six SAGA candidates. Read depleted ranks down the first column then the second; enriched hits follow. Use identical score axes and semibold SAGA names. |
| d | Full 48-candidate STRING network with circular nodes | Display 26 connected proteins and all 66 functional links at score ≥0.4, distinguishing physical support. List all 22 isolated candidates alongside. No extra proteins or enrichment claim. |
| e | DDX3X, SREBF1, TADA2B and TAF6L guide profiles | Retain the existing four-guide endpoint displays and their sampling definitions. |

The author selected ranked dot plots on 14 October 2026. Panels a/b occupy the upper left and c the upper right; the full network spans the middle and guides occupy the bottom. The main is 9.6 × 10.3 inches. Genotype blue/red colours remain distinct from FDR-class colours. The earlier compact seven-node network and leader-line labels are archived.

'''+s[end:]
s=s.replace('Use the same CPM definition as main d.', 'Use the same CPM definition as main e.')
s=s.replace('Extend main d to SREBF2', 'Extend main e to SREBF2')
old='The standalone version is for layout review before replacing the compact seven-node panel b; existing discovery exports and inline Results remain unchanged.'
assert old in s
s=s.replace(old,'The network is now incorporated into main panel d with circular nodes; the standalone view is retained at full size. Central and inline panel descriptions match the new layout.')
p.write_text(s,encoding='utf-8')

edit('draftv1/chapter4/chapter4_results_source_notes.md',[
('Main panels a–e show all-gene drugZ ranking (all nine enriched candidates labelled), the saved supporting PPI network, the four established gametologs, selected SAGA/SREBF candidates and endpoint guide profiles for the four follow-up genes. Panels a/b share the top row and the figure is approximately 12% shorter.',
 'The author subsequently selected horizontal ranked dot plots. Main panels a–e now show the compact all-gene drugZ ranking, the four established gametologs, named hit scores/FDRs, the full drugZ-candidate STRING network and endpoint guide profiles. The named set contains 28 depleted and two enriched FDR <0.05 hits plus TAF5L, retaining all six SAGA candidates. The 9.6 × 10.3-inch layout gives the full network a middle row with circular nodes.'),
('The requested PPI panel is redrawn from saved Metascape data with STRING/BioGRID evidence. Its original input exactly matches the 23 MAGeCK depleted candidates at FDR <0.25, not the 48 drugZ candidates. All seven connected proteins also meet the drugZ threshold. All 11 saved links are retained; dashed versus solid distinguishes STRING physical score <0.7 versus ≥0.7. Do not call this a complete PPI analysis of the primary drugZ list or a new enrichment test.',
 'The active PPI panel uses a frozen STRING v12.5 query of all 48 drugZ depleted candidates. All inputs map one-to-one; 26 connected proteins have 66 functional associations at score ≥0.4, and 22 isolated candidates remain visible. Physical-network evidence is distinguished by line type/weight. No proteins were added and no enrichment test was performed. This replaces the earlier seven-node network based on 23 MAGeCK candidates. Exact query provenance and the legacy input audit are retained in the figure source notes.')])

edit('Chapter4_working_plan.md',[(
'The compact main places drugZ ranking beside a small saved PPI network, labels all nine enriched genes, and retains the gametolog/candidate scores and follow-up guides. The PPI input is the supporting MAGeCK candidate set, explicitly identified.',
'The main now separates a compact genome-wide ranking from horizontal ranked dot plots with FDRs for 28 depleted and two enriched FDR <0.05 hits plus TAF5L, retaining all six SAGA candidates. A middle row contains the full 48-candidate STRING network with circular nodes, including 22 isolated candidates. The four-gametolog comparison and follow-up guide profiles are retained.')])
print('Synced captions, Results figure description, provenance and panel plan.')
