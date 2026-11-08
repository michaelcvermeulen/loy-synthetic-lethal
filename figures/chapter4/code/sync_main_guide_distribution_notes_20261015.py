"""One-time synchronisation of the reviewed main-figure display descriptions."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]


def replace(path, old, new):
    p=ROOT/path
    text=p.read_text(encoding='utf-8')
    assert text.count(old)==1,(path,old)
    p.write_text(text.replace(old,new),encoding='utf-8')


caption_path=ROOT/'Chapter4_figure_captions.md'
text=caption_path.read_text(encoding='utf-8')
heading='## main — HT1080 screening recovers DDX3X and nominates additional dependencies'
start=text.index('**a,**',text.index(heading))
end=text.index('\n\n[PDF]',start)
assert '**c,** Ranked dot plots' in text[start:end]
caption=(
    '**a,** Genome-wide drugZ normZ scores for 19,112 genes, ordered by depletion rank on a logarithmic axis. '
    'At FDR <0.25, 48 genes are depleted and nine enriched, including 28 and two, respectively, at FDR <0.05. '
    'Small grey labels identify DDX3X, TIMM8A, SREBF2 and TAF6L. The pooled non-targeting-control category is omitted; original scores, ranks and FDRs are retained. '
    '**b,** Guide-level contrasts for the four established X-linked gametolog dependencies. '
    '**c,** The same guide display for the 28 depleted and two enriched hits at directional FDR <0.05, plus TAF5L (depletion FDR 0.171), retaining all six SAGA candidates. '
    'Genes follow depletion rank across the two columns, followed by the enriched hits; bold names identify SAGA components. '
    'In b–c, equal-height black ticks show saved MAGeCK endpoint LOY-versus-CTRL guide log₂ fold changes, over the same grey empirical density of 999 non-targeting guides with endpoint reads. '
    'Negative values indicate relative depletion in LOY. Axes are shared, without horizontal jitter; nearby guides can overlap. '
    'One non-targeting guide with no endpoint reads is excluded from the reference curve. MAGeCK used non-targeting normalisation, so centring of this curve is not independent QC. '
    'Values beside genes are the original drugZ depletion FDRs, or enrichment FDRs for ATP6V0E1 and USP7; guide fold changes are distinct from gene-level normZ. '
    '**d,** STRING v12.5 network of all 48 depleted candidates, without added proteins. At functional score ≥0.4, 26 proteins have 66 connections; 22 unconnected candidates are listed alongside. '
    'Circular nodes are coloured by drugZ depletion FDR class. Dark solid edges have physical score ≥0.7; lighter solid edges have physical score 0.4–<0.7; dashed edges are other functional associations at score ≥0.4. '
    'Shading identifies the six SAGA components. Positions are illustrative; no enrichment test is shown, and associations do not establish direct binding or LOY-specific synthetic lethality. '
    '**e,** Endpoint abundance of all four Brunello guides per follow-up gene, as log₂(CPM + 1) after sample-wise normalisation across all 77,441 guides. '
    'Blue, CTRL; red, LOY; circles and triangles, outgrowth cultures R1 and R2 from one established population per genotype. '
    '†, DDX3X guide 4 and MED12 guide 3 have zero reads in all six samples; their fold-change ticks are omitted in b–c. '
    'Guide marks are not biological replicates. The endpoint contrasts do not adjust for baseline guide differences.'
)
caption_path.write_text(text[:start]+caption+text[end:],encoding='utf-8')

replace('draftv1/chapter4/chapter4_results.md',
    '*a, Genome-wide drugZ ranking; b, the four established gametolog dependencies; c, ranked dot plots for directional FDR <0.05 hits plus the SAGA candidates; d, the complete drugZ-candidate STRING network with circular nodes; e, endpoint guide abundance for the four follow-up targets. Numbering is provisional.*',
    '*a, Genome-wide drugZ ranking with four labelled hits; b, guide fold changes for the four established gametologs over the non-targeting reference distribution; c, the same guide display for directional FDR <0.05 hits plus the SAGA candidates, retaining drugZ FDRs; d, the complete drugZ-candidate STRING network with circular nodes; e, endpoint guide abundance for the four follow-up targets. Numbering is provisional.*')

replace('figures/chapter4/ch04_ht1080_candidate_figures_plan.md',
    '| a | Compact all-gene drugZ rank plot | Retain 19,112 genes and the original directional FDR colours. Named hits are presented separately in c, avoiding crowded callouts. |',
    '| a | Compact all-gene drugZ rank plot | Retain 19,112 genes and the original directional FDR colours. Small grey labels identify DDX3X, TIMM8A, SREBF2 and TAF6L. |')
replace('figures/chapter4/ch04_ht1080_candidate_figures_plan.md',
    '| b | Four X-linked gametologs: DDX3X, EIF1AX, ZFX and RPS4X | Preserve the positive benchmark and nonsignificant prespecified genes, with saved normZ and depletion FDR. |',
    '| b | Guide distributions for DDX3X, EIF1AX, ZFX and RPS4X | Equal-height ticks show saved MAGeCK guide log₂FC over the empirical non-targeting reference. Retain original drugZ depletion FDRs beside each gene. |')
replace('figures/chapter4/ch04_ht1080_candidate_figures_plan.md',
    '| c | Two columns of ranked dot plots with horizontal names and FDRs | Show 28 depleted and two enriched FDR <0.05 hits plus TAF5L, retaining all six SAGA candidates. Read depleted ranks down the first column then the second; enriched hits follow. Use identical score axes and semibold SAGA names. |',
    '| c | Two columns of guide distributions with horizontal names and drugZ FDRs | Keep the same 28 depleted and two enriched FDR <0.05 hits plus TAF5L, retaining all six SAGA candidates. Preserve rank order; use common guide-LFC axes and the same NTC density as b. SAGA names remain semibold. |')
replace('figures/chapter4/ch04_ht1080_candidate_figures_plan.md',
    'The author selected ranked dot plots on 14 October 2026. Panels a/b occupy the upper left and c the upper right; the full network spans the middle and guides occupy the bottom. The main is 9.6 × 10.3 inches. Genotype blue/red colours remain distinct from FDR-class colours. The earlier compact seven-node network and leader-line labels are archived.',
    'On 15 October 2026 the author requested replacing the b/c dot plots with guide distributions while retaining the genes. Panels a/b occupy the upper left and c the upper right; the full network spans the middle and endpoint abundance remains at the bottom. The main is 9.6 × 11.3 inches to keep the denser guide rows readable. Neutral guide ticks have equal height and no horizontal jitter; the grey curve uses 999 non-targeting guides with endpoint reads. Genotype blue/red colours remain distinct from FDR-class colours. The previous dot-plot main is archived in `archive/ht1080_main2_guide_distributions_20261015`.')

replace('figures/chapter4/ch04_ht1080_candidate_discovery_notes.md',
    '- [Main PDF](ch04_main_ht1080_candidate_discovery.pdf): a, compact genome-wide drugZ ranking; b, four established gametologs; c, ranked dot plots with horizontal gene labels and FDRs for 29 depleted and two enriched candidates; d, the full 48-candidate STRING network with circular nodes; e, endpoint guides for the four follow-up targets. Physical size 9.6 × 10.3 inches.',
    '- [Main PDF](ch04_main_ht1080_candidate_discovery.pdf): a, genome-wide drugZ ranking with four small grey gene labels; b, guide distributions for the four established gametologs; c, guide distributions for the same 29 depleted and two enriched candidates; d, the full 48-candidate STRING network with circular nodes; e, endpoint abundance for the four follow-up targets. Panels b/c retain original drugZ FDRs beside saved MAGeCK guide fold changes. Physical size 9.6 × 11.3 inches (15 October 2026 revision).')

notes=ROOT/'figures/chapter4/ch04_ht1080_candidate_discovery_notes.md'
with notes.open('a',encoding='utf-8') as f:
    f.write('''

## Current main guide distributions, 15 October 2026

At the author's request, b and c now show saved MAGeCK endpoint guide log₂ fold changes over a repeated empirical NTC density. All prior genes and their order are retained: DDX3X/EIF1AX/ZFX/RPS4X in b; 29 depleted and two enriched candidates in c. Original drugZ FDRs remain beside the rows. Panel a has small grey labels for DDX3X, TIMM8A, SREBF2 and TAF6L. The full STRING network and endpoint-abundance measurements are retained. The preceding dated dot-plot description is historical.

The display contains 34 unique genes and 136 unique library guides. DDX3X guide 4 and MED12 guide 3 have no reads in all six samples; their software-reported LFCs are preserved in the source table but omitted from the figure. There are 15 guide ticks in b and 122 in c (134 unique observed guides, with three DDX3X values repeated between panels). Other guides are not filtered by abundance or significance. The ticks have equal height; height has no biological meaning, and nearby x positions may overlap. There is no horizontal jitter. Shared limits of −8.5 to +5.8 contain every observed selected-guide and NTC value.

The common grey curve uses the 999 of 1,000 non-targeting guides with at least one endpoint read. It is a Gaussian KDE with bandwidth 0.0750263136, using `0.9 × min(sample SD, IQR/1.34) × n^(-1/5)`, scaled to a common peak. The one zero-endpoint-read NTC is excluded only from this display; original analyses are unchanged. This empirical reference is not a newly fitted statistical null, and its centring is not independent QC because MAGeCK used NTC normalisation. Guide marks are observations from the same experiment, not additional biological replicates.

Reproduce with `make_ht1080_candidates.py --main-only`, then `validate_ht1080_candidates.py`. The main renderer calls [ht1080_main_guide_distributions.py](code/ht1080_main_guide_distributions.py). [Source tables and provenance](source_data/ht1080_main_guide_distributions/) include all selected guide counts and saved statistics, all NTCs including the excluded one, the KDE grid, exact displayed tick coordinates, original input hashes and display rules. Validation checks the retained gene sets, exact guide/gene statistics, original count joins, equal tick heights, all 137 displayed marks, 999 eligible NTCs, axis inclusion and density normalisation, alongside existing export/network checks. Open Sans PDF/SVG/600-dpi PNG exports and colour/greyscale previews were inspected; no text overlaps or clipping were found. Previous main exports, source scripts and documentation are in `archive/ht1080_main2_guide_distributions_20261015`.
''')

replace('draftv1/chapter4/chapter4_results_source_notes.md',
    'The author subsequently selected horizontal ranked dot plots. Main panels a–e now show the compact all-gene drugZ ranking, the four established gametologs, named hit scores/FDRs, the full drugZ-candidate STRING network and endpoint guide profiles.',
    'The author subsequently selected horizontal ranked dot plots, then requested guide distributions in b/c on 15 October 2026. Main panels a–e now show the compact all-gene drugZ ranking with four small grey labels; guide distributions for the four established gametologs; the same guide display for the named hits with original drugZ FDRs; the full drugZ-candidate STRING network; and endpoint guide abundance. Equal-height ticks show saved MAGeCK log₂ fold changes over the empirical density of 999 represented non-targeting guides, with shared axes and zero-read guides flagged.')
replace('draftv1/chapter4/chapter4_results_source_notes.md',
    'The 9.6 × 10.3-inch layout gives the full network a middle row with circular nodes.',
    'The current 9.6 × 11.3-inch layout gives the full network a middle row with circular nodes. Guide-display provenance is in `figures/chapter4/source_data/ht1080_main_guide_distributions`; saved analyses are unchanged.')

replace('Chapter4_working_plan.md',
    'The main now separates a compact genome-wide ranking from horizontal ranked dot plots with FDRs for 28 depleted and two enriched FDR <0.05 hits plus TAF5L, retaining all six SAGA candidates.',
    'Updated 15 October 2026: the main separates a compact genome-wide ranking, with four small grey labels, from guide-level fold changes over the non-targeting reference distribution in b/c. These panels retain the same four gametologs and 28 depleted plus two enriched FDR <0.05 hits plus TAF5L, including all six SAGA candidates; original drugZ FDRs remain beside each row.')
print('Updated central caption, Results figure summary, source notes and working plans.')
