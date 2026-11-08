"""Record the 19 October 2026 adoption of paired networks into main Fig. 2d."""
from pathlib import Path
T=Path(__file__).resolve().parents[3]

def replace_once(text,old,new):
    assert text.count(old)==1,(old[:90],text.count(old))
    return text.replace(old,new,1)

def save(path,text):
    path.write_text(text,encoding='utf-8',newline='\n')

path=T/'Chapter4_figure_captions.md'
text=path.read_text(encoding='utf-8')
start=text.index('**d,** STRING v12.5 network of all 48 depleted candidates')
end=text.index(' Endpoint contrasts use two outgrowth cultures',start)
new='''**d,** Matched DepMap (left) and STRING (right) networks of all 57 primary candidates at directional FDR <0.25: 48 depleted and nine enriched. Circular node area is proportional to absolute drugZ normZ; purple filled circles denote depleted candidates and teal open circles enriched candidates. Positions and sizes match across the two views. The network area contains 39 candidates with a connection in at least one source; the identical lower strips retain the other 18. DepMap Public 25Q3 Pearson correlations at |r| >0.2 provide 41 edges across 1,186 models: 34 positive edges are solid blue and seven negative edges dashed rust, with width proportional to |r| plus a fixed minimum. STRING v12.5 functional associations at combined score ≥0.4 provide 82 grey edges, with width increasing with score; no connector proteins are added. Pale shading marks three communities obtained by weighted Louvain modularity optimisation of positive DepMap edges (resolution 1; seed 17), rather than complex membership. ZCCHC16* lacks a DepMap profile and is shown grey on the left. DepMap includes all recorded sexes and lineages, without a LOY restriction. The accompanying supplementary evidence matrix retains all candidate pairs, including 14 supported by both sources. These descriptive associations do not establish direct binding or LOY-specific synthetic lethality.'''
text=text[:start]+new+text[end:]
text=text.replace('\\*, DDX3X guide 4 and MED12 guide 3 have zero reads','In b–c, \\* marks DDX3X guide 4 and MED12 guide 3, which have zero reads',1)
start=text.index('## supp draft — Matched DepMap and STRING networks of HT1080 screen candidates')
end=text.index('\n## ',start+4)
new='''## supp draft — DepMap and STRING evidence matrix for HT1080 screen candidates

Companion to the matched networks in main Fig. 2d; final supplementary numbering pending.

![DepMap and STRING evidence matrix](figures/chapter4/ch04_supp_ht1080_depmap_string_evidence_matrix_draft.png)

Evidence for all 1,596 unique pairs among the 57 primary drugZ candidates (48 depleted and nine enriched at directional FDR <0.25). The lower triangle shows each pair once; the upper triangle and diagonal are omitted. Colours distinguish both sources (14 pairs), DepMap only (27), STRING only (68), neither above threshold (1,431), and unavailable DepMap data (56; hatched grey). Thresholds match Fig. 2d: DepMap Public 25Q3 Pearson |r| >0.2 and STRING v12.5 functional combined score ≥0.4. Minus signs identify the seven negative DepMap correlations. Genes are ordered by positive DepMap community, then remaining connected candidates and shared-strip candidates, ranked within groups by |normZ|. Block outlines mark the three positive-edge communities; the top strip shows signed drugZ normZ. ZCCHC16* lacks a DepMap profile, so its comparisons remain missing rather than being scored as zero. All retained DepMap edges use 1,186 models across recorded sexes and lineages, without a LOY restriction. These descriptive associations do not establish direct binding or LOY-specific mechanisms.

[PDF](figures/chapter4/ch04_supp_ht1080_depmap_string_evidence_matrix_draft.pdf) · [Editable SVG](figures/chapter4/ch04_supp_ht1080_depmap_string_evidence_matrix_draft.svg) · [Source and interpretation notes](figures/chapter4/ch04_ht1080_depmap_string_network_notes.md)
'''
text=text[:start]+new+text[end:]
save(path,text)

path=T/'draftv1/chapter4/chapter4_results.md'
text=path.read_text(encoding='utf-8')
text=replace_once(text,'d, the complete drugZ-candidate STRING network with circular nodes.',
    'd, matched DepMap (left) and STRING (right) networks for all 57 primary candidates, using identical circular-node positions and sizes.')
start=text.index('To place the candidate list in a functional context, we examined STRING associations among all 48 depleted candidates.')
end=text.index('\n\n',start)
new='''To place the candidate list in a functional context, we compared DepMap co-dependency correlations with STRING associations across all 57 primary candidates (**Fig. 2d**). DepMap yielded 41 correlations at |r| >0.2, including 34 positive and seven negative connections, whereas STRING yielded 82 functional associations at score ≥0.4. Five SAGA candidates—TAF6L, TADA2B, SUPT20H, TADA1 and TAF5L—grouped with MED12, STAG2 and CREBBP in one positive co-dependency community; ATXN7L3 remained outside this group. Fourteen candidate pairs met both source thresholds (**network-evidence supplementary figure**). These relationships contextualised the screen candidates but did not establish interactions in HT1080 or independently test LOY-specific dependency.

**Supplementary figure draft — DepMap and STRING evidence matrix.**

![DepMap and STRING evidence matrix](../../figures/chapter4/ch04_supp_ht1080_depmap_string_evidence_matrix_draft.png)

*Evidence for all 1,596 unique pairs among the 57 primary candidates, using the thresholds in Fig. 2d. The lower triangle shows each pair once; minus signs mark negative DepMap correlations, and hatching marks unavailable DepMap data. Supplementary numbering remains pending.* [Full caption](../../Chapter4_figure_captions.md#supp-draft--depmap-and-string-evidence-matrix-for-ht1080-screen-candidates) · [PDF](../../figures/chapter4/ch04_supp_ht1080_depmap_string_evidence_matrix_draft.pdf) · [Editable SVG](../../figures/chapter4/ch04_supp_ht1080_depmap_string_evidence_matrix_draft.svg).'''
text=text[:start]+new+text[end:]
save(path,text)

path=T/'figures/chapter4/ch04_ht1080_candidate_discovery_notes.md'
text=path.read_text(encoding='utf-8')
text=replace_once(text,'d, the full 48-candidate STRING network with circular nodes.',
    'd, paired DepMap and STRING networks for all 57 primary candidates, with matching circular-node positions and sizes.')
text=replace_once(text,'Physical size 9.6 × 8.8 inches (15 October 2026 revision).',
    'Physical size 9.6 × 10.6 inches (19 October 2026 revision); panels a–c retain their physical size and position relative to the page top.')
text=replace_once(text,'Network node colours use current drugZ depletion FDR classes.',
    'In main d, network node colour encodes depleted versus enriched direction, while circle area is proportional to absolute drugZ normZ; the FDR colour key applies to panel a.')
text=replace_once(text,'## Current full-candidate STRING network',
    '''## Current matched DepMap and STRING networks

On 19 October 2026, the author adopted the two matched networks into main Fig. 2d. Both show all 57 primary candidates: 48 depleted and nine enriched. Saved coordinates, areas and community assignments are reused from `source_data/ht1080_depmap_string_network_57/node_layout_and_modules.csv`; no analysis or layout is refitted. The left view retains 41 DepMap edges at |r| >0.2, and the right view 82 STRING v12.5 functional associations at score ≥0.4. Positions and circle areas match, with area 15 points² per unit |normZ|. Purple filled nodes denote depleted candidates; teal open nodes enriched candidates. Eighteen candidates with no edge in either source remain in both lower strips. Shading denotes positive DepMap communities, not SAGA membership. The [complete evidence matrix](ch04_supp_ht1080_depmap_string_evidence_matrix_draft.pdf) remains a separate supplementary figure. See [network notes](ch04_ht1080_depmap_string_network_notes.md) for provenance and interpretation. Reviewed originals were archived in `archive/ht1080_networks_into_fig2d_20261019` under the thesis root.

## Historical 48-candidate STRING network (superseded 19 October 2026)''')
text=replace_once(text,'Main panel d now incorporates all 48 primary drugZ depleted candidates',
    'The previous main panel d incorporated all 48 primary drugZ depleted candidates')
save(path,text)

path=T/'figures/chapter4/ch04_ht1080_depmap_string_network_notes.md'
text=path.read_text(encoding='utf-8')
start=text.index('Drafted 19 October 2026')
end=text.index('## Input and evidence')
text=text[:start]+'''Adopted 19 October 2026: the two matched networks now occupy main Fig. 2d, with DepMap on the left and STRING on the right. The complete evidence matrix remains a separate supplement for Results 4.3.1.2; final supplementary numbering is pending. The frozen candidate set, correlations, STRING scores, community assignments and coordinates are unchanged.

[Main Fig. 2 PDF](ch04_main_ht1080_candidate_discovery.pdf) · [Evidence-matrix PDF](ch04_supp_ht1080_depmap_string_evidence_matrix_draft.pdf) · [Evidence-matrix SVG](ch04_supp_ht1080_depmap_string_evidence_matrix_draft.svg)

![Matched networks in main Fig. 2d](review/ch04_main_ht1080_candidate_discovery_preview.png)

The [original three-panel draft](ch04_supp_ht1080_depmap_string_network_draft.pdf) is retained as a historical display. Reviewed files and captions were archived before replacement in `archive/ht1080_networks_into_fig2d_20261019` under the thesis root.

'''+text[end:]
text=replace_once(text,'Panel c shows each of the 1,596 unordered pairs once.',
    'The standalone supplementary matrix (formerly panel c) shows each of the 1,596 unordered pairs once in the lower triangle; the upper triangle and diagonal are omitted.')
text=replace_once(text,'- [Figure script](code/make_ht1080_depmap_string_network.py): reproducible communities, shared coordinates and editable plots.',
    '''- [Original figure script and shared drawing helpers](code/make_ht1080_depmap_string_network.py): reproducible communities, shared coordinates and editable plots.
- [Main Fig. 2 renderer](code/make_ht1080_candidate_layouts.py): paired networks in d, reusing frozen coordinates.
- [Standalone matrix renderer](code/make_ht1080_depmap_string_evidence_matrix.py): frozen evidence and ordering.''')
text=replace_once(text,'The thesis-style draft uses Open Sans, 11.2 ×14.4-inch pages, editable PDF/SVG text and a 600-dpi PNG.',
    'The current thesis-style main figure uses Open Sans on a 9.6 ×10.6-inch page; the standalone matrix is 11.2 ×8.0 inches. Both have editable PDF/SVG text and 600-dpi PNG exports. The original three-panel draft was 11.2 ×14.4 inches.')
save(path,text)
print('Updated captions, Results, and the two figure notes.')
