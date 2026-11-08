"""One-time caption/source-note update for linear main panel e."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]


def replace(rel,old,new):
    p=ROOT/rel;txt=p.read_text(encoding='utf-8')
    assert txt.count(old)==1,(rel,old)
    p.write_text(txt.replace(old,new),encoding='utf-8')


replace('Chapter4_figure_captions.md',
    '**e,** Endpoint abundance of all four Brunello guides per follow-up gene, as log₂(CPM + 1) after sample-wise normalisation across all 77,441 guides.',
    '**e,** Endpoint abundance of all four Brunello guides per follow-up gene, as CPM + 1 on a linear axis after sample-wise normalisation across all 77,441 guides. All four panels share the 0–50 scale; the faint reference line at 1 marks zero raw reads.')
replace('Chapter4_figure_captions.md',
    'Normalisation, colours and replicate symbols match main panel e; all endpoint guide panels share the same abundance axis.',
    'CPM normalisation, colours and replicate symbols match main panel e. These supplementary profiles retain a shared log₂(CPM + 1) axis; main e uses linear CPM + 1.')
replace('figures/chapter4/ch04_ht1080_candidate_discovery_notes.md',
    'For each sample, CPM = raw count / sum of all 77,441 guide counts × 10⁶. Plot log₂(CPM + 1).',
    'For each sample, CPM = raw count / sum of all 77,441 guide counts × 10⁶. Main panel e now plots linear CPM + 1 (15 October 2026); supplementary abundance displays retain log₂(CPM + 1).')
replace('figures/chapter4/ch04_ht1080_candidate_discovery_notes.md',
    'All endpoint panels share limits −0.17 to 6.25 (labelled ticks 0, 2, 4, 6); observed maximum 5.7184.',
    'Supplementary endpoint panels share log₂(CPM + 1) limits −0.17 to 6.25 (ticks 0, 2, 4, 6); observed maximum 5.7184. Main panel e uses shared linear CPM + 1 limits 0–50 (ticks every 10), with zero raw reads at 1 and a maximum plotted value of 44.1086. Its 64 values are in `source_data/ht1080_candidate_discovery/main_endpoint_cpm_plus_1.csv`.')
replace('figures/chapter4/ch04_ht1080_candidate_figures_plan.md',
    '| e | DDX3X, SREBF1, TADA2B and TAF6L guide profiles | Retain the existing four-guide endpoint displays and their sampling definitions. |',
    '| e | DDX3X, SREBF1, TADA2B and TAF6L guide profiles | Linear CPM + 1 with a common 0–50 axis, retaining all four guides and sampling definitions; zero raw reads plot at 1. |')
notes=ROOT/'figures/chapter4/ch04_ht1080_candidate_discovery_notes.md'
with notes.open('a',encoding='utf-8') as f:
    f.write('\n\n## Linear main-panel abundance, 15 October 2026\n\nAt the author’s request, main e now shows CPM + 1 on a linear axis. All four target panels share 0–50, and a faint line at 1 marks zero raw reads. CPM denominators, points, colours, replicate markers, guides and primary statistical results are unchanged. The new scale makes absolute abundance differences more prominent while compressing the lowest counts. The supplementary abundance plots retain their log₂ scales. Validation checks all 64 plotted values against CPM + 1 and the prior log₂ transform, the no-read values at 1, common limits and exported axis text. The previous main figure and scripts are preserved in `archive/ht1080_main2_cpm_plus1_20261015`.\n')
print('Updated main/supp captions, source notes and panel plan.')
