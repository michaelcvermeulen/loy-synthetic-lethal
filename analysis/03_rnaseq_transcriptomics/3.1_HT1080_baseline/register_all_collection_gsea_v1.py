"""Register new all-collection figures without replacing the three-collection drafts."""
from pathlib import Path
from datetime import datetime
import json,hashlib,shutil
import pandas as pd
ROOT=Path(__file__).resolve().parent;BASE=ROOT/'outputs/programme_analysis_v1';OUT=BASE/'all_collections'
THESIS=ROOT.parents[3]/'SEPT_11_2026_THESIS';FIG=THESIS/'figures/chapter3'
s=json.loads((OUT/'summary.json').read_text());stamp=datetime.now().strftime('%Y%m%d_%H%M%S')
archive=THESIS/'archive'/f'ht1080_all_collection_gsea_{stamp}';archive.mkdir(parents=True)
caption=THESIS/'Chapter3_figure_captions.md';manifest=FIG/'figure_manifest.tsv'
for p in [caption,manifest]:shutil.copy2(p,archive/p.name)
current=caption.read_text(encoding='utf-8')
main='ch03_main_ht1080_rna_programmes_regulators_all_collections';atlas='ch03_supp_ht1080_gsea_all_collections_atlas'
assert f'**Figure:** `{main}`' not in current
current+=f'''

## main - HT1080 programme and regulator patterns across all MSigDB collections

**Figure:** `{main}` | **Status:** all-collection alternative main draft; three-collection version retained

[PDF](figures/chapter3/{main}.pdf) | [SVG](figures/chapter3/{main}.svg) | [PNG](figures/chapter3/{main}.png)

**a,** Selected GSEA patterns from the full frozen MSigDB 2026.1.Hs H and C1-C9 resource. Within each collection and enrichment direction, the pathway with the smallest overall global FDR among directionally consistent discoveries is displayed ({s['displayed_main_pathways']} rows). Selection uses FDR then lexical pathway order to break ties; it is a display rule, not biological prioritisation. Colours indicate normalised enrichment scores (NES), and stars indicate global FDR <0.05 in the indicated contrast. Overall FDR is printed beside each row. **b,** The existing 12 smallest-overall-FDR directionally consistent CollecTRI regulator patterns; this panel is unchanged by the expanded GSEA screen and uses its own TF testing family and score scale.

GSEA uses signed moderated t rankings from the current overall batch + condition model and the August/May effects of the condition-by-period model. The primary ranking comprises 25,461 unique autosomal/non-PAR X symbols. Sets with 15-500 represented symbols were tested; BH correction spans all {s['tested_overall']:,} eligible sets across every collection, separately per contrast. There were {s['significant_overall']:,} overall discoveries, of which {s['directionally_consistent']:,} had matching NES directions in both periods and {s['both_periods_significant']:,} were also significant in each period. Positional, perturbation, disease and cell-type signatures describe their source definitions, not demonstrated mechanisms or identities in HT1080. Long labels are abbreviated for display; exact identifiers are retained in source tables. The three-collection CAMERA sensitivity does not constitute a correlation-aware test of this expanded family. GSEA remains exploratory.

Collection codes: H, Hallmark; C1, positional; C2, curated; C3, regulatory targets; C4, computational; C5, ontology; C6, oncogenic; C7, immunologic; C8, cell type; C9, cell perturbation.

## supp - Complete HT1080 all-collection GSEA significant-set atlas

**Figure:** `{atlas}` | **Status:** {s['atlas_pages']}-page supplementary atlas

[PDF](figures/chapter3/{atlas}.pdf)

All {s['significant_overall']:,} pathways passing overall global FDR <0.05 are shown, ordered by collection and overall FDR. Each row displays overall, August and May NES on a common colour scale across pages. Stars mark global FDR <0.05 in each contrast, and the overall FDR is printed beside the row. All overall discoveries are retained, including those without period-direction agreement; per-period statistics are therefore shown without a persistence filter. Full pathway identifiers, leading-edge genes, collection metadata, and page assignments are in the source tables. All tested pathways, including nonsignificant sets, are available in the searchable all-collection report and full CSV results.
'''
caption.write_text(current,encoding='utf-8')
m=pd.read_csv(manifest,sep='\t').fillna('');rows=[]
for record in s['exports']:
 name=record['file'];role='supp' if name.startswith('ch03_supp') else 'main'
 rows.append(dict(file=name,role=role,title=role+' - HT1080 all-collection GSEA',guide_set='',
 status='all-collection extension; original three-collection figures retained',
 original_source='RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/plot_all_collection_gsea_v1.py',sha256=record['sha256']))
assert not m.file.isin([r['file'] for r in rows]).any()
pd.concat([m,pd.DataFrame(rows)],ignore_index=True).to_csv(manifest,sep='\t',index=False)
for stem in [main,atlas]:
 (FIG/f'{stem}_notes.md').write_text(f'''# {stem}

Active captions: [central caption file](../../Chapter3_figure_captions.md).

Reproduce with Rscript --vanilla run_all_collection_gsea_v1.R, then Python plot_all_collection_gsea_v1.py in RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE. Source files remain unchanged; source_data/ht1080_rna_all_collections contains full overall results, displayed rows, atlas page index and export hashes. Raw tests for five contrasts, global/within-collection BH, warnings and validation are in outputs/programme_analysis_v1/all_collections. Main figure is 183 x 245 mm, Open Sans, PDF/SVG and 600-dpi PNG; atlas is 210 x 260 mm per page with editable vector text. Long labels are abbreviated in both figures; full names remain in tables and the report.
''',encoding='utf-8')
for p in [ROOT/'README.md',BASE/'README.md']:
 old=p.read_text(encoding='utf-8');shutil.copy2(p,archive/(('analysis_' if p.parent==ROOT else 'programme_')+p.name))
 link='outputs/programme_analysis_v1/all_collections/report.html' if p.parent==ROOT else 'all_collections/report.html'
 p.write_text(f'## Expanded all-collection GSEA\n\n[All-collection results and figures]({link}): {s["tested_overall"]:,} eligible sets; {s["significant_overall"]:,} overall global-FDR discoveries. The earlier three-collection analysis is retained.\n\n'+old,encoding='utf-8')
plan=THESIS/'research/molecular_characterisation/ht1080_rna_results_figure_selection.md'
old=plan.read_text(encoding='utf-8');shutil.copy2(plan,archive/plan.name)
plan.write_text(old+f'\n\n## All-collection GSEA extension\n\nAuthor requested expansion to all MSigDB pathways. New alternative main figure: {main}; supplementary atlas: {atlas} ({s["atlas_pages"]} pages). All {s["tested_overall"]:,} eligible H/C1-C9 sets use global BH per contrast. Earlier three-collection figures and their Results numbers remain labelled to that scope; do not mix counts between screens. The TF panel is unchanged. Review the expanded signatures before replacing the pathway narrative.\n',encoding='utf-8')
files=['run_all_collection_gsea_v1.R','plot_all_collection_gsea_v1.py','register_all_collection_gsea_v1.py']
(OUT/'code_hashes.json').write_text(json.dumps({f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},indent=2))
print('All-collection figures, captions, manifest and report entry points registered.')
