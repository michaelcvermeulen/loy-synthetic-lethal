"""Draft atlas of all selected H3K27/H3K4 signatures, preserving test scopes."""
from pathlib import Path
import html
import io
import json
import sys
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list

THESIS = Path(__file__).resolve().parents[3]
ANALYSIS = THESIS.parent / 'RNAseq/analysis/h3k27_h3k4_draft_2026_09_28'
SRC = THESIS / 'figures/chapter3/source_data/h3k27_h3k4_draft'
FIG = THESIS / 'figures/chapter3'
STEM = 'ch03_supp_draft_h3k27_h3k4_all_comparisons'
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style

SRC.mkdir(parents=True, exist_ok=True)
catalog = pd.read_csv(ANALYSIS / 'outputs/signature_catalogue.csv').fillna('')
manifest = pd.read_csv(ANALYSIS / 'outputs/contrast_manifest.csv').fillna('')
saved = pd.read_csv(ANALYSIS / 'outputs/saved_results.csv')
focused = pd.read_csv(ANALYSIS / 'outputs/focused_results.csv')
all_results = pd.concat([saved, focused], ignore_index=True)
assert not all_results.duplicated(['key', 'pathway']).any()
grid = pd.MultiIndex.from_product([catalog.gs_name, manifest.key], names=['pathway','key']).to_frame(index=False)
data = grid.merge(all_results, on=['key','pathway'], how='left', validate='one_to_one')
data = data.merge(manifest, on='key', how='left', validate='many_to_one')
data['status'] = data.status.fillna('not_in_saved_results')
data['display_score'] = data.NES
data['display_metric'] = 'NES'
cam = data.method == 'saved_CAMERA_0.01'
# CAMERA direction can differ from the descriptive mean logFC; use its own direction.
data.loc[cam, 'display_score'] = np.where(data.loc[cam,'camera_direction']=='Up', 1, -1) * -np.log10(data.loc[cam,'FDR'].clip(lower=1e-300))
data.loc[cam,'display_metric'] = 'CAMERA signed -log10(FDR)'
data.loc[data.status != 'tested', 'display_score'] = np.nan
data['display_limit'] = np.where(cam, 5., 4.)
data['display_clipped'] = data.display_score.abs() > data.display_limit
data['significant'] = (data.status == 'tested') & (data.FDR < .05)
assert len(data) == 87 * 52
assert data.loc[data.status=='tested','FDR'].between(0,1).all()
assert np.isfinite(data.loc[data.status=='tested','display_score']).all()
data.to_csv(SRC / 'complete_results.tsv', sep='\t', index=False)
catalog.to_csv(SRC / 'signature_catalogue.tsv', sep='\t', index=False)
manifest.to_csv(SRC / 'comparison_manifest.tsv', sep='\t', index=False)
summary = data.groupby(['key','label','block'],sort=False).agg(tested=('status',lambda s:sum(s=='tested')),
               significant=('significant','sum'), clipped_display_values=('display_clipped','sum')).reset_index()
summary.to_csv(SRC/'coverage_summary.tsv',sep='\t',index=False)

# Display-only filter: retain original statistics and their original FDR families.
manifest = manifest.loc[manifest.block != '2 Add-backs'].copy()
data = data.loc[data.key.isin(manifest.key)].copy()
coverage = data.groupby('pathway', sort=False).agg(
    tested_comparisons=('status', lambda s: sum(s == 'tested')),
    total_comparisons=('status', 'size')).reset_index()
coverage['tested_fraction'] = coverage.tested_comparisons / coverage.total_comparisons
coverage['included_in_display'] = coverage.tested_fraction >= .5
coverage.to_csv(SRC/'signature_display_coverage.tsv', sep='\t', index=False)
catalog = catalog.loc[catalog.gs_name.isin(coverage.loc[coverage.included_in_display, 'pathway'])].copy()
data = data.loc[data.pathway.isin(catalog.gs_name)].copy()
data.to_csv(SRC/'displayed_results.tsv', sep='\t', index=False)
manifest.to_csv(SRC/'displayed_comparison_manifest.tsv', sep='\t', index=False)

CMAP = LinearSegmentedColormap.from_list('signed', ['#3C5488','#FFFFFF','#E64B35'])
CMAP.set_bad('#D9D9D9')
BLOCK_LABEL = {'1 Thesis models':'Thesis models', '2 Add-backs':'Add-backs',
               '3 External Y comparisons':'External Y comparisons', '4 Early RPTEC RNA calls':'Early RPTEC',
               '5 Established autosomal loss':'Autosomal loss', '6 p53 induction controls':'p53 induction'}
positions = {}; centers = []; x = 0
for block, frame in manifest.groupby('block',sort=False):
    start = x
    for key in frame.key: positions[key] = x; x += 1
    centers.append((start, x-1, BLOCK_LABEL[block])); x += 1
ncols = x-1
apply_style(8)


def figure_for(names, title, stem, pdf):
    height = max(7.0, 3.6 + len(names)*.23)
    fig = plt.figure(figsize=(18, height))
    bottom = 2.20/height
    ax = fig.add_axes([.315,bottom,.67,1-bottom-.105])
    matrix = np.full((len(names),ncols),np.nan)
    missing = []; points = []
    lookup = data.set_index(['pathway','key'])
    for i,name in enumerate(names):
        for key,pos in positions.items():
            r=lookup.loc[(name,key)]
            if r.status == 'tested':
                matrix[i,pos]=np.clip(r.display_score/r.display_limit,-1,1)
                if r.significant: points.append((pos,i))
            else: missing.append((pos,i))
    ax.imshow(matrix,cmap=CMAP,vmin=-1,vmax=1,aspect='auto',interpolation='none')
    for start,end,label in centers:
        if start: ax.axvspan(start-1.5,start-.5,color='white',zorder=2)
        ax.text((start+end)/2,1.013,label,transform=ax.get_xaxis_transform(),
                ha='center',va='bottom',fontsize=8,weight='semibold')
    if points:
        px,py=zip(*points); ax.scatter(px,py,s=5,c='#111111',zorder=3)
    if missing:
        mx,my=zip(*missing);ax.scatter(mx,my,s=6,c='#888888',marker='x',linewidth=.3)
    ax.set_xticks([positions[k] for k in manifest.key],manifest.label,rotation=65,ha='right',fontsize=7.8)
    labels=['\n'.join(textwrap.wrap(n.replace('_',' '),65)) for n in names]
    ax.set_yticks(range(len(names)),labels,fontsize=7.6)
    ax.tick_params(length=0)
    for sp in ax.spines.values():sp.set_visible(False)
    fig.suptitle('Draft | '+title, x=.02, y=.975,ha='left',fontsize=12,weight='semibold')
    fig.text(.02,.939,f'{len(names)} signatures × {len(manifest)} comparisons; rows clustered by enrichment pattern (average linkage; 1 − Pearson r)',fontsize=9)
    cax1=fig.add_axes([.025,.085,.105,.012]);cb=fig.colorbar(ScalarMappable(norm=Normalize(-4,4),cmap=CMAP),cax=cax1,orientation='horizontal',ticks=[-4,0,4])
    cb.ax.set_title('GSEA: NES',fontsize=8)
    fig.text(.315,.044,'Dot: FDR < 0.05 within the stated testing family. Grey ×: untested / unavailable. Colour limits: NES ±4.',fontsize=8)
    fig.text(.315,.023,'Saved broad-screen FDR: thesis models. New focused-panel FDR: external Y, RPTEC and autosomal/p53 comparisons.',fontsize=8)
    pdf.savefig(fig)
    for ext in ['png','svg']:
        buffer=io.BytesIO();fig.savefig(buffer,format=ext,dpi=180)
        (FIG/f'{stem}.{ext}').write_bytes(buffer.getvalue())
    plt.close(fig)


pages = [(['01 H3K27 marks and functions'],'H3K27 marks and functions','h3k27'),
         (['02 H3K4 marks and functions'],'H3K4 marks and functions','h3k4'),
         (['03 Both H3K4 and H3K27'],'Signatures containing both H3K4 and H3K27','both_marks'),
         (['04 PRC2 and H3K27-related signatures','05 H3K4 regulator signatures'],'PRC2 and H3K4 regulator signatures','regulators')]
clustered_catalog = []
cluster_pairs = []
cluster_order = []
for groups, title, suffix in pages:
    part = catalog.loc[catalog.panel_group.isin(groups)].copy()
    names = part.gs_name.tolist()
    values = data.pivot(index='pathway', columns='key', values='display_score').reindex(index=names, columns=manifest.key).to_numpy()
    distances = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            shared = np.isfinite(values[i]) & np.isfinite(values[j])
            assert shared.sum() >= 3, 'Insufficient shared comparisons for correlation'
            a, b = values[i, shared], values[j, shared]
            assert np.std(a) > 0 and np.std(b) > 0, 'Constant enrichment profile'
            distance = float(np.clip(1 - np.corrcoef(a, b)[0, 1], 0, 2))
            distances.append(distance)
            cluster_pairs.append(dict(panel=suffix, pathway_1=names[i], pathway_2=names[j],
                                      shared_comparisons=int(shared.sum()), distance=distance))
    tree = linkage(np.asarray(distances), method='average', optimal_ordering=True)
    pd.DataFrame(tree, columns=['child_1', 'child_2', 'distance', 'leaf_count']).to_csv(SRC/f'row_linkage_{suffix}.tsv', sep='\t', index=False)
    part = part.iloc[leaves_list(tree)].copy()
    # Keep the combined regulator panel in its clustered order in HTML too.
    part['panel_group'] = groups[0] if len(groups) == 1 else '04 PRC2 and H3K4 regulator signatures'
    clustered_catalog.append(part)
    cluster_order.extend(dict(panel=suffix, row_order=i + 1, pathway=n) for i, n in enumerate(part.gs_name))
catalog = pd.concat(clustered_catalog, ignore_index=True)
pd.DataFrame(cluster_pairs).to_csv(SRC/'row_clustering_distances.tsv', sep='\t', index=False)
pd.DataFrame(cluster_order).to_csv(SRC/'clustered_row_order.tsv', sep='\t', index=False)
with PdfPages(FIG/f'{STEM}.pdf') as pdf:
    for groups,title,suffix in pages:
        names=[r['pathway'] for r in cluster_order if r['panel']==suffix]
        figure_for(names,title,f'{STEM}_{suffix}',pdf)


def fmt(v):
    if pd.isna(v):return 'not available'
    return f'{float(v):.4g}'


def colour(v,limit):
    rgb=CMAP((np.clip(v/limit,-1,1)+1)/2)[:3]
    return '#'+''.join(f'{int(round(z*255)):02x}' for z in rgb)


head = '''<!doctype html><html><head><meta charset="utf-8"><title>Draft H3K27/H3K4 comparison atlas</title>
<style>body{font:14px Arial,sans-serif;color:#222;margin:20px}h1{font-size:24px}p{max-width:1100px;line-height:1.5}
input,select{font:inherit;padding:7px;margin:4px 12px 12px 0}.scroll{overflow:auto;max-height:74vh;border:1px solid #bbb}
table{border-collapse:separate;border-spacing:0}th,td{border-bottom:1px solid #eee}thead th{position:sticky;top:0;background:#f7f7f7;z-index:3;height:190px;min-width:33px;vertical-align:bottom;font-weight:normal;font-size:11px}
thead th span{writing-mode:vertical-rl;transform:rotate(180deg);white-space:nowrap;padding:8px}
thead th:first-child{left:0;z-index:4;min-width:460px;text-align:left;height:auto;vertical-align:bottom;padding:8px}
tbody th{position:sticky;left:0;z-index:2;background:white;min-width:460px;max-width:550px;padding:5px 8px;text-align:left;font-size:11px;font-weight:normal;overflow-wrap:anywhere}
td{text-align:center;height:26px;min-width:33px;cursor:pointer}td.missing{background:#ddd;color:#888}.boundary{border-left:3px solid #222}
.group th{background:#eee;font-weight:bold;padding-top:10px}.group td{background:#eee}.legend{display:flex;gap:25px;flex-wrap:wrap;margin:12px 0}.gradient{display:inline-block;width:160px;height:13px;background:linear-gradient(to right,#3C5488,white,#E64B35)}
#detail{white-space:pre-wrap;padding:12px;background:#f5f5f5;min-height:55px;line-height:1.5}.hidden{display:none}</style></head><body>
<h1>Draft H3K27 / H3K4 gene-set atlas</h1>
<p>DISPLAY_COUNTS. These are expression-enrichment results, not measurements of histone marks.
Positive values refer to the first group named in each comparison. Add-backs are excluded; signatures must be tested in at least half of the remaining comparisons (23 of 46). This display filter does not change the original FDR values or testing families.
ARPE-19 remains pending and is not included.</p>
<p>Thesis models retain their original broad-screen FDR. The other comparisons use newly calculated exploratory FDR across the original focused panel within each contrast. Statistical significance is not directly equivalent between those families. Shared clones, samples and controls are not independent experiments.</p>
<label>Find a signature <input id="search" placeholder="e.g. MIKKELSEN, H3K4, EED"></label>
<label>Comparisons <select id="block"><option value="">All comparisons</option>BLOCK_OPTIONS</select></label>
<label>Signatures <select id="category"><option value="">All groups</option>CATEGORY_OPTIONS</select></label>
<div class="legend"><span>GSEA: −4 <i class="gradient"></i> +4 NES</span><span>● FDR &lt; 0.05; × untested / unavailable</span></div>
<p>Rows are clustered within each of the four panels using average linkage and Pearson correlation distance (1 − r), calculated from the uncapped NES values in jointly tested comparisons. Untested values are excluded, not replaced by zero. Column order is unchanged; clustering is descriptive, with each displayed comparison equally weighted.</p>
<p>Hover or click a cell for exact values, source, test status and correction scope. Display scales are capped; exact values are retained. HCP/ICP/LCP denote high/intermediate/low-CpG promoter sets. Names and descriptions preserve the source resource, including the H327ME3 spelling.</p>
<div id="detail">Select a cell to inspect its result.</div><div class="scroll"><table><thead><tr><th>Signature</th>'''
head=head.replace('DISPLAY_COUNTS', f'{len(catalog)} signatures across {len(manifest)} comparisons')
head=head.replace('BLOCK_OPTIONS',''.join(f'<option>{html.escape(b)}</option>' for b in manifest.block.unique()))
head=head.replace('CATEGORY_OPTIONS',''.join(f'<option>{html.escape(b)}</option>' for b in catalog.panel_group.unique()))
body=[];last=None
for m in manifest.itertuples():
    boundary='boundary' if m.block!=last else ''
    body.append(f'<th class="{boundary}" data-block="{html.escape(m.block)}"><span>{html.escape(m.label)}</span></th>');last=m.block
body.append('</tr></thead><tbody>')
lookup=data.set_index(['pathway','key']);last_group=None
for c in catalog.itertuples():
    if c.panel_group!=last_group:
        body.append(f'<tr class="group" data-category="{html.escape(c.panel_group)}"><th>{html.escape(c.panel_group)}</th>')
        body.extend(f'<td data-block="{html.escape(m.block)}"></td>' for m in manifest.itertuples());body.append('</tr>');last_group=c.panel_group
    body.append(f'<tr data-name="{html.escape(c.gs_name.lower())}" data-category="{html.escape(c.panel_group)}"><th title="{html.escape(str(c.gs_description),quote=True)}">{html.escape(c.gs_name)}</th>')
    last=None
    for m in manifest.itertuples():
        r=lookup.loc[(c.gs_name,m.key)]
        boundary=' boundary' if m.block!=last else '';last=m.block
        tip=f'{c.gs_name}\n{m.label}\nStatus: {r.status}\nMetric: {r.display_metric}\nValue: {fmt(r.display_score)}\nNES: {fmt(r.NES)}\nP: {fmt(r.pval)}\nFDR: {fmt(r.FDR)}\nTested members: {fmt(r["size"])}\nCorrection: {m.fdr_scope}\nRanking: {m.rank}\nSource: {m.source}'
        missing=r.status!='tested'
        style='' if missing else f'background:{colour(r.display_score,r.display_limit)}'
        marker='×' if missing else ('●' if r.significant else '')
        body.append(f'<td class="{"missing" if missing else "tested"}{boundary}" data-block="{html.escape(m.block)}" title="{html.escape(tip,quote=True)}" style="{style}">{marker}</td>')
    body.append('</tr>')
end='''</tbody></table></div><p><a href="source_data/h3k27_h3k4_draft/complete_results.tsv">Complete numerical table</a> · <a href="source_data/h3k27_h3k4_draft/signature_catalogue.tsv">Signature definitions</a> · <a href="ch03_supp_draft_h3k27_h3k4_all_comparisons.pdf">Four-page PDF</a></p>
<script>function filter(){const q=document.querySelector('#search').value.toLowerCase(),b=document.querySelector('#block').value,c=document.querySelector('#category').value;
document.querySelectorAll('tbody tr').forEach(r=>{r.style.display=(!c||r.dataset.category===c)&&(!r.dataset.name||r.dataset.name.includes(q))?'':'none'});
document.querySelectorAll('[data-block]').forEach(e=>{e.style.display=!b||e.dataset.block===b?'':'none'});}
document.querySelectorAll('input,select').forEach(e=>e.addEventListener('input',filter));document.querySelectorAll('td[title]').forEach(e=>e.addEventListener('click',()=>document.querySelector('#detail').textContent=e.title));</script></body></html>'''
(FIG/f'{STEM}.html').write_text(head+''.join(body)+end,encoding='utf-8')
previous_validation = json.loads((SRC/'validation.json').read_text(encoding='utf-8')) if (SRC/'validation.json').exists() else {}
(SRC/'validation.json').write_text(json.dumps({**previous_validation, 'comparisons':len(manifest),'signatures':len(catalog),
        'cells':len(data),'tested':int((data.status=='tested').sum()),'missing_statuses':data.status.value_counts().to_dict(),
        'clipped_display_values':int(data.display_clipped.sum()),'pages':len(pages),
        'new_focused_comparisons':int((manifest.method=='new_focused_fgsea').sum()),
        'source_statistics_preserved':True,'addbacks_in_display':False,
        'minimum_tested_fraction':0.5,'removed_signatures':int((~coverage.included_in_display).sum()),
        'row_clustering':'within panel; average linkage; 1 - Pearson r on jointly tested uncapped NES; optimal leaf ordering',
        'minimum_pairwise_shared_comparisons':min(r['shared_comparisons'] for r in cluster_pairs),
        'FDR_recomputed_after_display_filter':False},indent=2)+'\n',encoding='utf-8')
print(json.dumps({'comparisons':len(manifest),'signatures':len(catalog),'tested_cells':int((data.status=='tested').sum()),'pdf_pages':len(pages)}))
