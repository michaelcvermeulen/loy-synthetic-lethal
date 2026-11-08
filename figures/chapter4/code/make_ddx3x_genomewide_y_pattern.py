"""Compact, annotated XY/X0/XX heatmap for the frozen DDX3X pattern set."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import pymupdf as fitz

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS/'figures/chapter4'
DATA = OUT/'source_data/ddx3x_genomewide_y_pattern'
STEM = 'ch04_supp_ddx3x_xy_x0_xx_partner_heatmap_draft'
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style


def main():
    apply_style(9)
    source = pd.read_csv(DATA/'selected_partners.csv')
    manifest = json.loads((DATA/'manifest.json').read_text())
    assert len(source)==manifest['selected_genes'] and source.gene.is_unique
    assert source.display_rank.tolist()==list(range(1,len(source)+1))
    mat = source[['xy_r','loy_r','female_r']].to_numpy(float)
    assert np.isfinite(mat).all()
    limit = .5
    assert np.max(np.abs(mat))<=limit, 'Increase symmetric colour limits; do not silently clip.'
    width,height = 4.6,8.5
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_axes([.255,.055,.49,.86])
    cmap = plt.get_cmap('RdBu_r')
    norm = Normalize(-limit,limit)
    im = ax.imshow(mat,cmap=cmap,norm=norm,aspect='auto',interpolation='none')
    np.testing.assert_allclose(np.asarray(im.get_array()),mat,rtol=0,atol=0)
    labels = []
    for name,field in [('XY','xy_n_models'),('X0','loy_n_models'),('XX','female_n_models')]:
        counts = source[field]
        ntext = str(int(counts.iloc[0])) if counts.nunique()==1 else f'{int(counts.min())}–{int(counts.max())}'
        labels.append(f'{name}\nn = {ntext}')
    ax.set(xticks=np.arange(3),xticklabels=labels,yticks=np.arange(len(source)),yticklabels=source.gene)
    ax.xaxis.tick_top()
    ax.tick_params(axis='both',length=0,pad=5)
    ax.tick_params(axis='y',labelsize=9)
    ax.set_xticks(np.arange(-.5,3,1),minor=True)
    ax.set_yticks(np.arange(-.5,len(source),1),minor=True)
    ax.grid(which='minor',color='white',linewidth=1)
    ax.tick_params(which='minor',bottom=False,left=False,top=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    reference = set(source.loc[source.reference_partner,'gene'])
    for label,gene in zip(ax.get_yticklabels(),source.gene):
        if gene in reference:
            label.set_fontweight('semibold')
    annotations = []
    for row in range(len(source)):
        for col in range(3):
            val = mat[row,col]
            rounded = round(float(val),2)
            label = f'{0. if rounded==0 else rounded:.2f}'
            rgb = cmap(norm(val))[:3]
            luminance = sum(w*c for w,c in zip([.2126,.7152,.0722],rgb))
            colour = 'white' if luminance<.49 else '#252525'
            ax.text(col,row,label,ha='center',va='center',fontsize=8.5,color=colour)
            annotations.append({'gene':source.gene.iloc[row],'group':['XY','X0','XX'][col],
                                'value':float(val),'printed':label})
    cax = fig.add_axes([.80,.33,.028,.30])
    cb = fig.colorbar(im,cax=cax,ticks=[-.5,-.25,0,.25,.5])
    cb.set_label('DDX3X correlation (adjusted r)',fontsize=9,labelpad=5)
    cb.ax.tick_params(labelsize=8,width=.6,length=2)
    cb.outline.set_linewidth(.5)
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{STEM}.{ext}',dpi=600 if ext=='png' else 300)
    fig.savefig(OUT/'review'/f'{STEM}_preview.png',dpi=180)
    plt.close(fig)
    with fitz.open(OUT/f'{STEM}.pdf') as doc:
        page = doc[0]
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b
                 for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        outside = [s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        overlaps = []
        for i,a in enumerate(spans):
            for b in spans[i+1:]:
                ra,rb = fitz.Rect(a['bbox']),fitz.Rect(b['bbox'])
                area = (ra&rb).get_area()
                if area>1 and area/min(ra.get_area(),rb.get_area())>.1:
                    overlaps.append([a['text'],b['text']])
        assert not outside,outside
        assert not overlaps,overlaps
        page.get_pixmap(matrix=fitz.Matrix(2,2),colorspace=fitz.csGRAY).save(
            OUT/'review'/f'{STEM}_greyscale.png')
    record = dict(stem=STEM,dimensions_inches=[width,height],font='Open Sans',
                  n_genes=len(source),n_cells=len(annotations),cell_values=annotations,
                  colour_scale={'palette':'RdBu_r','limits':[-limit,limit],'clipped_values':0,
                                'row_standardisation':False},
                  row_order=manifest['order'],bold_genes=sorted(reference),
                  validation={'all_matrix_values_match':True,'text_outside_page':outside,'text_overlaps':overlaps})
    (DATA/'render_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k!='cell_values'},indent=2))


if __name__=='__main__':
    main()
