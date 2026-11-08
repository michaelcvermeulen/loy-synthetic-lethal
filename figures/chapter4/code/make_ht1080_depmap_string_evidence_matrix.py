"""Standalone evidence matrix retained after adopting the networks in Fig. 2d."""
from pathlib import Path
import json, hashlib
import make_ht1080_depmap_string_network as N
import matplotlib.pyplot as plt
import pymupdf as fitz

STEM='ch04_supp_ht1080_depmap_string_evidence_matrix_draft'


def main():
    N.apply_style(8)
    groups,module,core,strip,order,pos=N.load_saved_structure()
    fig=plt.figure(figsize=(11.2,8.0))
    N.matrix_panel(fig,order,groups,rect=[.12,.145,.83,.801],bar_rect=[.12,.955,.83,.018])
    fig.text(.12,.024,'* ZCCHC16: no DepMap profile. Upper triangle and diagonal omitted; each gene pair is shown once.',
             fontsize=8,color='#555555')
    for ext in ['pdf','svg','png']:
        fig.savefig(N.OUT/f'{STEM}.{ext}',dpi=600 if ext=='png' else 300)
    fig.savefig(N.OUT/'review'/f'{STEM}_preview.png',dpi=140)
    plt.close(fig)
    with fitz.open(N.OUT/f'{STEM}.pdf') as doc:
        page=doc[0]
        assert all(g in page.get_text() for g in N.NODES.index)
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans=[s for block in page.get_text('dict')['blocks'] if 'lines' in block
               for line in block['lines'] for s in line['spans'] if s['text'].strip()]
        outside=[s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        assert not outside,outside
        page.get_pixmap(matrix=fitz.Matrix(1.25,1.25),colorspace=fitz.csGRAY).save(N.OUT/'review'/f'{STEM}_greyscale.png')
    record=dict(stem=STEM,page_inches=[11.2,8.0],png_dpi=600,genes=57,unique_pairs=1596,
                both=14,depmap_only=27,string_only=68,negative_correlations=7,
                missing_depmap_pairs=56,paired_network_destination='Main Fig. 2d',
                source='Frozen pairwise evidence and matrix order; no refit',text_outside_page=outside)
    record['sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        [Path(__file__),N.DATA/'pairwise_evidence.csv',N.DATA/'node_layout_and_modules.csv',N.DATA/'panel_c_evidence_codes.csv']+
        [N.OUT/f'{STEM}.{ext}' for ext in ['pdf','svg','png']]}
    (N.DATA/'matrix_figure_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
