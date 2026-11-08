"""Layout checks and colour-vision previews for the standalone XY partner draft."""
from pathlib import Path
import json
import numpy as np
from PIL import Image
import pymupdf as fitz

OUT=Path(__file__).resolve().parents[1]
STEM='ch04_supp_xy_paralog_codependency_partners_draft'
with fitz.open(OUT/f'{STEM}.pdf') as doc:
    page=doc[0]
    spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
    overlaps=[]
    for i,a in enumerate(spans):
        for b in spans[i+1:]:
            ra,rb=fitz.Rect(a['bbox']),fitz.Rect(b['bbox'])
            area=(ra&rb).get_area()
            if area>1 and area/min(ra.get_area(),rb.get_area())>.10:
                overlaps.append([a['text'],b['text'],area])
    assert not overlaps,overlaps
    page.get_pixmap(matrix=fitz.Matrix(1.8,1.8),clip=fitz.Rect(0,page.rect.height-120,page.rect.width,page.rect.height)).save(OUT/'review'/f'{STEM}_footer.png')

# Same Machado severity-100 matrices and linear-light transform used by the
# bundled thesis-style preview_condition_colours.py; display-gamut clipping.
matrices={
    'deuteranopia':[[.367322,.860646,-.227968],[.280085,.672501,.047413],[-.011820,.042940,.968881]],
    'protanopia':[[.152286,1.052583,-.204868],[.114503,.786281,.099216],[-.003882,-.048116,1.051998]],
    'tritanopia':[[1.255528,-.076749,-.178779],[-.078411,.930809,.147602],[.004733,.691367,.303900]],
}
rgb=np.asarray(Image.open(OUT/'review'/f'{STEM}_preview.png').convert('RGB'),dtype=float)/255
linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
for label,matrix in matrices.items():
    transformed=np.clip(linear@np.array(matrix).T,0,1)
    encoded=np.where(transformed<=.0031308,12.92*transformed,1.055*transformed**(1/2.4)-.055)
    Image.fromarray(np.rint(encoded*255).astype('uint8')).save(OUT/'review'/f'{STEM}_{label}.png')
record={'pdf_text_overlaps':overlaps,'colour_vision_previews':list(matrices),'simulation':'Machado severity 100, linear-light sRGB, display-gamut clipping','greyscale_limitation':'Magnitude is visible; correlation sign relies on colour or source values.'}
(OUT/'source_data/xy_paralog_codependency/visual_validation.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps(record,indent=2))
