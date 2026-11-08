"""Freeze the seven author-selected add-back figures and update active links.

8 October 2026. Copies existing artwork; no statistical analysis. Existing
selections and unrelated manuscript sections are preserved. Selection notes
link to the authoritative central captions rather than duplicating them.
"""
from pathlib import Path
import csv
import hashlib
import io
import json
import re
import shutil

raise SystemExit('Historical seven-figure selection script: S03 was consolidated into S02 on 9 October 2026. Use render_addback_signature_final.py --supp-only and current selection notes.')

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
SELECTED = FIG / 'selected'
ENTRIES = [
    ('F01','main','ch03_main_ht1080_addback_design_expression','ht1080_addback_design_expression','Five panels: derivation, Western blots, raw construct counts, plasmid-population PCA and residual Y-expression heatmap.'),
    ('F02','main','ch03_main_ht1080_addback_signature_recovery','ht1080_addback_signature_recovery','Two panels: 150-gene signature scores and percentage of the vector-to-sgCTRL score gap recovered.'),
    ('F03','main','ch03_main_ht1080_addback_by_atac_class','ht1080_addback_by_atac_class','Two panels: expression differences by promoter class and adjusted estimates, with vector comparisons first.'),
    ('S01','supp','ch03_supp_ht1080_addback_pca_populations','ht1080_addback_pca_populations','All six established populations, three collection occasions each. PDF re-exported from its saved SVG to correct Arial substitution to embedded Open Sans; plot geometry/data unchanged.'),
    ('S02','supp','ch03_supp_ht1080_addback_baseline_reversal','ht1080_addback_baseline_reversal','Four panels: 150-gene heatmap and fry tests (a,b), plus frozen-100-gene scores and contrasts (c,d).'),
    ('S03','supp','ch03_supp_ht1080_addback_baseline_reversal_frozen100','ht1080_addback_baseline_reversal_frozen100','Original four-panel frozen-100-gene figure. Its a,b are also present in S02 as c,d; retain this complete figure as requested.'),
    ('S04','supp','ch03_supp_ht1080_addback_atac_overlap','ht1080_addback_atac_overlap','Expression responses of the 39 genes with baseline RNA and promoter-accessibility loss.'),
]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def encode_like(text, raw):
    return text.replace('\n','\r\n' if b'\r\n' in raw else '\n').encode('utf-8')

md = ROOT / 'Results_edit.md'
cp = ROOT / 'Chapter3_figure_captions.md'
md_raw, cp_raw = md.read_bytes(), cp.read_bytes()
result, central = md_raw.decode('utf-8').replace('\r\n','\n'), cp_raw.decode('utf-8').replace('\r\n','\n')
index = SELECTED / 'selected_figure_index.tsv'
index_raw = index.read_bytes()
index_lines = index_raw.decode('utf-8').splitlines(keepends=True)
index_header = index_lines[0].strip().split('\t')
old_rows = list(csv.DictReader(io.StringIO(index_raw.decode('utf-8')),delimiter='\t'))
selection = []
new_rows = []
parts = re.split(r'(?=^## )',central,flags=re.M)

for suffix,role,stem,slug,description in ENTRIES:
    working_id = f'CH3-S3.3.5-{suffix}'
    folder = SELECTED / f'CH3_S3_3_5_{suffix}_{slug}'
    folder.mkdir(exist_ok=True)
    for ext in ['pdf','svg','png']:
        assert (FIG/f'{stem}.{ext}').is_file()
    # Preserve figure wording; replace only the provisional identifier.
    pattern = rf'\*\*Figure (?:{re.escape(stem)}|{re.escape(working_id)}) \({role}\)\. (.*?)\*\*'
    match = re.search(pattern,result)
    assert match,stem
    title = match.group(1)
    result = result[:match.start()] + f'**Figure {working_id} ({role}). {title}**' + result[match.end():]
    hashes = []
    for ext in ['pdf','svg','png']:
        src, dst = FIG/f'{stem}.{ext}', folder/f'{stem}.{ext}'
        shutil.copy2(src,dst)
        assert sha(src)==sha(dst)
        hashes.append({'file':dst.name,'source':src.relative_to(ROOT).as_posix(),'sha256':sha(dst)})
    selected_pdf = (folder/f'{stem}.pdf').relative_to(ROOT).as_posix()
    selected_png_abs = (folder/f'{stem}.png').as_posix()
    # The document may point at a hash-named preview; redirect that exact block.
    image_pattern = rf'^\[!\[(?:{re.escape(stem)}|{re.escape(working_id)})\]\([^\n]+$'
    result,n = re.subn(image_pattern,f'[![{working_id}]({selected_png_abs})]({(folder/f"{stem}.pdf").as_posix()})',result,flags=re.M)
    assert n==1,(stem,n)
    old_pdf = f'figures/chapter3/{stem}.pdf'
    # Text citations receive the same working IDs, retaining panel qualifiers.
    def citation(match):
        label = match.group(1)
        panel = re.search(r',\s*panels?\s+.+',label)
        return f'[{"Figure "+working_id+(panel.group(0) if panel else "")}]({selected_pdf})'
    result = re.sub(r'\[([^\]\n]+)\]\('+re.escape(old_pdf)+r'\)',citation,result)
    result = result.replace(']('+old_pdf+')',']('+selected_pdf+')')

    found = False
    for i,part in enumerate(parts):
        if f'**Figure:** `{stem}`' not in part:
            continue
        assert not found
        found=True
        metadata = re.search(r'^\*\*Figure:\*\*.*$',part,re.M).group(0)
        if '**Working ID:**' not in metadata:
            updated = metadata.replace(f'**Figure:** `{stem}`',f'**Figure:** `{stem}` | **Working ID:** {working_id} | **Results:** 3.3.5')
            part=part.replace(metadata,updated,1)
        else:
            assert working_id in metadata
        selected_links=' | '.join(f'[Selected {ext.upper()}]({(folder/f"{stem}.{ext}").relative_to(ROOT).as_posix()})' for ext in ['pdf','svg','png'])
        if '[Selected PDF]' not in part:
            location=part.index('\n\n[PDF]')
            part=part[:location]+'\n\n'+selected_links+part[location:]
        parts[i]=part
    assert found,stem
    note = f'''# {working_id} ({role})

Selected from the current Results_edit.md add-back section on 8 October 2026.

{description}

- [Current Results](../../../../Results_edit.md)
- [Authoritative central caption](../../../../Chapter3_figure_captions.md) (search {working_id})
- [Working source/reproduction notes](../../{stem}_notes.md)
- Working source stem: `{stem}`. Selected PDF/SVG/PNG are byte-identical to the current corresponding source exports.
- `selection_manifest.json` records the copied files and SHA-256 hashes. Numerical source data and renderers remain in the parent chapter directory.

These identifiers are provisional. Main/supplementary designations and figure content have been preserved; the other chapter selections were not changed by this operation.
'''
    (folder/'selection_notes.md').write_text(note,encoding='utf-8',newline='\n')
    (folder/'selection_manifest.json').write_text(json.dumps({'working_id':working_id,'section':'3.3.5','role':role,'selected_date':'2026-10-08','source_stem':stem,'files':hashes},indent=2),encoding='utf-8')
    row={'working_id':working_id,'chapter':'3','section':'3.3.5','role':role,'content':title.rstrip('.'),'status':'selected; current Results_edit.md, 2026-10-08','selected_pdf':(folder/f'{stem}.pdf').relative_to(SELECTED).as_posix(),'source_stem':stem,'sha256_pdf':sha(folder/f'{stem}.pdf')}
    existing=[r for r in old_rows if r['working_id']==working_id or r['source_stem']==stem]
    assert len(existing)<=1,(stem,existing)
    if existing:
        assert existing[0]==row,'A differing selection already exists; archive/review it first'
    else:
        new_rows.append(row)
    selection.append({'working_id':working_id,'role':role,'source_stem':stem,'selected_pdf':row['selected_pdf'],'title':title,'hashes':hashes})

assert md.read_bytes()==md_raw and cp.read_bytes()==cp_raw and index.read_bytes()==index_raw,'Concurrent document edit detected'
md.write_bytes(encode_like(result,md_raw))
cp.write_bytes(encode_like(''.join(parts),cp_raw))
out=io.StringIO()
writer=csv.DictWriter(out,fieldnames=index_header,delimiter='\t',lineterminator='\r\n')
writer.writerows(new_rows)
index.write_bytes(index_raw+out.getvalue().encode('utf-8'))
(ROOT/'.tmp/addback_selected_records.json').write_text(json.dumps(selection,indent=2),encoding='utf-8')
print(f'Selected {len(selection)} figures ({len(new_rows)} new index rows), copied 21 artwork files and updated Results/caption links.')
