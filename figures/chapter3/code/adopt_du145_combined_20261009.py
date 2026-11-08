"""Adopt approved DU145 figure merger with targeted, byte-preserving text edits."""
from pathlib import Path
import re, json, hashlib, shutil
ROOT=Path(__file__).resolve().parents[3]
FIG=ROOT/'figures/chapter3'
ARC=ROOT/'archive/DU145_F01_F03_combined_20261009'
SRC=FIG/'source_data/du145_global_expression_combined'
main='CH3_S3_4_6_3_F01_DU145_RNA_global'
old='CH3_S3_4_6_3_F03_DU145_RNA_expression'
supp='CH3_S3_4_6_3_S02_DU145_RNA_copy_number'
def rel(folder):return f'figures/chapter3/selected/{folder}/{folder}.pdf'
def read(p):return p.read_bytes().decode('utf-8')
def save(p,text):
    dest=ARC/'pre_adoption'/p.relative_to(ROOT)
    if not dest.exists():dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
    original=read(p)
    if '\r\n' in original:text=text.replace('\r\n','\n').replace('\n','\r\n')
    p.write_bytes(text.encode('utf-8'))
def once(text,old,new):
    assert text.count(old)==1,(old[:90],text.count(old))
    return text.replace(old,new,1)

title='Transcriptional variation and expression differences in DU145 cells.'
caption='''**a,** PCA of 18,767 expressed autosomal and non-PAR X-linked genes. Shapes indicate passage; dashed lines connect matched samples. PC1 and PC2 explain 42.3% and 18.4% of variance; R² values describe associations with passage or condition. **b,** Passage-blocked voom–limma analysis (~ passage + condition; residual df = 2), with Benjamini–Hochberg correction across these genes. Red marks FDR <0.05, |log₂ fold change| ≥0.5 and matching directions across all three passage pairs; blue, other significant copy-number-concordant genes; orange, other significant genes; grey, nonsignificant genes. Dashed lines mark the effect threshold. **c,** Minimum and maximum log₂ CPM differences across passage pairs P12–P14. Shading indicates consistent direction; red marks the highlighted subset, orange other FDR <0.05 genes and grey nonsignificant genes. **d,** Mean RNA versus median lpWGS log₂ differences for autosomal/X arms with ≥20 mapped genes. Point area reflects gene count; inset values give descriptive Pearson correlations at gene, 10-Mb-window and arm resolution. **e,** Copy-number classes among significant genes and the highlighted subset. Stable denotes an absolute lpWGS segment difference <0.2; concordant/discordant changes meet this threshold and match/oppose RNA direction. **f,** GRHL2 CPM, with horizontal marks showing means. Blue and red in a and f denote sgCTRL and Y-depleted cells. Samples represent three passages of one established population per condition. Differences are Y-depleted minus sgCTRL; lpWGS uses the separately sampled LOY-2X population. DNA summaries use saved QDNAseq segments and CHM13v2.0 coordinates.'''
supp_title='DU145 copy-number context, model sensitivity and sample structure.'
addition='''**f,** PC1 scores across passages P12–P14, with lines connecting samples within each condition; R² treats passage as categorical. **g,** Male-specific Y-linked reads as a percentage of assigned reads, with the ratio of condition means annotated. **h,** Pearson correlations of TMM-normalised log₂ CPM (r = 0.973–0.988), with average-linkage clustering based on 1 − r. C, sgCTRL; Y, Y-depleted; numbers indicate passage. Blue and red in f–h denote sgCTRL and Y-depleted cells. RNA samples represent three passages of one established population per condition.'''
previews={}
for folder in [main,supp]:
    png=FIG/'selected'/folder/(folder+'.png')
    name=folder+'_combined_'+hashlib.sha256(png.read_bytes()).hexdigest()[:12]+'.png'
    shutil.copy2(png,png.with_name(name));previews[folder]=name

p=ROOT/'Results_edit.md'; original=read(p); text=original.replace('\r\n','\n')
start=text.index('#### 3.3.4.3 ');end=text.index('#### 3.3.4.4 ',start)
section=text[start:end]
def block(code):
    pattern=r'\*\*Figure '+re.escape(code)+r' \(main\)\.[\s\S]*?(?=\n---)'
    found=list(re.finditer(pattern,section));assert len(found)==1
    return found[0]
f03=block('CH3-S3.3.6.3-F03');section=section[:f03.start()]+section[f03.end():]
section=section.replace('\n---\n\n\n---\n','\n---\n')
f01=block('CH3-S3.3.6.3-F01')
abs_dir=(FIG/'selected'/main).as_posix()
main_block=f'**Figure CH3-S3.3.6.3-F01 (main). {title}**\n\n[![CH3-S3.3.6.3-F01]({abs_dir}/{previews[main]})]({abs_dir}/{main}.pdf)\n\n{caption}\n'
section=section[:f01.start()]+main_block+section[f01.end():]
# Rewrite only the link labels to the moved panels, together with their destinations.
section=once(section,f'[Figure CH3-S3.3.6.3-F01, panels a and b, top]({rel(main)})',
             f'[Figure CH3-S3.3.6.3-F01, panel a]({rel(main)}); [Figure CH3-S3.3.6.3-S02, panel f]({rel(supp)})')
section=once(section,f'[Figure CH3-S3.3.6.3-F01, panel c]({rel(main)})',f'[Figure CH3-S3.3.6.3-S02, panel h]({rel(supp)})')
section=once(section,f'[Figure CH3-S3.3.6.3-F01, panel b, bottom]({rel(main)})',f'[Figure CH3-S3.3.6.3-S02, panel g]({rel(supp)})')
for oldpanel,newpanel in [('panels a and b','panels b and c'),('panel c','panel d'),('panel d','panel e'),('panel e','panel f')]:
    a=f'[Figure CH3-S3.3.6.3-F03, {oldpanel}]({rel(old)})'
    assert a in section,a
    section=section.replace(a,f'[Figure CH3-S3.3.6.3-F01, {newpanel}]({rel(main)})')
section=once(section,'**Figure CH3-S3.3.6.3-S02 (supp). DU145 copy-number context and model sensitivity.**',f'**Figure CH3-S3.3.6.3-S02 (supp). {supp_title}**')
section=once(section,f'{supp}/{supp}.png',f'{supp}/{previews[supp]}')
tail='A sensitivity refit restricted to 16,719 genes'
section=once(section,tail,addition+'\n\n'+tail)
text=text[:start]+section+text[end:]
assert 'CH3-S3.3.6.3-F03' not in section
save(p,text)

p=ROOT/'Chapter3_figure_captions.md';text=read(p).replace('\r\n','\n')
def central_span(stem):
    k=text.index('**Figure:** `'+stem+'`');begin=text.rfind('\n## ',0,k)+1
    end=text.find('\n## ',k);return begin,len(text) if end<0 else end
s,e=central_span('ch03_main_du145_rna_global_structure');part=text[s:e]
part=part.replace('## main - DU145 RNA-seq design and global expression structure','## main - DU145 transcriptional variation and expression differences')
cs=part.index('**Global transcriptional variation');ce=part.index('[Source and reproduction notes]',cs)
part=part[:cs]+'**'+title+'** '+caption+'\n\n'+part[ce:]
part=part.replace(f'{main}/{main}.png',f'{main}/{previews[main]}').replace('source_data/du145_section_revision/README.md','source_data/du145_global_expression_combined/README.md')
text=text[:s]+part+text[e:]
s,e=central_span('ch03_main_du145_rna_expression_changes');part=text[s:e]
part=part.replace('**Status:** Selected','**Status:** Retired — merged into F01, panels b–f (9 October 2026)')
part=part.replace('## main - DU145 autosomal expression differences and their copy-number context','## main - DU145 autosomal expression differences and their copy-number context (retired standalone figure)')
text=text[:s]+part+text[e:]
s,e=central_span('ch03_supp_du145_rna_copy_number');part=text[s:e]
part=part.replace('**DU145 copy-number context and model sensitivity.**','**'+supp_title+'**')
part=part.replace('## supp - DU145 copy-number context and model sensitivity','## supp - DU145 copy-number context, model sensitivity and sample structure')
part=once(part,tail,addition+'\n\n'+tail)
part=part.replace(f'{supp}/{supp}.png',f'{supp}/{previews[supp]}').replace('source_data/du145_section_revision/README.md','source_data/du145_global_expression_combined/README.md')
text=text[:s]+part+text[e:];save(p,text)

p=FIG/'selected/selected_figure_index.tsv';text=read(p).replace('\r\n','\n');lines=text.splitlines()
header=lines[0].split('\t')
for i,line in enumerate(lines[1:],1):
    vals=line.split('\t'); row=dict(zip(header,vals)); code=row['working_id']
    if code not in ['CH3-S3.4.6.3-F01','CH3-S3.4.6.3-F03','CH3-S3.4.6.3-S02']:continue
    if code.endswith('F03'):row['status']='retired_merged_into_F01'
    else:
        folder=main if code.endswith('F01') else supp
        row['content']=(title if code.endswith('F01') else supp_title).rstrip('.')
        row['sha256_pdf']=hashlib.sha256((FIG/'selected'/folder/(folder+'.pdf')).read_bytes()).hexdigest()
    lines[i]='\t'.join(row[h] for h in header)
save(p,'\n'.join(lines)+'\n')

for name in ['make_du145_rna_global_structure.py','make_du145_rna_expression_changes.py']:
    p=FIG/'code'/name
    save(p,'"""Compatibility entry point for the approved combined DU145 main figure."""\nfrom make_du145_global_expression_combined import draw_main\n\nif __name__ == "__main__":\n    draw_main()\n')
p=FIG/'code/make_du145_rna_supp.py';text=read(p).replace('\r\n','\n')
text=text[:text.index('# ================= Supplementary 2:')]+'''# S02 now includes the former F01 supporting panels; retain the current combined layout.
from make_du145_global_expression_combined import draw_supp
draw_supp()
''';save(p,text)

notes='''# DU145 F01/F03 merge and S02 extension — 9 October 2026

Author-approved layout change only. The six-panel main figure retains F01a as a and F03a–e as b–f. S02a–e retain their original physical dimensions; f, g and h add the former F01 PC1-by-passage, MSY-read and sample-correlation panels. No panels or statistical results were discarded. F03 exports and its retired central caption remain historical records.

Reproduce: `python figures/chapter3/code/make_du145_global_expression_combined.py` from the thesis root. `--only main` and `--only supp` are available. Legacy entry points use the combined renderer. All 15 numerical input files are frozen in `inputs/` with original paths and hashes in `inputs.json`. Display tables and export hashes are in `main/` and `supp/`. Original scripts, figures and source tables are archived in `archive/DU145_F01_F03_combined_20261009`.

Current captions: [Chapter3_figure_captions.md](../../../../Chapter3_figure_captions.md). Results uses provisional 3.3 figure IDs while selected filenames and central metadata retain 3.4. Preview filenames include content hashes to avoid stale images.
'''
(SRC/'README.md').write_text(notes,encoding='utf-8')
for folder in [main,supp,old]:
    p=FIG/'selected'/folder/'selection_notes.md'
    if folder==old:note='Retired standalone figure, 9 October 2026. Panels a–e are retained as F01b–f. Historical exports are preserved.'
    else:note='Updated 9 October 2026 after the approved F01/F03 merge. '+('Six-panel main figure.' if folder==main else 'Eight-panel supplement; original a–e retained, with former F01 support added as f–h.')
    save(p,note+'\n\nSee [source and reproduction notes](../../source_data/du145_global_expression_combined/README.md) and [central captions](../../../../Chapter3_figure_captions.md).\n')
p=ROOT/'AGENTS.md';text=read(p).replace('\r\n','\n')
text+='''

# DU145 F01/F03 merge

Author approval, 9 October 2026: combine DU145 F01 and F03 under F01, with a, PCA; b, volcano; c, passage consistency; d, arm-level DNA/RNA; e, copy-number classes; f, GRHL2 CPM. Retain all six agreed panels and remove the standalone F03 block from Results_edit.md. Add former F01b top, b bottom and c to existing S02 as f (PC1 by passage), g (MSY reads) and h (sample correlations), retaining S02a–e. Use make_du145_global_expression_combined.py and source_data/du145_global_expression_combined. Legacy renderers delegate to the current layout. Captions, selected index and active Results citations are updated; F03 exports and caption remain historical. No model, numerical input or statistical family changed. Review remaining DU145 recommendations one at a time; do not implement them without further author direction. A549 parental S03 and ssGSEA S06 were subsequently retained unchanged by the author.
''';save(p,text)
(ARC/'adoption.json').write_text(json.dumps({'previews':previews,'main_panels':'a:F01a; b-f:F03a-e','supp_panels':'a-e:unchanged; f:F01b top; g:F01b bottom; h:F01c'},indent=2),encoding='utf-8')
print('Updated Results, central captions, index, notes and legacy entry points.')
