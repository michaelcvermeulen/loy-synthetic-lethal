"""One-time integration of the refined molecular display and approved numbering."""
from pathlib import Path
import hashlib,json,re,shutil
ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'figures/chapter4/source_data/ht1080_molecular_integration_supp_20261028'
ARCHIVE=Path((Path(__file__).parent/'archive_path.txt').read_text(encoding='utf-8'))
ORIGINAL=json.loads((ARCHIVE/'original_hashes.json').read_text(encoding='utf-8'))
STEM='ch04_supp_ht1080_molecular_integration'
TITLE='Transcriptomic, chromatin-accessibility, and add-back profiles of HT1080 screen candidates and their target programmes.'
CAPTION='''**a,** Pathway-level comparison of the CRISPR screen, baseline molecular profiles and Y-gene add-backs. Screen columns show preranked fgsea of drugZ and D0-adjusted IsogenicZ rankings (−normZ) across 5,880 Hallmark, Reactome, GO:BP and GO:CC sets. Baseline columns show Chapter 3 LOY-versus-CTRL RNA-seq fgsea and trend-corrected promoter ATAC-seq fgsea. Add-back columns show CAMERA tests for UTY, KDM5D and UTY+KDM5D against the empty-vector (EV) control. Rows are grouped as screen sensitizer complexes (SAGA, SWI/SNF), downstream regulated target programmes (sterol, MYC and proliferation), and screen resistor pathways (mitochondrial translation, respiratory chain and cytosolic translation). Cell values give FDRs, and colour intensity encodes signed −log₁₀(FDR). **b,** Gene-level screen scores (−normZ) beside baseline RNA-seq, promoter ATAC-seq and add-back RNA-seq log₂ fold changes for representative candidates and MYC. In both panels, the two screen columns use purple for sensitizers (relative depletion in LOY) and teal for resistors (relative enrichment in LOY). Baseline and add-back columns use red for higher and blue for lower expression/accessibility in LOY versus CTRL or add-back versus EV, respectively. Separate horizontal colourbars beneath each panel define these scales. Bold values indicate FDR <0.05 within the corresponding analysis; grey cells labelled n.t. were not tested. Colour limits are ±8 for signed −log₁₀(FDR) and screen −normZ, and ±1.5 for gene-level log₂ fold change; values beyond these limits retain their printed values. Tests and FDR correction families differ between columns, so colour intensity does not provide a common effect-size scale across analyses.'''


def read_current(rel):
    p=ROOT/rel
    assert hashlib.sha256(p.read_bytes()).hexdigest()==ORIGINAL[rel],f'Changed since archive: {rel}'
    return p.read_text(encoding='utf-8')


def renumber(text):
    def number(supp,n):
        if supp:return True,n+1 if 7<=n<=11 else n
        if n==4:return True,7
        return False,n-1 if 5<=n<=9 else n
    def ref(m):
        prefix=m['prefix'];supp,n=number('Supp' in prefix,int(m['n']))
        name='Supp Fig' if supp else 'Fig'
        if prefix.startswith('**'):name='**'+name+'**'
        return name+' 4.'+str(n)
    text=re.sub(r'(?P<prefix>\*\*(?:Supp )?Fig\*\*|(?:Supp )?Fig) 4\.(?P<n>\d+)(?!\d)',ref,text)
    def anchor(m):
        supp,n=number(bool(m['supp']),int(m['n']))
        return ('supp-' if supp else '')+'fig-4-'+str(n)
    return re.sub(r'(?<![\w-])(?P<supp>supp-)?fig-4-(?P<n>\d+)(?!\d)',anchor,text)


def main():
    fig=ROOT/'figures/chapter4'
    preview=fig/'review'/f'{STEM}_preview.png'
    digest=hashlib.sha256(preview.read_bytes()).hexdigest()[:12]
    frozen=preview.with_name(f'{STEM}_{digest}.png');shutil.copy2(preview,frozen)
    central=read_current('Chapter4_figure_captions.md')
    blocks=re.split(r'(?=^## )',central,flags=re.M)
    updated=[]
    for block in blocks:
        if block.startswith('## main - Fig 4.4.'):
            block=(f'## supp - Supp Fig 4.7. {TITLE}\n\n<a id="supp-fig-4-7"></a>\n\n'
                   f'![HT1080 screen and molecular profiles]({frozen.as_posix()})\n\n'
                   f'**Supp Fig 4.7. {TITLE}** {CAPTION}\n\n'
                   f'[PDF](figures/chapter4/{STEM}.pdf) · [Editable SVG](figures/chapter4/{STEM}.svg) · '
                   f'[600-dpi PNG](figures/chapter4/{STEM}.png) · '
                   '[Source and reproduction notes](research/ht1080_molecular_integration_supp_20261028/README.md)\n\n')
        elif re.match(r'## (main|supp) - (?:Supp )?Fig 4\.\d+\.',block):
            block=renumber(block)
        updated.append(block)
    (ROOT/'Chapter4_figure_captions.md').write_text(''.join(updated),encoding='utf-8')

    results=read_current('draftv2/chapter4/chapter4_results.md')
    # Results has already been renumbered by the author: only replace this figure block.
    marker='**Supp Fig 4.7. '+TITLE+'**'
    assert results.count(marker)==1
    caption_start=results.index('\n\n',results.index(marker))+2
    caption_end=results.index('\n\n---',caption_start)
    links=(f' [Full caption](../../Chapter4_figure_captions.md#supp-fig-4-7) · '
           f'[PDF]({(fig/(STEM+".pdf")).as_posix()}) · [Editable SVG]({(fig/(STEM+".svg")).as_posix()}) · '
           f'[Gene evidence table]({(fig/"source_data/ch4_molecular_integration_20261021/focused_gene_evidence.csv").as_posix()}).')
    results=results[:caption_start]+CAPTION+links+results[caption_end:]
    results=results.replace('A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/figures/chapter4/review/ch04_main_ht1080_molecular_integration_draft_54c568fc3a53.png',frozen.as_posix())
    (ROOT/'draftv2/chapter4/chapter4_results.md').write_text(results,encoding='utf-8')

    # Remaining chapter sources still have the older numbering; preserve all prose.
    for name in ['chapter4_discussion.md','chapter4_methods.md','chapter4_supplementary_methods.md','chapter4_introduction.md']:
        rel='draftv2/chapter4/'+name
        text=read_current(rel);changed=renumber(text)
        if changed!=text:(ROOT/rel).write_text(changed,encoding='utf-8')

    register=read_current('draftv2/chapter4/chapter4_figure_numbering.md')
    rows=[]
    for line in register.splitlines():
        if not re.match(r'\| (?:Supp )?Fig 4\.\d+ ',line):continue
        fields=[s.strip() for s in line.strip('|').split('|')]
        fields[0]=renumber(fields[0])
        if fields[4]=='ch04_main_ht1080_molecular_integration_draft':
            fields[1]='supp';fields[3]=TITLE;fields[4]=STEM
        rows.append(fields)
    rows.sort(key=lambda f:('Supp' in f[0],int(f[0].split('.')[-1])))
    register='''# Chapter 4 figure numbering

Updated 28 October 2026. Use **Fig** 4.1 and **Supp Fig** 4.1 with separate sequences and lowercase panel suffixes. Central caption anchors and all current draft 2 references follow this register; descriptive asset filenames remain provenance keys.

| Number | Role | Section | Caption title | Asset stem |
|---|---|---|---|---|
'''+''.join('| '+' | '.join(row)+' |\n' for row in rows)+'''
The established numbered sequence contains 8 main and 12 supplementary figures. The next unused numbers are Fig 4.9 and Supp Fig 4.13. The provisional public-analysis figures labelled Fig 4.0 and Supp Fig 4.0a–d in Results retain their author-assigned labels and are outside this register.

**Current notes — 28 October 2026**

- The molecular-integration figure is now Supp Fig 4.7. Former main figures 4.5–4.8 are 4.4–4.7, and the Discussion figure moves from 4.9 to 4.8. Former supplements 4.7–4.11 are 4.8–4.12.
- Fig 4.3 contains endpoint candidates (a), endpoint gametologs (b), and SAGA organisation/scores (c). Supp Fig 4.3 contains method comparisons (a,b) and D0 heatmaps (c,d). Reproduce them with `figures/chapter4/code/make_ht1080_support_compact.py`.
- Supp Fig 4.1 includes the chromosome-shift strip; Supp Fig 4.2 combines guide non-detection and overrepresentation; Supp Fig 4.6 combines Hallmark/expanded GSEA and Metascape analyses; Supp Fig 4.9 combines A549 trimmed-P3 and IsogenicZ sensitivity. Superseded standalone figures remain historical.
- Fig 4.6h–i retain genomic PCR and clone-ddPCR panels; Fig 4.7 remains the interim competition time course.
- Current molecular-integration reproduction and validation: [source notes](../../research/ht1080_molecular_integration_supp_20261028/README.md). Earlier figure-split history: [record](../../research/fig4_3_split_20261028/README.md).
'''
    (ROOT/'draftv2/chapter4/chapter4_figure_numbering.md').write_text(register,encoding='utf-8')
    record={'figure_number':'Supp Fig 4.7','stem':STEM,'title':TITLE,'caption':CAPTION,
            'preview':frozen.as_posix(),'archive':ARCHIVE.as_posix(),
            'main_number_map':{'4':'Supp Fig 4.7','5':'4','6':'5','7':'6','8':'7','9':'8'},
            'supp_number_map':{str(n):str(n+1) for n in range(7,12)},'results_prose_preserved':True}
    (DATA/'integration.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('Updated captions, figure register, target Results block and remaining chapter references.')


if __name__=='__main__':main()
