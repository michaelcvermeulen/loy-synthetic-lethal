"""One-time adoption of the four approved DU145 Metascape supplements."""
from pathlib import Path
from datetime import datetime
import hashlib
import json
import re
import shutil

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
CODE = FIG / 'code'
SOURCE = FIG / 'source_data'
DATE = '2026-10-05'
SPECS = [
    ('S06', 'metascape_up_network', 'DU145_Metascape_UP', 'Pathway associations among genes with higher expression in Y-depleted DU145 cells.', 'DU145 Metascape network for genes higher in Y-depleted cells'),
    ('S07', 'metascape_down_network', 'DU145_Metascape_DOWN', 'Pathway associations among genes with lower expression in Y-depleted DU145 cells.', 'DU145 Metascape network for genes lower in Y-depleted cells'),
    ('S08', 'metascape_terms_and_gene_modules', 'DU145_Metascape_terms_genes', 'Enriched pathways and recurring contributing genes in DU145.', 'DU145 Metascape terms and recurrent contributing genes'),
    ('S09', 'metascape_ppi_modules', 'DU145_Metascape_PPI', 'Protein-interaction modules and selected pairs among differentially expressed genes in DU145.', 'DU145 Metascape protein-interaction modules and selected pairs'),
]
MAP = {}
for key, suffix, folder_suffix, title, heading in SPECS:
    folder = f'CH3_S3_4_6_3_{key}_{folder_suffix}'
    MAP[key] = dict(code=f'CH3-S3.4.6.3-{key}', stem=f'ch03_supp_du145_{suffix}',
                    folder=folder, title=title, heading='supp - ' + heading, role='supp',
                    path=f'figures/chapter3/selected/{folder}/{folder}')


def read(path):
    return path.read_bytes().decode('utf-8')


def write(path, text):
    # Preserve each existing file's newline convention.
    old = path.read_bytes() if path.exists() else b''
    nl = '\r\n' if b'\r\n' in old else '\n'
    text = text.replace('\r\n', '\n').replace('\r', '\n')
    path.write_bytes(text.replace('\n', nl).encode('utf-8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cite(key, panels=''):
    spec = MAP[key]
    return f'[{"Figure " + spec["code"] + panels}]({spec["path"]}.pdf)'


def anchor(heading):
    return re.sub(r'[^a-z0-9 -]', '', heading.lower()).replace(' ', '-')


def main():
    index = FIG / 'selected/selected_figure_index.tsv'
    assert all(v['code'] not in read(index) for v in MAP.values()), 'Already adopted'
    paths = [ROOT/n for n in ['Results_edit.md', 'Chapter3_results.md', 'Chapter3_figure_captions.md', 'Chapter3_methods.md']]
    paths += [index, FIG/'figure_manifest.tsv', CODE/'du145_section_figure_map.json',
              CODE/'make_du145_metascape_figures.py', CODE/'make_du145_metascape_ppi_modules.py',
              SOURCE/'du145_metascape_current/README.md', SOURCE/'du145_metascape_current/figure_provenance.json',
              SOURCE/'du145_metascape_terms_and_gene_modules/provenance.json',
              SOURCE/'du145_metascape_ppi_modules/README.md', SOURCE/'du145_metascape_ppi_modules/provenance.json']
    archive = ROOT / 'archive' / ('du145_metascape_supp_adoption_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
    for path in paths:
        dest = archive/path.relative_to(ROOT)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
    # Retain the approved figure bytes; only the two working network stems change role.
    for direction in ['up', 'down']:
        old = f'ch03_main_du145_metascape_{direction}_network'
        new = f'ch03_supp_du145_metascape_{direction}_network'
        for ext in ['pdf', 'svg', 'png']:
            shutil.copy2(FIG/f'{old}.{ext}', FIG/f'{new}.{ext}')
        shutil.copy2(FIG/'previews'/f'{old}_preview.png', FIG/'previews'/f'{new}_preview.png')
    for spec in MAP.values():
        folder = FIG/'selected'/spec['folder']
        folder.mkdir()
        for ext in ['pdf', 'svg', 'png']:
            shutil.copy2(FIG/f'{spec["stem"]}.{ext}', folder/f'{spec["folder"]}.{ext}')
        notes = (f'Selected 5 October 2026 as {spec["code"]} (supp). Working source: `{spec["stem"]}`. '
                 'The approved figure is retained without changing data or layout. '
                 'See [central captions](../../../../Chapter3_figure_captions.md) and '
                 '[adoption notes](../../source_data/du145_metascape_adoption/README.md).\n')
        write(folder/'selection_notes.md', notes)
    # Keep reproduction filenames and status consistent with the selected figures.
    p = CODE/'make_du145_metascape_figures.py'
    s = read(p).replace("stem=f'ch03_main_du145_metascape_{direction}_network'", "stem=f'ch03_supp_du145_metascape_{direction}_network'")
    s = s.replace('draft counterparts of HT1080 F04A/F04B and S04', 'selected as CH3-S3.4.6.3-S06, S07 and S08 on 2026-10-05')
    write(p, s)
    p = CODE/'make_du145_metascape_ppi_modules.py'
    write(p, read(p).replace('Draft for review; not selected', 'Selected as CH3-S3.4.6.3-S09 on 2026-10-05'))
    # Update the four central caption entries, preserving their scientific content.
    p = ROOT/'Chapter3_figure_captions.md'
    caps = read(p).replace('\r\n', '\n')
    bodies = {}
    for key, spec in MAP.items():
        old_heading = spec['heading'].replace('supp - ', 'main - ', 1) if key in ['S06', 'S07'] else spec['heading']
        pattern = r'^## ' + re.escape(old_heading) + r'\n.*?(?=^## |\Z)'
        m = re.search(pattern, caps, re.M | re.S)
        assert m, old_heading
        block = m.group()
        body_start = block.index('**' + spec['title'] + '**')
        body = block[body_start:].strip()
        bodies[key] = body
        intro = f'## {spec["heading"]}\n\n**Working ID:** {spec["code"]} | **Role:** supp | **Status:** Selected | **Results:** 3.4.6.3\n\n'
        intro += ' | '.join(f'[Selected {ext.upper()}]({spec["path"]}.{ext})' for ext in ['pdf', 'svg', 'png'])
        intro += f'\n\n[Working PDF](figures/chapter3/{spec["stem"]}.pdf)\n\n'
        caps = caps[:m.start()] + intro + body + '\n\n' + caps[m.end():]
    write(p, caps)
    # The new account follows individual-gene findings and precedes ranked GSEA.
    paragraphs = [
        'To examine the functional themes represented by these expression changes, we analysed the 31 higher- and 83 lower-expression genes separately using Metascape, with the expressed-gene background supplied for enrichment testing. Higher-expression genes were enriched for the Hallmark epithelial–mesenchymal transition signature (six genes; q ≈ 0.016), collagen biosynthesis and modifying enzymes (four genes; q ≈ 0.020), and Reactome extracellular-matrix organisation (six genes; q ≈ 0.025). The collagen-biosynthesis overlap comprised *COL4A1*, *COL4A2*, *COL12A1* and *P3H2*, while *BMP7* and *ADAM12* additionally contributed to the extracellular-matrix term (' + cite('S06') + ').',
        'Lower-expression genes were enriched for epithelial- and CDH1-associated reference signatures. Nineteen genes overlapped ONDER_CDH1_TARGETS_2_DN (q ≈ 2.5 × 10⁻⁹), and 15 overlapped CHARAFE_BREAST_CANCER_LUMINAL_VS_MESENCHYMAL_UP (q ≈ 2.5 × 10⁻⁵), including *CLDN7*, *ST14*, *AP1M2* and *GRHL2*. Immune- and cytokine-associated annotations were also represented, including positive regulation of interleukin-12 production (four genes; q ≈ 0.005) and cell chemotaxis (six genes; q ≈ 0.020; ' + cite('S07') + '). These overlaps extended the epithelial- and immune-associated patterns seen at individual genes.',
        'The highest-ranked enrichment terms in both directions are shown in ' + cite('S08', ', panels a and b') + '. Several higher-expression genes contributed repeatedly to related terms: *COL4A2*, *COL12A1* and *COL4A1* occurred in 24, 16 and 14 terms, respectively, within the UP network (' + cite('S08', ', panel c') + '). This recurrence shows that several enrichment results reflect overlapping annotations of the same genes.',
        'Protein-interaction analysis provided a complementary view of these gene groups. The higher-expression network contained a collagen-associated module comprising *COL4A1*, *COL4A2*, *COL12A1*, *P3H2* and *TAF9*. Two modules among lower-expression genes contained *CSF1*, *SIRPA*, *SYK* and *B2M*, and *CCR9*, *CCL22* and *CXCL1*, respectively (' + cite('S09', ', panels a–c') + '). Selected pairs outside these modules included *BMP7–NPTX1* and *APLP1–CNTN1* among higher-expression genes, and *CLDN7–TACSTD2* among lower-expression genes (' + cite('S09', ', panels d–f') + '). These links represent database-derived protein associations rather than interactions measured in the DU145 populations.'
    ]
    for name in ['Chapter3_results.md', 'Results_edit.md']:
        p = ROOT/name
        s = read(p).replace('\r\n', '\n')
        start = re.search(r'^#+ 3\.4\.6\.3\b', s, re.M).start()
        end = re.search(r'^#+ 3\.4\.6\.4\b', s, re.M).start()
        section = s[start:end]
        addition = '##### Pathway and protein-interaction analysis of differentially expressed genes\n\n'
        for key, paragraph in zip(MAP, paragraphs):
            spec = MAP[key]
            addition += paragraph + '\n\n'
            if name == 'Results_edit.md':
                absolute = ROOT.as_posix() + '/' + spec['path']
                body = bodies[key].replace('**' + spec['title'] + '** ', '', 1)
                # Keep source links as metadata in the central caption rather than duplicated inline prose.
                body = re.sub(r'\n\n\[Source data[^\n]+\Z', '', body)
                addition += f'---\n\n**Figure {spec["code"]} (supp). {spec["title"]}**\n\n[![{spec["code"]}]({absolute}.png)]({absolute}.pdf)\n\n{body}\n\n---\n\n'
            else:
                addition += f'[![{spec["title"]}]({spec["path"]}.png)]({spec["path"]}.pdf)\n\n**Figure {spec["code"]} (supp). {spec["title"]}** [Full caption](Chapter3_figure_captions.md#{anchor(spec["heading"])}).\n\n'
        marker = '##### Coordinated transcriptional programmes\n\n'
        assert section.count(marker) == 1
        section = section.replace(marker, addition + marker, 1)
        old = 'To examine coordinated expression differences, we used preranked GSEA on mean log2CPM differences'
        new = 'To assess coordinated expression differences across the full ranked gene list, we used preranked GSEA on mean log2CPM differences'
        assert old in section
        section = section.replace(old, new, 1)
        old = 'Reduced epithelial-associated expression was not accompanied by clear evidence of a mesenchymal programme.'
        new = 'Although matrix-associated genes and the Hallmark epithelial–mesenchymal transition signature were enriched among higher-expression genes, changes in canonical mesenchymal markers were limited.'
        assert old in section
        section = section.replace(old, new, 1)
        section = section.replace('The findings therefore support a partial reduction in epithelial-associated expression without establishing an epithelial-to-mesenchymal transition.', 'Together, the findings support reduced epithelial-associated expression accompanied by increased matrix-associated transcription, without establishing a complete epithelial-to-mesenchymal transition.', 1)
        old = 'Histone- and G2/M-associated expression increases provided additional, more selective features of the edited population.'
        new = 'Increased collagen- and extracellular-matrix-associated expression accompanied these reductions, with histone- and G2/M-associated increases providing additional features of the edited population.'
        assert old in section
        section = section.replace(old, new, 1)
        write(p, s[:start] + section + s[end:])
    # Methods needed to interpret the newly adopted analysis.
    p = ROOT/'Chapter3_methods.md'
    s = read(p)
    marker = '#### DU145 gene-set sensitivity screens'
    method = '''#### DU145 Metascape enrichment and protein-interaction networks

The 31 higher- and 83 lower-expression genes meeting FDR <0.05, |log₂ fold change| ≥0.5 and matching directions across all three passage pairs were submitted separately to Metascape with a custom expressed-gene background (18,760 unique symbols representing 18,767 tested Ensembl genes). All 114 input symbols were recognised. The saved reports, using the 1 July 2026 database release, record 16,520 background genes; the enrichment tables use 15,651 genes as the library denominator. The mapped background list was not included in the exports, so individual background mapping losses could not be reconstructed. Enrichment used cumulative hypergeometric tests with Benjamini–Hochberg correction; reported terms passed nominal P <0.01, at least three overlapping genes and an enrichment factor >1.5.

Results and figures use the complete 31- or 83-gene tables (`GO_AllLists.csv`), rather than the smaller evidence-filtered subsets in the final Metascape summaries. Reported q values were reconstructed from rounded exported log₁₀ q values and are approximate. Term-network displays retain all 102 UP and 105 DOWN nodes and 472 and 689 unique connections, respectively. Display communities were identified using Louvain clustering with seed 42. Node area encodes full-list −log₁₀ P; community labels summarise constituent annotations. The recurrence display counts membership of each higher-expression gene across the 102 UP network terms. The top-term panels show 20 terms per direction ranked by nominal P, including terms that did not meet q <0.05.

The protein-interaction figure retains the three exported MCODE modules, comprising 12 proteins and 17 links. Three selected disconnected pairs from the complete networks (BMP7–NPTX1, APLP1–CNTN1 and CLDN7–TACSTD2) add six proteins and three links. Module membership, edges and evidence attributes were retained without rerunning network detection. Protein associations are database-derived and include physical, functional and inferred evidence; they were not measured in the DU145 populations. Source reports, full-list enrichment statistics and plotted network tables are retained with the [figure source data](figures/chapter3/source_data/du145_metascape_current/README.md).

'''
    assert marker in s and '#### DU145 Metascape enrichment' not in s
    write(p, s.replace(marker, method + marker, 1))
    # Selected index and source mapping.
    s = read(index).rstrip('\r\n') + '\n'
    for spec in MAP.values():
        s += '\t'.join([spec['code'], '3', '3.4.6.3', 'supp', spec['title'].rstrip('.'), 'selected',
                        f'{spec["folder"]}/{spec["folder"]}.pdf', spec['stem'], sha(FIG/f'{spec["stem"]}.pdf')]) + '\n'
    write(index, s)
    p = CODE/'du145_section_figure_map.json'
    mapping = json.loads(read(p)); mapping.update(MAP)
    write(p, json.dumps(mapping, ensure_ascii=False, indent=2) + '\n')
    # Manifest keeps old main-role working exports as superseded candidates.
    p = FIG/'figure_manifest.tsv'; lines = read(p).splitlines(); header = lines[0].split('\t')
    entries = [dict(zip(header, row.split('\t'))) for row in lines[1:]]
    for entry in entries:
        if entry['file'].startswith(('ch03_main_du145_metascape_up_network.', 'ch03_main_du145_metascape_down_network.')):
            entry['status'] = 'superseded 2026-10-05; selected under supplementary stem'
        for spec in MAP.values():
            if entry['file'].startswith(spec['stem'] + '.'):
                entry.update(role='supp', title=spec['heading'], status=f'selected {DATE} as {spec["code"]}')
    existing = {e['file'] for e in entries}
    for spec in MAP.values():
        for ext in ['pdf', 'svg', 'png']:
            name = f'{spec["stem"]}.{ext}'
            if name not in existing:
                entries.append(dict(file=name, role='supp', title=spec['heading'], guide_set='', status=f'selected {DATE} as {spec["code"]}',
                                    original_source='figures/chapter3/code/make_du145_metascape_figures.py', sha256=sha(FIG/name)))
    write(p, '\t'.join(header) + '\n' + '\n'.join('\t'.join(e.get(k, '') for k in header) for e in entries) + '\n')
    # Source records point to the selected figures and the revised reproduction script.
    for p in [SOURCE/'du145_metascape_current/figure_provenance.json', SOURCE/'du145_metascape_terms_and_gene_modules/provenance.json']:
        data = json.loads(read(p));data['status'] = 'selected as CH3-S3.4.6.3-S06, S07 and S08 on 2026-10-05'
        for out in data['exports']:
            out['file'] = out['file'].replace('ch03_main_du145_metascape_', 'ch03_supp_du145_metascape_')
        data['script_sha256'] = sha(CODE/'make_du145_metascape_figures.py')
        write(p, json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    p = SOURCE/'du145_metascape_ppi_modules/provenance.json'
    data = json.loads(read(p));data['status'] = 'Selected as CH3-S3.4.6.3-S09 on 2026-10-05'
    write(p, json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    p = SOURCE/'du145_metascape_current/README.md';s = read(p)
    s = s.replace('This review has not changed thesis Results prose, captions or selected figures.', 'The initial import review did not change thesis prose. The four figures were subsequently selected and incorporated into the DU145 Results on 5 October 2026; see the adoption notes below.')
    s = s.replace('## Figure drafts generated 5 October 2026', '## Figures selected 5 October 2026')
    s = s.replace('ch03_main_du145_metascape_', 'ch03_supp_du145_metascape_').replace('#main---du145-metascape-network', '#supp---du145-metascape-network')
    s = s.replace('They are drafts and have not been assigned selected figure IDs or inserted into Results prose.', 'The UP network, DOWN network and terms/gene composite are selected as CH3-S3.4.6.3-S06, S07 and S08. The PPI modules and selected pairs are S09. All four are cited in the DU145 Results and displayed with captions in Results_edit.md. The three-page review PDF remains an unnumbered viewing copy.')
    write(p, s)
    p = SOURCE/'du145_metascape_ppi_modules/README.md';s = read(p)
    s = s.replace('This draft combines', 'This figure combines').replace('It is not assigned a selected figure identifier or inserted into Results.', 'It was selected as CH3-S3.4.6.3-S09 on 5 October 2026 and incorporated into the DU145 Results.')
    write(p, s)
    notes = SOURCE/'du145_metascape_adoption'; notes.mkdir(exist_ok=True)
    write(notes/'README.md', '# DU145 Metascape supplementary-figure adoption\n\nSelected 5 October 2026 at the author\'s request. The four approved figures are S06 (UP network), S07 (DOWN network), S08 (terms and recurring genes) and S09 (PPI modules and selected pairs). Figure data and layouts are unchanged. The two network working stems now use `supp`.\n\nBoth Results sources contain the new account before ranked GSEA; Results_edit.md includes all four previews and full captions. Existing independently edited prose was preserved. The later EMT interpretation and final DU145 summary were adjusted to acknowledge the matrix-associated findings. Methods now document the input lists, mapped background, complete-list statistics, recurrence and PPI displays.\n\nActive captions are in [Chapter3_figure_captions.md](../../../../Chapter3_figure_captions.md).\n\nArchive: `' + str(archive.relative_to(ROOT)) + '`.\n')
    write(notes/'adoption.json', json.dumps(dict(date=DATE, archive=str(archive.relative_to(ROOT)), figures=MAP, paragraphs=paragraphs), ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(dict(archive=str(archive), selected=list(MAP), paragraphs_added=len(paragraphs)), indent=2))


if __name__ == '__main__':
    main()
