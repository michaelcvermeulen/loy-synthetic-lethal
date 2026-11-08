"""Register the reviewed-layout A549 draft figures without selecting them."""
from pathlib import Path
import csv, io, json, shutil, datetime
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'figures/chapter3'
stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
archive=ROOT/'archive'/('a549_rna_figure_drafts_'+stamp)
archive.mkdir(parents=True,exist_ok=True)
entries=[
('main_a549_rna_global_structure','Global transcriptional variation in A549 cells',
'**a,** PCA of 22,792 expressed autosomal and non-PAR X-linked genes in four sgCTRL and four Y-depleted libraries. Labels identify samples within each condition. **b,** Unpaired voom–limma expression differences; the vertical axis shows nominal P values. No gene met FDR <0.05 after Benjamini–Hochberg correction. Blue denotes sgCTRL and vermilion Y-depleted cells in a.'),
('main_a549_rna_sex_chromosomes','Y-linked and pseudoautosomal expression in Y-depleted A549 cells',
'**a,** Expression of eight Y-linked genes. **b,** Selected PAR1 genes, relative to mean sgCTRL expression. Annotations give the change in mean CPM and adjusted P from unpaired voom–limma tests across 31 PAR features; none met genome-wide correction. **c,** Mean log2CPM differences for the X-linked counterparts. NLGN4X abundance was below 0.04 CPM in both groups. Points in a and b show four samples per condition; horizontal marks show means of the displayed values. Blue circles denote sgCTRL and vermilion triangles Y-depleted cells. The dashed line in b marks 100% of mean sgCTRL expression.'),
('main_a549_rna_pathways','Pathway enrichment in A549 is sensitive to inter-gene correlation',
'**a,** The eight Hallmark sets meeting global adjusted P <0.05 in fgsea, ranked using mean log2CPM differences. Negative scores indicate enrichment towards lower expression in Y-depleted cells. **b,** Benjamini–Hochberg-adjusted P values across 23,457 MSigDB sets for fgsea and CAMERA using fixed (0.01) or estimated nonnegative inter-gene correlations. Bold values meet adjusted P <0.05. Cell shading represents −log10 adjusted P on a fixed 0–8 scale. MYC, E2F and G2M fgsea nominal P values reached the analysis floor of 10⁻¹⁰.'),
('supp_a549_rna_qc','A549 RNA-seq quality and parental context',
'**a,** Assigned reads as a percentage of total reads. **b,** Mitochondrial reads as a percentage of assigned reads. **c,** PCA including three parental samples, with normalisation performed separately from the eight-library primary comparison. **d,** Association between mitochondrial read fraction and PC1 in the primary comparison. Labels in c and d identify samples within each condition. Blue denotes sgCTRL, vermilion Y-depleted and grey parental samples.'),
('supp_a549_rna_hallmark_complete','Complete Hallmark enrichment and sensitivity analysis in A549',
'**a,** Normalised enrichment scores for all 50 Hallmark sets, ranked using mean log2CPM differences. Blue indicates negative and orange positive enrichment. **b,** Global Benjamini–Hochberg-adjusted P values across 23,457 MSigDB sets for fgsea and CAMERA with fixed (0.01) or estimated nonnegative inter-gene correlations. Bold values meet adjusted P <0.05; shading represents −log10 adjusted P on a fixed 0–8 scale. MYC, E2F and G2M fgsea nominal P values reached the analysis floor of 10⁻¹⁰.'),
('supp_a549_rna_par_expression','Pseudoautosomal expression across A549 samples',
'Expression of 13 displayed PAR genes relative to each gene’s mean sgCTRL CPM. SPRY3 and VAMP7 are in PAR2; the other displayed genes are in PAR1. Blue circles denote sgCTRL and vermilion triangles Y-depleted samples (four per condition); horizontal marks show means. Dashed lines mark 100% of mean sgCTRL expression. Annotations give Benjamini–Hochberg-adjusted P values from unpaired voom–limma tests across 31 PAR features. Vertical scales vary among genes; measurements do not distinguish X- and Y-derived transcripts.')]
captions=ROOT/'Chapter3_figure_captions.md'
shutil.copy2(captions,archive/captions.name)
text=captions.read_text(encoding='utf-8')
for short,title,caption in entries:
 stem='ch03_'+short
 if f'`{stem}`' in text: raise RuntimeError('Caption already registered: '+stem)
 role=short.split('_')[0]
 text+=f'\n\n## {role} - {title}\n\n**Figure:** `{stem}` | **Section:** 3.4.6.2 | **Status:** Draft for author review\n\n[PDF](figures/chapter3/{stem}.pdf) | [Editable SVG](figures/chapter3/{stem}.svg) | [PNG](figures/chapter3/{stem}.png)\n\n**{title}.** {caption}\n\n[Source and reproduction notes](figures/chapter3/a549_rna_thesis_notes.md)\n'
captions.write_text(text,encoding='utf-8')
results=ROOT/'Chapter3_results.md';shutil.copy2(results,archive/results.name)
text=results.read_text(encoding='utf-8')
oldlines=[line for line in text.splitlines() if 'displays to be assembled' in line and 'A549' in line]
assert len(oldlines)==3
for old,(main,supp) in zip(oldlines,[(entries[0],entries[3]),(entries[1],entries[5]),(entries[2],entries[4])]):
 stem='ch03_'+main[0];ss='ch03_'+supp[0]
 anchor=main[0].split('_')[0]+'---'+main[1].lower().replace(' ','-')
 new=f'[![{main[1]}](figures/chapter3/{stem}.png)](figures/chapter3/{stem}.pdf)\n\n**Figure (main; draft).** {main[1]}. [Full caption](Chapter3_figure_captions.md#{anchor}). [Supplementary figure](figures/chapter3/{ss}.pdf).'
 text=text.replace(old,new,1)
results.write_text(text,encoding='utf-8')
manifest=OUT/'figure_manifest.tsv';shutil.copy2(manifest,archive/manifest.name)
with manifest.open(encoding='utf-8-sig',newline='') as f:
 reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
provenance=json.loads((OUT/'source_data/a549_rna_thesis/provenance.json').read_text(encoding='utf-8'))
names={x['file']:x['sha256'] for x in provenance['exports']}
rows=[r for r in rows if r['file'] not in names]
for short,title,_ in entries:
 for ext in ['pdf','svg','png']:
  name='ch03_'+short+'.'+ext
  rows.append(dict(file=name,role=short.split('_')[0],title=title,guide_set='',status='draft; layout reviewed; awaiting author review',original_source='figures/chapter3/code/make_a549_rna_thesis.py',sha256=names[name]))
with manifest.open('w',encoding='utf-8',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
(OUT/'a549_rna_thesis_notes.md').write_text('''# A549 RNA thesis figure drafts

Source: `RNAseq/analysis/A549/SEPT_18_2026_REANALYZE/outputs/tables`, relative to the projects root. Original analyses remain unchanged. Copies of all input tables and SHA-256 hashes are in [source_data/a549_rna_thesis](source_data/a549_rna_thesis/). Reproduce with `python figures/chapter3/code/make_a549_rna_thesis.py` from the thesis root.

The six figures are drafts for author review, not selected figures. Active captions are held only in [Chapter3_figure_captions.md](../../Chapter3_figure_captions.md). The main figures and supplementary companions are linked in section 3.4.6.2.

Open Sans is loaded from the thesis-figure-style skill, with 8-pt labels at 7.1-inch width (6.8 pt for the full 50-pathway display). Outputs comprise PDF, editable SVG and 300-dpi PNG. No statistical models were refitted. Adjusted P values, expression estimates and PCA coordinates come from frozen source tables. The pathway heatmaps use the same fixed scale across main and supplementary outputs; exact adjusted P values are displayed.

Sampling remains four libraries per core condition with pairing and independent population derivation unresolved in the source manifest. Context PCA includes three parental samples and uses separate normalisation. These do not establish independent editing replicates. The volcano uses model log fold changes; the X-counterpart plot and enrichment ranking use descriptive mean logCPM differences. PAR percentages use arithmetic mean CPM, and the Y expression panel shows log2(CPM+1); these measures are not interchangeable.

The PC1 association panel reports only the retained Spearman coefficient; no P value was retained in the source table. Context PCA axes omit variance percentages because these were not included in the exported coordinate table. The source plots report approximately 39.0% and 19.6%, but these were not substituted as precise values. Expression panels use condition-specific shapes as well as colour. Pathway direction is also encoded by position relative to zero; adjusted P values remain readable without colour.
''',encoding='utf-8')
print('Registered six draft figures, captions, source notes and section links. Backup:',archive)
