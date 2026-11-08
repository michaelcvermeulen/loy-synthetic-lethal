"""Validate completed outputs and build a self-contained, searchable review report."""
from pathlib import Path
import base64,hashlib,json,html
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
ROOT=Path(__file__).resolve().parent;O=ROOT/'outputs/programme_analysis_v1';T=O/'tables'
def read(s):return pd.read_csv(T/f'{s}.csv')
summary=read('DE_summary');de=read('gene_consistency');g=read('programme_consistency');tf=read('TF_consistency')
allg=read('GSEA_primary_all');cam=read('CAMERA_overall_sensitivity')
integration=read('programme_integration_overview');short=read('provisional_theme_shortlist')
oldint=pd.read_csv(ROOT/'outputs/tables/batchDE_passage_period__MAY_minus_AUG.csv')
newint=read('DE_interaction');parent=read('parental_context_overlap')
# Independent numerical checks of exported inference and eligibility.
for r in summary.itertuples():
 q=read(f'DE_{r.contrast}')
 assert q.gene_id.is_unique and np.isfinite(q[['logFC','PValue','FDR']]).all().all()
 assert np.allclose(false_discovery_control(q.PValue),q.FDR,rtol=1e-10,atol=1e-12)
 assert ((q.FDR<.05)&(q.logFC.abs()>=1)).sum()==r.FDR05_LFC1
assert len(de)==25526 and de.directionally_consistent.sum()==662
assert np.array_equal(de.directionally_consistent,(de.FDR<.05)&(de.logFC.abs()>=1)&(np.sign(de.logFC)==np.sign(de.AUG_logFC))&(np.sign(de.logFC)==np.sign(de.MAY_logFC)))
assert np.allclose(de.MAY_logFC-de.AUG_logFC,de.interaction_logFC)
for k,q in allg.groupby('contrast'):
 assert q.pathway.is_unique and (q.test_status=='tested').all()
 assert np.allclose(false_discovery_control(q.pval),q.FDR,rtol=1e-9,atol=1e-12)
 universe=set(read(f'ranking_{k}').gene_name)
 assert all(set(str(le).split(';'))<=universe for le in q.leadingEdge)
assert np.array_equal(g.directionally_consistent,(g.FDR<.05)&(np.sign(g.NES)==np.sign(g.AUG_2025_NES))&(np.sign(g.NES)==np.sign(g.MAY_2026_NES)))
for r in pd.read_csv(O/'provenance/input_hashes.csv').itertuples():
 assert hashlib.sha256((ROOT/r.path).read_bytes()).hexdigest()==r.sha256,r.path
for r in json.loads((O/'provenance/ATAC_integration.json').read_text())['sources']:
 assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256']
validation={'validated_utc':datetime.now(timezone.utc).isoformat(),'checks':['Independent BH recalculation for every DE/GSEA contrast',
 'Unique gene identifiers and finite test statistics','Exact interaction coefficient arithmetic',
 'Gene and pathway directional-selection rules','Leading-edge membership in tested ranking universe',
 'Original scientific input hashes unchanged','ATAC source hashes unchanged'],
 'genes_consistent':int(de.directionally_consistent.sum()),'pathways_consistent':int(g.directionally_consistent.sum()),
 'pathways_both_periods':int(g.both_periods_FDR05.sum()),'TFs_consistent':int(tf.directionally_consistent.sum()),
 'GSEA_floor_pvalues':int((allg.pval<=1e-10).sum()),'GSEA_missing_error_estimates':int(allg.log2err.isna().sum())}
(O/'provenance/validation.json').write_text(json.dumps(validation,indent=2))
files=['run_programme_analysis_v1.R','build_programme_networks_v1.py','integrate_programmes_ATAC_v1.py',
 'check_programme_camera_v1.R','plot_programme_analysis_v1.py','report_programme_analysis_v1.py','PLAN_20260918.md']
(O/'provenance/code_hashes.json').write_text(json.dumps({f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},indent=2))
def table(d):
 return '<div class="scroll">'+d.to_html(index=False,escape=True,border=0,float_format=lambda v:f'{v:.4g}',classes='results')+'</div>'
def link(name):return f'<a href="tables/{name}.csv">Download {html.escape(name)}</a>'
parts=[]
def section(title,body):parts.append(f'<section><h2>{html.escape(title)}</h2>{body}</section>')
section('1. Question, design and agreed decisions',
 '<p>Which transcriptional programmes change in Y-depleted HT1080 cells? Lead with directionally consistent changes across August and May; investigate period-dependent responses separately.</p>'
 '<p>The main model is limma-voom batch + condition on 12 libraries, with existing core TMM exposures and 25,526 autosomal/non-PAR X genes. Period effects and interaction use batch × condition. The parental comparison is August-only, unpaired, with a separate six-library expression filter and TMM normalisation.</p>'
 '<p>Samples repeatedly measure one established population per condition. Model-based tests assume independent residuals; these data do not establish reproducibility across independently generated Y-depleted populations. Direction agreement does not imply equal magnitude.</p>'
 '<p><a href="../../PLAN_20260918.md">Full agreed analysis plan and implementation choices</a></p>')
section('2. Gene-level results',table(summary)+
 '<p>All 662 genes passing overall FDR &lt;0.05 and |log₂FC| ≥1 have the same fitted direction in both periods: 229 higher and 433 lower. Tests assess a nonzero effect; the fold-change requirement filters estimated effects rather than testing a minimum magnitude.</p>'+link('gene_consistency'))
plot=O/'figures/ch03_supp_ht1080_rna_programme_overview.png'
section('3. Visual overview','<img class="figure" alt="RNA programme overview" src="data:image/png;base64,'+base64.b64encode(plot.read_bytes()).decode()+'">'
 '<p><b>a</b>, Period-specific coefficients; coloured genes pass the overall and directional criteria. <b>b</b>, August sgCTRL minus parental, with the 12 qualifying genes labelled. <b>c</b>, Four lowest-FDR directionally consistent pathways per collection; colours show NES. <b>d</b>, Twelve lowest-FDR directionally consistent TF patterns; colours show target-pattern scores. Stars in c/d indicate FDR &lt;0.05 within the defined contrast-specific testing family. Complete numerical results are linked below.</p>')
camera_summary=cam.assign(significant=cam.FDR<.05).groupby('correlation_mode',as_index=False).significant.sum()
section('4. Programme enrichment and sensitivity',
 '<p>Hallmark, Reactome and GO Biological Process were screened equally using signed moderated t rankings of all tested unique symbols. The core ranking contains 25,461 symbols. Gene sets with 15–500 represented members were eligible. BH correction spans all three selected collections separately for each contrast. The primary screen tested 4,868 sets: 220 passed FDR &lt;0.05, all with matching NES directions across periods; 55 were also significant in each period.</p>'
 '<div class="note"><b>Method sensitivity:</b> the correlation-aware CAMERA check finds 17 pathways at fixed correlation 0.01 and none with estimated nonnegative correlation. The GSEA findings therefore remain exploratory ranked-expression patterns; cross-method support is limited.</div>'
 +table(camera_summary)+f'<p>Numerical accounting across all GSEA contrasts: {validation["GSEA_floor_pvalues"]} P values reached the 1e-10 computational floor; {validation["GSEA_missing_error_estimates"]} rows have unavailable Monte Carlo error estimates. All returned finite P values. Floor values are bounds, not exact probabilities.</p>'
 +link('GSEA_primary_all')+' | '+link('CAMERA_overall_sensitivity')+
 table(g[g.directionally_consistent].sort_values('FDR')[['pathway','analysis_collection','NES','FDR','AUG_2025_NES','AUG_2025_FDR','MAY_2026_NES','MAY_2026_FDR','both_periods_FDR05','period_leading_edge_jaccard']]))
section('5. Pathway groups and RNA–ATAC overview',
 '<p>Connect pathways sharing at least three leading-edge genes and Jaccard ≥0.25; weighted Louvain communities are formed separately by direction. All 220 directionally consistent pathways are retained, including isolated pathways, yielding 88 groups. Representative pathway names are provisional identifiers, not inferred cell identity. Full-membership and edge-threshold sensitivities are exported.</p>'
 '<p>Integration uses existing hg38/GRCh38 promoter-overlap links (TSS ±2 kb), unique exact gene symbols and current RNA effects. It covers 12,834 tested RNA symbols with promoter measurements; one ambiguous ATAC symbol is excluded. RNA/ATAC specimen and time matching remain unresolved: this is a group-effect comparison, not a sample-matched correlation. Minimum promoter-peak FDR is not a gene-level FDR. Opposing peaks are retained.</p>'
 +table(integration)+link('programme_integration_overview')+' | '+link('programme_network_sensitivity'))
section('6. Complementary STRING modules',
 '<p>The new directionally consistent gene set gives 500 exact preferred-symbol STRING mappings. Association score ≥700 and the existing MCODE settings identify five modules; the separate physical-only graph has no MCODE modules. These are reference protein-association patterns, not measured coexpression or binding.</p>'
 +table(pd.read_csv(T/'STRING/coverage.csv'))+table(pd.read_csv(T/'STRING/modules.csv'))+
 '<p>Module GO/Reactome annotation uses one-sided hypergeometric tests and BH across every module × term test, separately for the expressed STRING-mapped and direction-candidate backgrounds. Shared database knowledge limits independence of network discovery and annotation.</p>')
section('7. Inferred regulators',
 '<p>CollecTRI ULM uses model log-fold changes and at least five represented targets. The 60 prioritised factors have overall target-pattern FDR &lt;0.05 and consistent score direction across periods. TF–programme links preserve overlapping targets; at least three shared targets are required for the overview counts. These tests describe conditional target-expression patterns rather than independent population replication.</p>'
 '<p>Existing chromVAR and sequence-motif enrichment are joined through exact or uniquely case-insensitive TF labels; dimer members and motif families remain ambiguous about specific binding. Passage-adjusted chromVAR results are retained as sensitivity. These measurements do not establish TF occupancy or causality.</p>'
 +table(tf[tf.directionally_consistent].sort_values('FDR')[['source','score','FDR','AUG_2025score','AUG_2025FDR','MAY_2026score','MAY_2026FDR']])+link('TF_programme_links')+' | '+link('TF_RNA_chromVAR_evidence')+' | '+link('TF_RNA_sequence_enrichment'))
interaction_compare=newint[['gene_id','gene_name','logFC','FDR']].merge(oldint[['gene_id','logFC','FDR']],on='gene_id',suffixes=('_factorial','_passage_aware'))
interaction_compare.to_csv(T/'interaction_model_sensitivity.csv',index=False)
section('8. Period-dependent responses',
 '<p>The batch × condition model tests the May-minus-August difference in the Y-depleted-versus-sgCTRL effect. Two genes pass interaction FDR &lt;0.05: NRROS and SERPINB1. The earlier passage-aware model yields 451 interaction-significant genes; different residual models produce a major sensitivity that precludes a firm gene-count narrative. A nonsignificant interaction does not prove stability.</p>'
 +table(newint[newint.FDR<.05][['gene_name','logFC','FDR']])+link('interaction_model_sensitivity')+
 '<p>Interaction-ranked GSEA is available separately; positive enrichment describes a more positive condition effect in May, not necessarily higher expression in either period. Negative enrichment has the converse meaning.</p>'+link('GSEA_interaction'))
over=parent[parent.directionally_consistent & parent.parental_passes_cutoffs]
section('9. Parental context',
 '<p>August sgCTRL versus parental tests 23,866 genes and identifies 12 at FDR &lt;0.05 and |log₂FC| ≥1. Five also meet the main Y-depleted comparison criteria. Opposite signs below indicate that sgCTRL-to-parental and Y-depleted-to-sgCTRL differences oppose one another; they do not establish restoration to the parental state. This comparison does not identify nucleofection-specific effects.</p>'
 +table(over[['gene_name','logFC','FDR','logFC_CTRL_parental','FDR_CTRL_parental']])+link('DE_CTRL_vs_parental')+' | '+link('GSEA_CTRL_vs_parental'))
section('10. Provisional detailed-review groups',
 '<p>Five review candidates are selected from groups with ≥2 pathways, ≥1 directionally consistent DEG and ≥1 concordant significant promoter-associated gene. Sorting prioritises both-period pathway support, then promoter and motif support; groups with leading-gene Jaccard ≥0.5 to an earlier selection are skipped. This is a transparent review aid, not a statistical ranking or a final choice of biological themes.</p>'+table(short))
for r in short.itertuples():
 genes=read(f'dossier_{r.community}_genes')
 use=genes[genes.directionally_consistent.fillna(False)|genes.ATAC_any_significant_promoter.fillna(False)]
 cols=['gene_name','logFC','FDR','AUG_logFC','MAY_logFC','median_peak_logFC','minimum_peak_FDR','gained_peaks','lost_peaks']
 section(f'Dossier: {r.community} — {r.representative_pathway}',table(use[cols])+link(f'dossier_{r.community}_genes')+' | '+link(f'dossier_{r.community}_TFs'))
section('11. Status and next review',
 '<p>Implemented: agreed specification, core DE, period interaction, parental DE, three-collection GSEA, correlation sensitivity, pathway grouping, STRING modules, TF activity, current-RNA promoter/motif integration, all-group overview and five provisional gene/TF dossiers.</p>'
 '<p>Next interpretive work: review and name the biological themes; choose the final 3–5; inspect locus-level coverage and supporting regulatory links; add relevant broader signatures as follow-up; prepare the detailed thesis panels and Results prose. Direct RNA–ATAC sample matching and distal enhancer links remain unestablished. Existing Y/PAR dosage outputs remain separate and unchanged.</p>'
 '<p>Original source measurements and earlier analysis outputs are preserved. Independent numerical validation and SHA256 manifests are in provenance/.</p>')
style='''body{font-family:"Open Sans",Arial,sans-serif;color:#222;max-width:1250px;margin:32px auto;padding:0 24px;line-height:1.5}h1,h2{color:#263e62}section{margin:30px 0;border-top:1px solid #ddd;padding-top:12px}table{border-collapse:collapse;font-size:12px}th,td{padding:6px 9px;border-bottom:1px solid #ddd;text-align:left}th{background:#edf1f7;position:sticky;top:0}.scroll{overflow:auto;max-height:580px}img.figure{max-width:850px;width:100%}.note{background:#fff4dd;border-left:4px solid #ae7a22;padding:12px}input{padding:8px;width:95%;margin:12px 0}a{color:#315788}'''
doc='<!doctype html><html><head><meta charset="utf-8"><title>HT1080 RNA programme analysis</title><style>'+style+'</style></head><body><h1>HT1080 RNA programme analysis</h1><p>Applied interview plan · 18 September 2026</p><label>Filter table rows: <input id="filter" placeholder="Gene, pathway, TF or community name"></label>'+''.join(parts)+'''<script>document.getElementById('filter').addEventListener('input',function(){let q=this.value.toLowerCase();document.querySelectorAll('table tbody tr').forEach(r=>r.style.display=r.textContent.toLowerCase().includes(q)?'':'none');});</script></body></html>'''
(O/'report.html').write_text(doc,encoding='utf-8')
readme='''# HT1080 RNA programme analysis — current implementation

Open [the searchable report](report.html). Decisions are in [the agreed plan](../../PLAN_20260918.md).

## Completed

- Core batch + condition DE; coherent period-specific/interaction fit; separate August parental comparison.
- Gene consistency, Hallmark/Reactome/GO-BP GSEA, programme consistency and leading-edge groups.
- STRING/MCODE modules, CollecTRI activity, promoter and motif evidence from existing ATAC results.
- All-theme overview, five provisional detailed-review dossiers, thesis-style overview figure, independent numerical validation.

## Main results

- 662 directionally consistent genes (229 higher, 433 lower).
- 220 directionally consistent enriched pathways; 55 significant in both periods; 88 pathway groups including isolated terms.
- 60 directionally consistent TF patterns; five STRING association modules.
- 12 sgCTRL-versus-parental genes at both cutoffs; five overlap the primary gene set.
- Two interaction-significant genes in the factorial model versus 451 in the earlier passage-aware fit: major model sensitivity.
- CAMERA sensitivity: 17 discoveries at fixed correlation 0.01, zero with estimated nonnegative correlation. GSEA interpretation remains exploratory.

## Remaining interpretive work

Review theme labels and final 3–5 selections; inspect locus-level evidence and supported enhancer links; select broader signatures for follow-up; prepare detailed thesis figures and prose. RNA–ATAC specimen matching remains unconfirmed; current integration is at group level. The Y/PAR branch remains separate and unchanged.

## Reproduce

From this analysis directory, run with Rscript --vanilla and Python:

1. run_programme_analysis_v1.R
2. build_programme_networks_v1.py
3. check_programme_camera_v1.R
4. integrate_programmes_ATAC_v1.py
5. plot_programme_analysis_v1.py
6. report_programme_analysis_v1.py

Original inputs remain unchanged; input/code hashes and validation are in provenance/. All CSV tables retain full numerical precision. Parentals are independently filtered and normalized, and no PCA batch-adjusted values enter DE.
'''
(O/'README.md').write_text(readme,encoding='utf-8')
print(json.dumps(validation,indent=2));print('Report saved:',O/'report.html')
