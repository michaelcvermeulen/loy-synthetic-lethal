"""Build and verify the focused TFAP2C review without changing upstream results."""
from pathlib import Path
import base64
import hashlib
import html
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'motifs/TF_extension/TFAP2C_review'
read = lambda name: pd.read_csv(OUT/f'{name}.csv')
ap = read('AP2_cross_method_comparison')
tf = ap[ap.motif_name.eq('TFAP2C')].copy().sort_values('motif_id')
pk = read('TFAP2C_peak_summary')
union = pk[pk.motif_id.eq('TFAP2C_any')].iloc[0]
ss = read('TFAP2C_TOBIAS_site_summary')
ss['bound_count_percent_change'] = 100*(ss.LOY_bound/ss.CTRL_bound-1)
ss.to_csv(OUT/'TFAP2C_TOBIAS_site_summary_with_percent.csv',index=False)
hom = read('AP2_HOMER_primary_and_passage')
own = read('TFAP2C_own_promoter_all_peaks')
gs = read('TFAP2C_promoter_RNA_summary')
genes = read('TFAP2C_lost_promoter_genes_thresholds')
strict = genes[genes.ATAC_FDR_cutoff.eq(.05) & genes.ATAC_abs_logFC_cutoff.eq(1) & genes.rna_FDR.lt(.05)].sort_values('rna_FDR')
strict.to_csv(OUT/'TFAP2C_joint_promoter_genes_strict.csv',index=False)
robust = read('TFAP2C_sample_robustness')
over = read('AP2_peak_overlap')
table = lambda df: df.to_html(index=False,border=0,float_format=lambda x:f'{x:.4g}',escape=True)
def figure(name,caption):
    p=OUT/f'{name}.png'
    assert p.exists() and p.stat().st_size>1000
    b64=base64.b64encode(p.read_bytes()).decode()
    return f'<figure><img src="data:image/png;base64,{b64}" alt="{html.escape(caption)}"><figcaption>{html.escape(caption)}</figcaption></figure>'

summary = f'''TFAP2C-focused review

Conclusion: LOY is associated with selective loss of accessibility at AP-2 motif-containing regulatory regions, accompanied by lower TFAP2C RNA and a modest shift in pooled predicted occupancy. This supports TFAP2C as a candidate contributor, but sequence motifs do not distinguish it from related AP-2 proteins.

All three TFAP2C models agree in direction. HOMER lost-peak enrichment is 1.89-2.12-fold, with FDR 2.0e-15 to 2.8e-7; no primary gained enrichment passes 0.05. The earlier JASPAR2022 CMH screen supports the same direction, with different motif versions and target definitions.

All LOY chromVAR scores are below all CTRL scores for every TFAP2C model. Primary limma FDR is 0.058-0.067; passage sensitivity FDR is 0.000125-0.000299, while Welch FDR is 0.255-0.276. Omitting any one sample preserves negative mean effects. Direction is robust; the significance threshold depends on statistical model. Technical batch is unconfirmed.

Across the union of {int(union.n_peaks):,} motif-positive peaks, {int(union.lost_FDR005_LFC1)} lose and {int(union.gained_FDR005_LFC1)} gain accessibility at FDR <0.05 and |log2FC| >1. Median accessibility log2FC is {union.median_all_peak_logFC:.4f}; {100*union.fraction_all_peak_negative:.1f}% of peaks have negative effects. These results do not imply uniform closure of all AP-2 regions.

TOBIAS predicts 2.1-6.1% fewer bound sites, depending on motif model, with both losses and gains at individual sites. Aggregate corrected cleavage profiles largely overlap and do not show a clear condition-specific loss of central protection. Its pooled-site significance is not a biological-replicate test and should not be treated as an independent validation experiment. Footprint-derived scores depend on local accessibility as well as cleavage patterns.

{100*union.shared_other_AP2/union.n_peaks:.2f}% of TFAP2C-positive peaks also contain another AP-2 model somewhere in the same peak. This is peak-level overlap, not proof of identical binding positions. TFAP2A RNA is unchanged statistically; TFAP2B is absent from the selected RNA table. TFAP2C RNA log2FC is -0.747 (FDR 0.0117), approximately 40.4% lower. Both annotated TFAP2C promoter peaks are nonsignificant in primary ATAC.

At the strict DA cutoff, lost motif-containing promoter peaks link to 39 RNA-mapped genes: 17 have significantly lower RNA at FDR <0.05, none higher. With both assays at FDR <0.10 and no ATAC effect-size filter, 32 are RNA-down and 2 RNA-up (BMP7 and SLC9A2). These are positional candidate links, not experimentally established TFAP2C targets or pathway enrichment.

No copy-number adjustment has been applied. The assays are compared at group level; the three ATAC-based methods reuse the same data and are not independent biological experiments.
'''
(OUT/'CONCLUSIONS.md').write_text(summary,encoding='utf-8')
comparison=tf[['motif_id','HOMER_lost_fold_enrichment','HOMER_lost_BH_FDR_both_directions','chromVAR_delta_z','chromVAR_FDR','chromVAR_passage_FDR','chromVAR_Welch_FDR','LOY_CTRL_change']]
parts=['<!doctype html><html><head><meta charset="utf-8"><title>TFAP2C: evidence across ATAC and RNA methods</title>',
 '<style>body{font:16px/1.55 system-ui,sans-serif;max-width:1150px;margin:40px auto;padding:0 24px;color:#18242b}h1,h2{line-height:1.25}h2{margin-top:40px}table{border-collapse:collapse;display:block;overflow:auto;font-size:12px}td,th{padding:7px;border-bottom:1px solid #ddd;text-align:left}th{background:#edf3f5}img{max-width:100%}figcaption{font-size:13px;color:#46545c}.lead{background:#edf3f5;padding:18px}a{color:#006b95}</style></head><body>',
 '<h1>TFAP2C: evidence across ATAC and RNA methods</h1>',
 '<p>HT1080 LOY versus CTRL · Three independent ATAC replicates per condition · Ben thesis combined-batch RNA, design batch + treatment</p>',
 '<p class="lead"><b>Selective AP-2-associated accessibility losses are supported across methods.</b> Lower TFAP2C RNA makes TFAP2C a candidate contributor. The data do not establish a unique TFAP2C mechanism or wholesale loss of AP-2 binding.</p>',
 '<h2>1. Three TFAP2C motif models agree, with different statistical strength</h2>',table(comparison),
 '<p>HOMER uses gained/lost peaks at primary edgeR FDR &lt;0.05 and |log2FC| &gt;1, with GC/accessibility-matched non-DA backgrounds. Its BH correction covers all motifs and both directions. No primary gained TFAP2C enrichment passes FDR 0.05. chromVAR uses all eligible peaks and sample-level scores; its effect is a z-score difference, not log2FC. TOBIAS change is the native aggregate score, not percent occupancy.</p>',
 figure('TFAP2C_enrichment_and_binding','HOMER model sensitivity and native predicted bound-site counts. Passage-adjusted DA is a sensitivity arm, not a confirmed technical-batch adjustment.'),
 '<details><summary>HOMER target and matched-background fractions</summary>',table(hom[hom.motif_name.eq('TFAP2C')][['motif_id','arm','direction','target_motif_count','target_total','target_fraction','background_fraction','fold_enrichment','BH_FDR_both_directions']]),'</details>',
 '<details><summary>Original JASPAR2022 stratified CMH screen</summary>',table(read('AP2_original_JASPAR2022_screen')),'<p>Odds ratios are not HOMER fold enrichment. Earlier motif versions and target definitions differ; this is a robustness comparison within the same experiment.</p></details>',
 '<h2>2. Replicate consistency and statistical sensitivity</h2>',
 figure('TFAP2C_sample_scores','Every LOY sample has a lower chromVAR score than every CTRL sample for each TFAP2C model. Scores are standardized against chromVAR backgrounds, not across the six samples.'),
 table(robust),
 '<p>Primary limma FDR is 0.058–0.067, passing 0.10 but not 0.05. Passage adjustment strengthens significance while retaining the effect direction. Unmoderated Welch tests do not pass FDR 0.10 (FDR 0.255–0.276). Leave-one-sample-out mean differences stay negative; these are descriptive checks, not refitted significance tests. The two mapped AP-2 group unions both have primary chromVAR FDR 0.0633.</p>',
 '<details><summary>Group-union chromVAR statistics and membership</summary>',table(read('AP2_group_chromVAR')),table(read('AP2_group_members')),'</details>',
 '<h2>3. Losses concentrate in a subset of AP-2-containing peaks</h2>',table(pk),
 '<p>The union contains 20,306 motif-positive peaks: 298 lost and 14 gained at FDR &lt;0.05 and |log2FC| &gt;1. Without the effect-size filter the counts are 326 lost and 26 gained. The median log2FC across all motif-positive peaks is approximately zero. Thus the strong lost-peak enrichment should not be interpreted as uniform closure of the entire AP-2 program.</p>',
 figure('TFAP2C_peak_volcano','Accessibility effects across the union of the three TFAP2C motif models; peaks are counted once.'),
 table(read('TFAP2C_feature_counts')),
 '<h2>4. Footprinting: a modest shift with many sites retained</h2>',table(ss),
 '<p>For MA0524.3, CTRL has 1,839 predicted bound sites and LOY 1,726 (6.1% fewer). There are 1,528 sites bound in both, 311 CTRL-only and 198 LOY-only. The other models show 2.1% and 5.2% fewer bound sites. Approximately 56–58% of native site log2FC values are negative; the shift is not universal. The three site sets overlap and should not be added together.</p>',
 figure('TFAP2C_aggregate_footprints','Mean cleavage profiles at the same full motif-site set in both conditions, before and after Tn5 bias correction. Y-axis scales match within each motif row. These pooled profiles have no biological-replicate confidence interval.'),
 '<p><b>The aggregate corrected cleavage profiles largely overlap and do not show a clear condition-specific loss of central protection.</b> The profiles contain complex central cleavage structure rather than a simple protected dip. Thus the native score direction provides limited supporting evidence, and the aggregate plot does not independently demonstrate reduced TFAP2C binding. Lower pooled scores can reflect changes in surrounding accessibility as well as occupancy. Very small pooled-site P-values do not establish a large effect or replace biological-replicate inference.</p>',
 '<h2>5. How specific is the signal to TFAP2C?</h2>',
 f'<p>{int(union.shared_other_AP2):,} of {int(union.n_peaks):,} peaks ({100*union.shared_other_AP2/union.n_peaks:.2f}%) also contain another AP-2 family motif model. This measures overlap within 501 bp peaks, not exact co-occupancy at the same base. All tested AP-2 family motif models shift in the same direction. RNA supports TFAP2C specifically as a candidate, but motif matching cannot assign the protein responsible.</p>',
 table(ap[['motif_id','motif_name','archetype','HOMER_lost_fold_enrichment','chromVAR_FDR','LOY_CTRL_change','rna_log2FC','rna_FDR','RNA_statistic_status']]),
 '<h2>6. TFAP2C RNA, its own promoter, and candidate downstream genes</h2>',
 '<p>TFAP2C RNA is approximately 40% lower (log2FC −0.747; FDR 0.0117). Its two annotated promoter peaks have ATAC FDR 0.988 and 0.307, so there is no significant primary promoter accessibility loss at TFAP2C itself. RNA downregulation need not arise from closure of these assayed promoter peaks.</p>',table(own[['peak_id','logFC','FDR','TFAP2C_motif_positive','rna_logFC','rna_FDR']]),
 table(gs),
 '<p>The strict subset uses ATAC FDR &lt;0.05 and log2FC &lt;−1 plus RNA FDR &lt;0.05: 17 genes have lower RNA, none higher. The relaxed 0.10/0.10 subset without an ATAC fold-change filter has 32 RNA-down and two RNA-up genes. Promoter overlap is a positional association; these are candidate genes, not validated direct TFAP2C targets. No new pathway enrichment claim is made.</p>',
 figure('TFAP2C_promoter_RNA_genes','Up to 25 genes with lost motif-containing promoter peaks and significant RNA at FDR 0.10, ranked by RNA FDR. The ATAC value summarizes only selected lost promoter peaks, not all regulatory peaks for the gene.'),
 '<details open><summary>Strict candidate genes</summary>',table(strict[['rna_symbol','rna_logFC','rna_FDR','minimum_peak_FDR','median_lost_peak_logFC','peaks']]),'</details>',
 '<h2>Interpretation and reproducibility</h2>',
 '<p><b>Best-supported conclusion:</b> LOY is associated with selective accessibility losses at AP-2 motif-containing regulatory regions, accompanied by reduced TFAP2C transcription and a modest decrease in pooled predicted occupancy. TFAP2C is a plausible contributor; a TFAP2C-specific causal mechanism remains untested.</p>',
 '<p>Technical batch is unconfirmed, passage sensitivity is substantial, and no copy-number adjustment has been applied. ATAC methods reuse the same underlying experiment; RNA comparisons are at group level. See METHODS.txt and verification.txt in the output folder for exact accounting.</p>',
 '<p><a href="../motifs/TF_extension/TFAP2C_review/">CSV tables, PNG/PDF figures, native aggregates and methods</a> · <a href="ATACseq_report.html">Main ATAC report</a></p>',
 '</body></html>']
report=ROOT/'report/TFAP2C_review.html'
report.write_text('\n'.join(parts),encoding='utf-8')
snippet='<section id="tfap2c-focused-review"><h2>Focused TFAP2C review</h2><p>All three TFAP2C motif models support selective AP-2-associated accessibility losses. Primary chromVAR FDR is 0.058–0.067; pooled predicted bound-site counts decrease by 2–6%. TFAP2C RNA is lower, while its annotated promoter peaks are not significantly less accessible. <a href="TFAP2C_review.html">Open the detailed cross-method review, replicate plots, footprint profiles and candidate genes.</a></p></section>'
main=ROOT/'report/ATACseq_report.html'
text=main.read_text(encoding='utf-8')
if 'id="tfap2c-focused-review"' not in text:
    assert text.count('</body>')==1
    main.write_text(text.replace('</body>',snippet+'\n</body>'),encoding='utf-8')
# Independent output accounting; no rerunning or changing upstream statistics.
positive=read('TFAP2C_motif_positive_peaks')
assert positive.peak_id.is_unique and len(positive)==union.n_peaks
assert ((positive.FDR<.05)&(positive.logFC< -1)).sum()==union.lost_FDR005_LFC1
assert ((positive.FDR<.05)&(positive.logFC>1)).sum()==union.gained_FDR005_LFC1
assert all(ss.CTRL_bound==ss.both_bound+ss.CTRL_only)
assert all(ss.LOY_bound==ss.both_bound+ss.LOY_only)
assert all(ss.total_sites==ss.both_bound+ss.CTRL_only+ss.LOY_only+ss.neither)
assert len(strict)==17 and (strict.rna_logFC<0).all()
for row in gs.itertuples():
    g=genes[genes.ATAC_FDR_cutoff.eq(row.FDR_cutoff)&genes.ATAC_abs_logFC_cutoff.eq(row.ATAC_abs_logFC_cutoff)]
    assert g.rna_ensembl_id.is_unique and len(g)==row.linked_RNA_genes
    assert ((g.rna_FDR<row.FDR_cutoff)&(g.rna_logFC<0)).sum()==row.RNA_down_same_FDR
    assert ((g.rna_FDR<row.FDR_cutoff)&(g.rna_logFC>0)).sum()==row.RNA_up_same_FDR
rna=pd.read_csv(ROOT/'config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv')
g=genes.merge(rna[['gene_id','logFC','FDR']],left_on='rna_ensembl_id',right_on='gene_id',validate='many_to_one')
assert len(g)==len(genes) and np.allclose(g.rna_logFC,g.logFC) and np.allclose(g.rna_FDR,g.FDR)
curves=read('TFAP2C_aggregate_curves')
assert len(curves)==360 and np.isfinite(curves.select_dtypes('number')).all().all()
for name in ['TFAP2C_sample_scores','TFAP2C_enrichment_and_binding','TFAP2C_peak_volcano','TFAP2C_promoter_RNA_genes','TFAP2C_aggregate_footprints']:
    assert base64.b64encode((OUT/f'{name}.png').read_bytes()).decode() in report.read_text()
with (OUT/'verification.txt').open('a') as handle:
    handle.write('PASS: Python independently verified unique peak counts, native bound-state conservation, all threshold gene counts, RNA source values and embedded current figure bytes.\n')
paths=[ROOT/'code'/n for n in ['09_report.py','34_tfap2c_review.R','35_tfap2c_footprints.py','36_tfap2c_report.py']]
paths += [ROOT/'motifs/TF_extension/chromVAR/motif_matches.rds',ROOT/'motifs/TF_extension/concordance/TF_ATAC_footprint_RNA.csv',ROOT/'config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv']
paths += [ROOT/p for p in ['differential/primary_TMM_all.csv','annotation/annotated_differential_peaks.csv','integration/RNAseq/ATAC_RNA_promoter_peak_links.csv','motifs/TF_extension/chromVAR/peak_accounting.csv','motifs/TF_extension/HOMER/known_motif_all_results.csv','motifs/TF_extension/TOBIAS/BINDetect/bindetect_results.txt']]
pd.DataFrame([dict(path=str(p.relative_to(ROOT)),SHA256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest()) for p in paths]).to_csv(OUT/'source_code_checksums.csv',index=False)
(OUT/'ANALYSIS_COMPLETE.json').write_text(json.dumps(dict(completed_UTC=datetime.now(timezone.utc).isoformat(),report=str(report.relative_to(ROOT)),data_checks='passed',primary_design='independent replicates',technical_batch='unconfirmed',copy_number_adjustment=False),indent=2)+'\n')
print(f'Created and verified {report}',flush=True)
