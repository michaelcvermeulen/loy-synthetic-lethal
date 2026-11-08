"""Validate interval joins, effect exports and report artifacts independently."""
from pathlib import Path
from html.parser import HTMLParser
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'integration/RNA_interaction_CN'
model=json.loads((OUT/'model_verification.json').read_text())
assert model['status']=='passed_with_reported_boundary_gene'
b=pd.read_csv(OUT/'WGS_hg38_100kb_bins.csv')
p=pd.read_csv(OUT/'ATAC_promoter_peak_CN_context.csv')
j=pd.read_csv(OUT/'ATAC_RNA_CN_joint_gene_table.csv',dtype={'entrez_id':str})
rna=pd.read_csv(OUT/'RNA_interaction_all_genes.csv',dtype={'entrezid':str})
old=pd.read_csv(ROOT/'integration/feature_scatters/feature_gene_effects.csv',dtype={'entrez_id':str})
old=old[old.feature.eq('promoter')].set_index('entrez_id')
assert len(j)==11898 and not j.entrez_id.duplicated().any()
for c in ['atac_median_logFC','rna_logFC','rna_FDR','minimum_peak_FDR']:
    np.testing.assert_allclose(j[c],old.loc[j.entrez_id,c],equal_nan=True)
uniquerna=rna[~rna.entrezid.duplicated(False)].set_index('entrezid')
mapped=j[j.entrez_id.isin(uniquerna.index)]
for c in ['interaction_log2FoldChange','interaction_padj','August_log2FoldChange','May_log2FoldChange']:
    np.testing.assert_allclose(mapped[c],uniquerna.loc[mapped.entrez_id,c],equal_nan=True)
assert j.Both_ATAC010_RNA005.sum()==82 and j.down_down.sum()==70 and j.up_up.sum()==9
assert b.loc[~b['use'],'delta_segment_log2'].isna().all()
for chrom,g in b.groupby('chrom'):
    assert (g.end0.to_numpy()[:-1]<=g.start0.to_numpy()[1:]).all()
    assert ((g.end0-g.start0)>0).all()
valid=b[b['use']].copy();assert len(valid)==24994
# Brute-force interval intersections for every exported promoter peak. This is
# independent of the binary-search routine in the report generator.
bychrom={c:g for c,g in valid.groupby('chrom')}
for r in p.itertuples():
    g=bychrom.get(r.chrom,pd.DataFrame(columns=b.columns))
    overlap=g[(g.start0<r.end)&(g.end0>r.start)]
    if len(overlap):
        w=np.minimum(overlap.end0,r.end)-np.maximum(overlap.start0,r.start)
        v=overlap.delta_segment_log2
        np.testing.assert_allclose([r.CN_covered_fraction,r.CN_delta_mean,r.CN_delta_min,r.CN_delta_max,r.CN_delta_max_abs],[w.sum()/(r.end-r.start),np.average(v,weights=w),v.min(),v.max(),abs(v).max()],rtol=1e-10,atol=1e-12)
    else:assert r.CN_covered_fraction==0 and np.isnan(r.CN_delta_mean)
links=pd.read_csv(OUT/'promoter_gene_peak_CN_links.csv',dtype={'entrez_id':str})
bygene=j.set_index('entrez_id');known_genes=set(j.entrez_id)
for gid,g in links.groupby('entrez_id'):
    if gid not in known_genes:continue
    row=bygene.loc[gid]
    median=g.CN_delta_mean.dropna().median() if g.CN_delta_mean.notna().any() else np.nan
    np.testing.assert_allclose([row.CN_promoter_min_coverage,row.CN_promoter_delta_median,row.CN_promoter_max_abs_delta],[g.CN_covered_fraction.min(),median,g.CN_delta_max_abs.max()],equal_nan=True)
for t in [.1,.2,.3]:
    expected=np.where(j.CN_promoter_max_abs_delta>=t,'D0 dosage difference',np.where(j.CN_promoter_min_coverage>=.9,'No large D0 difference','Insufficient coverage'))
    assert (j[f'CN_flag_{t:.2f}']==expected).all()
key=j[j.rna_symbol.isin(['SAMD9L','TNFAIP2','MX1','BST2'])]
assert len(key)==4 and (key.CN_promoter_min_coverage==1).all() and (key.CN_promoter_max_abs_delta<.03).all()
chr19=j[j.down_down&j.gene_body_chrom.eq('chr19')]
assert len(chr19)==19 and chr19['CN_flag_0.20'].eq('No large D0 difference').sum()==13 and chr19['CN_flag_0.20'].eq('Insufficient coverage').sum()==6
lrt=pd.read_csv(OUT/'RNA_interaction_LRT_sensitivity.csv')
finite=np.isfinite(lrt.LRT_pvalue)
np.testing.assert_allclose(lrt.loc[finite,'LRT_pvalue'],chi2.sf(lrt.loc[finite,'LRT_stat'],df=1),rtol=1e-9,atol=1e-14)
both=pd.read_csv(OUT/'RNA_interaction_supported_both_tests.csv')
assert len(both)==52 and 'LENG1' not in set(both.gene_name)
class Links(HTMLParser):
    def handle_starttag(self,tag,attrs):
        for k,v in attrs:
            if k in ['href','src'] and v and not v.startswith(('https:','http:','#','data:')):
                assert (ROOT/'report'/v.split('#')[0]).resolve().is_file(),v
Links().feed((ROOT/'report/RNA_interaction_CN.html').read_text(encoding='utf-8'))
for stem in ['RNA_interaction_volcano','RNA_key_gene_interactions','WGS_genome_dosage_context','WGS_key_gene_regions','CN_pathway_overlap_stratification']:
    for ext in ['pdf','png','svg']:assert (OUT/'figures'/f'{stem}.{ext}').stat().st_size>1000
completion=json.loads((OUT/'completion.json').read_text())
for rel,digest in completion['source_sha256'].items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest
result={'status':'passed','promoter_peaks_independently_overlapped':len(p),'joint_genes':len(j),'original_ATAC_RNA_effects_unchanged':True,'all_CN_threshold_flags_verified':True,'local_report_links_verified':True,'visual_review':'Volcano, interaction intervals, genome profile and key-gene WGS panels reviewed','model_verification':model}
(OUT/'report_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
completion.update(status='complete',verification=result)
completion['source_sha256'][str(Path(__file__).relative_to(ROOT))]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(OUT/'completion.json').write_text(json.dumps(completion,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
