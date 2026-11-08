"""Group-level promoter and motif evidence for current RNA programme themes."""
from pathlib import Path
import json,hashlib,re
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'outputs/programme_analysis_v1';T=OUT/'tables'
ATAC=ROOT.parents[3]/'ATACseq/SEPT_08_2026_REPROCESS'
def save(x,name):x.to_csv(T/f'{name}.csv',index=False)
sources=[ATAC/'annotation/promoter_peak_gene_links.csv',ATAC/'annotation/promoter_gene_summary.csv',
 ATAC/'motifs/TF_extension/chromVAR/motifs_independent_limma.csv',
 ATAC/'motifs/TF_extension/chromVAR/motifs_passage_sensitivity_limma.csv',ATAC/'motifs/motif_enrichment.csv',ATAC/'config/samples.tsv']
rna=pd.read_csv(T/'gene_consistency.csv')
reps=pd.read_csv(ROOT/'outputs/tables/symbol_mapping_representatives.csv')
rna=rna[rna.gene_id.isin(reps.gene_id)].copy()
assert rna.gene_name.is_unique
prom=pd.read_csv(sources[0]); summary=pd.read_csv(sources[1])
# Exclude ambiguous symbol-to-Entrez mappings; retain one-to-many peak links.
mapping=summary[['symbol','gene_id']].drop_duplicates()
amb=set(mapping.loc[mapping.symbol.duplicated(keep=False),'symbol'])
save(mapping[mapping.symbol.isin(amb)],'ATAC_ambiguous_symbols')
summary=summary[~summary.symbol.isin(amb)].copy()
summary=summary.rename(columns={'gene_id':'ATAC_entrez_id','symbol':'gene_name'})
assert summary.gene_name.is_unique
joined=rna.merge(summary,on='gene_name',how='left',validate='one_to_one')
joined['ATAC_promoter_tested']=joined.tested_promoter_peaks.notna()
joined['ATAC_any_significant_promoter']=(joined.gained_peaks.fillna(0)+joined.lost_peaks.fillna(0))>0
joined['ATAC_concordant_peak']=((joined.logFC>0)&(joined.gained_peaks>0))|((joined.logFC<0)&(joined.lost_peaks>0))
joined['ATAC_discordant_peak']=((joined.logFC>0)&(joined.lost_peaks>0))|((joined.logFC<0)&(joined.gained_peaks>0))
save(joined,'RNA_ATAC_promoter_gene_evidence')
peaks=prom[~prom.symbol.isin(amb)].merge(rna[['gene_name','gene_id','logFC','FDR','directionally_consistent']],
 left_on='symbol',right_on='gene_name',suffixes=('_ATAC','_RNA'),validate='many_to_one')
save(peaks,'RNA_ATAC_promoter_peak_evidence')

tf=pd.read_csv(T/'TF_consistency.csv');motif=pd.read_csv(sources[2]);paired=pd.read_csv(sources[3])
paired=paired[['id','LOY_minus_CTRL_zscore','adj.P.Val']].rename(columns={'LOY_minus_CTRL_zscore':'passage_adjusted_effect','adj.P.Val':'passage_adjusted_FDR'})
motif=motif.merge(paired,on='id',how='left',validate='one_to_one')
known=set(tf.source);upper={}
for k in known:upper.setdefault(k.upper(),[]).append(k)
links=[]
for r in motif.to_dict('records'):
 for member in re.split(r'::',str(r['motif_name'])):
  hits=[member] if member in known else upper.get(member.upper(),[])
  if len(hits)!=1:continue
  links.append({**r,'source':hits[0],'mapping_type':'exact' if member in known else 'case_insensitive_unique','motif_member':member})
ml=pd.DataFrame(links).merge(tf,on='source',how='left',suffixes=('_ATAC','_RNA'))
ml['ATAC_motif_FDR05']=ml['adj.P.Val']<.05
ml['RNA_ATAC_direction_agrees']=np.sign(ml.score)==np.sign(ml.LOY_minus_CTRL_zscore)
save(ml,'TF_RNA_chromVAR_evidence')
enrichment=pd.read_csv(sources[4]);em=[]
for r in enrichment.to_dict('records'):
 for member in re.split(r'::',str(r['motif_name'])):
  hits=[member] if member in known else upper.get(member.upper(),[])
  if len(hits)==1:em.append({**r,'source':hits[0]})
save(pd.DataFrame(em).merge(tf,on='source',how='left',suffixes=('_motif','_RNA')),'TF_RNA_sequence_enrichment')

communities=pd.read_csv(T/'programme_communities.csv');drivers=pd.read_csv(T/'programme_gene_drivers.csv')
tf_links=pd.read_csv(T/'TF_programme_links.csv')
gene_evidence=drivers.merge(joined,on='gene_name',how='left',validate='many_to_one');save(gene_evidence,'programme_gene_ATAC_evidence')
rows=[]
for c in communities.itertuples():
 g=gene_evidence[gene_evidence.community==c.community]
 ts=set(tf_links.loc[(tf_links.community==c.community)&(tf_links.target_overlap>=3),'TF'])
 mt=ml[ml.source.isin(ts)&ml.directionally_consistent]
 supported=mt[mt.ATAC_motif_FDR05 & mt.RNA_ATAC_direction_agrees]
 rows.append(dict(community=c.community,pathways=c.pathways,representative_pathway=c.representative_pathway,direction=c.direction,
   both_periods_supported=c.both_periods_supported,leading_genes=len(g),consistent_DE_genes=int(g.directionally_consistent.fillna(False).sum()),
   promoter_tested_genes=int(g.ATAC_promoter_tested.fillna(False).sum()),promoter_concordant_genes=int(g.ATAC_concordant_peak.fillna(False).sum()),
   promoter_discordant_genes=int(g.ATAC_discordant_peak.fillna(False).sum()),consistent_linked_TFs=len(ts),
   concordant_chromVAR_TFs=supported.source.nunique(),concordant_chromVAR_TF_names=';'.join(sorted(set(supported.source)))))
review=pd.DataFrame(rows)
review=review.sort_values(['both_periods_supported','promoter_concordant_genes','concordant_chromVAR_TFs','consistent_DE_genes','community'],ascending=[False,False,False,False,True])
save(review,'programme_integration_overview')
# Candidate dossiers are a review aid, not an automatically endorsed biological ranking.
eligible=review[(review.promoter_concordant_genes>0)&(review.consistent_DE_genes>0)&(review.pathways>=2)]
selected=[];selected_sets=[]
for r in eligible.itertuples():
 geneset=set(drivers.loc[drivers.community==r.community,'gene_name'])
 if any(len(geneset&s)/len(geneset|s)>=.5 for s in selected_sets):continue
 selected.append(r.community);selected_sets.append(geneset)
 if len(selected)==5:break
short=eligible[eligible.community.isin(selected)].copy()
short['status']='provisional_review_candidate_not_final_theme_selection'
save(short,'provisional_theme_shortlist')
for c in short.community:
 save(gene_evidence[gene_evidence.community==c],f'dossier_{c}_genes')
 save(tf_links[tf_links.community==c],f'dossier_{c}_TFs')

manifest={'sources':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources],
 'assembly':'RNA GENCODE v49 GRCh38; ATAC hg38/UCSC knownGene; unique exact gene-symbol joins',
 'link':'Existing tested peak overlaps TSS +/-2kb; all links retained. No distal enhancer assignment inferred.',
 'sample_relation':'Group-level integration; direct RNA/ATAC specimen and time matching unresolved.',
 'ATAC_model':'Existing primary ~condition comparison; preserve passage-adjusted chromVAR sensitivity.',
 'gene_statistics':'Promoter minimum FDR is a peak FDR, not a gene-level FDR; gained and lost peaks can coexist.',
 'TF_mapping':'Exact or uniquely case-insensitive motif-member symbols; motif families and dimers do not identify unique binding.',
 'shortlist_rule':'Review candidates only: >=2 pathways, >=1 consistent DEG, >=1 concordant promoter gene; sort by both-period-supported pathways then promoter support, motif support and DE gene count; skip leading-gene Jaccard >=0.5 to earlier choices. Linked TF evidence requires >=3 programme targets. No minimum evidence is a new statistical test.',
 'counts':{'RNA_symbols':len(rna),'promoter_tested':int(joined.ATAC_promoter_tested.sum()),'ambiguous_ATAC_symbols':len(amb),'themes':len(review)}}
(OUT/'provenance/ATAC_integration.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps(manifest['counts']))
