"""Freeze descriptive Figure 2 inputs without refitting drugZ or MAGeCK."""
from pathlib import Path
import sys, json, hashlib
THESIS = Path(__file__).resolve().parents[3]
DEPS = THESIS.parent / '.figure_python_deps'
if DEPS.is_dir(): sys.path.insert(0, str(DEPS))
import numpy as np
import pandas as pd
SCREEN = THESIS.parent / 'crispr_screen/ht1080_loy_screen_2'
OUT = THESIS / 'figures/chapter4/source_data/ht1080_candidate_discovery'
SOURCES = {
    'drugz': SCREEN / 'results/drugZ/LOY_vs_CTRL_unfiltered.txt',
    'mageck': SCREEN / 'test/POST_LOY_vs_CTRL.gene_summary.txt',
    'mageck_guides': SCREEN / 'test/POST_LOY_vs_CTRL.sgrna_summary.txt',
    'counts': SCREEN / 'count/ht1080_screen2_R1.count.txt',
    'library': SCREEN / 'brunello/broadgpp-brunello-library-corrected.txt',
    'previous_log_cpm': THESIS / 'figures/chapter4/source_data/ht1080_figure1_qc/guide_log2_cpm.csv',
    'ppi_network': SCREEN / 'results/metascape/depleted/Enrichment_PPI/PPINetwork.js',
    'ppi_input': SCREEN / 'results/metascape/depleted/metascape_result.xlsx',
}
SAMPLES = ['D0_CTRL','D0_LOY','POST_CTRL_REP1','POST_CTRL_REP2','POST_LOY_REP1','POST_LOY_REP2']
FOLLOWUP = ['DDX3X','SREBF1','TADA2B','TAF6L']
# Supp2b D0 heatmap: ten top-ranked drugZ sensitizers plus TADA2B (author request, 22 Oct 2026).
PANEL_B = ['DDX3X','TIMM8A','MED12','SCLY','WTAP','SREBF2','SREBF1','LLGL1','TAF6L','UBE2G2','TADA2B']
EXTENDED = ['SREBF2','SUPT20H','ATXN7L3','TADA1','TIMM8A','MED12','EIF4A3']
OTHER = ['EIF1AX','ZFX','RPS4X','KDM5C','KDM6A']
GAMETOLOGS = ['DDX3X','EIF1AX','RPS4X','ZFX']
SAGA = {
    'Core':['TAF5L','SUPT20H','TAF6L','TAF9B','TAF10','TAF12','SUPT7L','TADA1','SUPT3H'],
    'HAT':['KAT2A','TADA2B','TADA3','SGF29'],
    'DUB':['ATXN7','ATXN7L3','USP22','ENY2'],
    'TRRAP':['TRRAP'],
    'Splicing':['SF3B3','SF3B5'],
    'Alternatives':['KAT2B','TAF9','ATXN7L1','ATXN7L2'],
}

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    dz = pd.read_csv(SOURCES['drugz'], sep='\t')
    assert dz.GENE.is_unique and len(dz)==19113
    dz = dz.loc[dz.GENE.ne('Non-Targeting-Control')].copy()
    assert len(dz)==19112 and np.isfinite(dz.normZ).all()
    dz['directional_fdr'] = np.where(dz.normZ.lt(0), dz.fdr_synth, dz.fdr_supp)
    dz['display_class'] = np.select([dz.directional_fdr.lt(.05),dz.directional_fdr.lt(.25)],['FDR < 0.05','0.05–0.25'],default='FDR ≥ 0.25')
    dz.to_csv(OUT/'all_gene_scores.csv',index=False)
    depleted=dz.loc[dz.fdr_synth.lt(.25)].sort_values('rank_synth')
    enriched=dz.loc[dz.fdr_supp.lt(.25)].sort_values('rank_supp')
    assert len(depleted)==48 and len(enriched)==9
    assert depleted.fdr_synth.lt(.05).sum()==28 and enriched.fdr_supp.lt(.05).sum()==2
    depleted.to_csv(OUT/'depleted_candidates.csv',index=False)
    enriched.to_csv(OUT/'enriched_candidates.csv',index=False)
    mg=pd.read_csv(SOURCES['mageck'],sep='\t').rename(columns={'id':'GENE'})
    assert mg.GENE.is_unique
    join=dz.merge(mg, on='GENE', how='inner',validate='one_to_one')
    join.to_csv(OUT/'method_comparison.csv',index=False)
    supports=join.loc[join.fdr_synth.lt(.25)]
    assert len(supports)==48 and supports['neg|fdr'].lt(.25).sum()==21
    assert supports['neg|fdr'].lt(.05).sum()==9
    rank_r=join[['rank_synth','neg|rank']].rank().corr().iloc[0,1]

    counts=pd.read_csv(SOURCES['counts'],sep='\t')
    assert len(counts)==77441 and counts.sgRNA.is_unique
    cpm=counts[SAMPLES].div(counts[SAMPLES].sum(),axis=1)*1e6
    log=np.log2(cpm+1)
    saved=pd.read_csv(SOURCES['previous_log_cpm'])
    assert counts.sgRNA.equals(saved.sgRNA)
    delta=float(np.max(np.abs(log.to_numpy()-saved[SAMPLES].to_numpy())))
    assert delta<1e-12
    selected=counts.loc[counts.Gene.isin(FOLLOWUP+EXTENDED+OTHER+PANEL_B)].copy()
    selected=selected.sort_values(['Gene','sgRNA'])
    selected['guide_number']=selected.groupby('Gene').cumcount()+1
    selected['all_samples_zero']=selected[SAMPLES].sum(axis=1).eq(0)
    library=pd.read_csv(SOURCES['library'],sep='\t')
    library_pairs=set(zip(library.gene,library.Seq))
    assert selected.groupby('Gene').size().eq(4).all()
    assert all((r.Gene,r.sgRNA.rsplit('_',1)[-1]) in library_pairs for r in selected.itertuples())
    selected[['Gene','sgRNA','guide_number','all_samples_zero']].to_csv(OUT/'guide_key.csv',index=False)
    records=[]
    for idx,row in selected.iterrows():
        for sample in SAMPLES:
            records.append(dict(gene=row.Gene,sgrna=row.sgRNA,guide_number=int(row.guide_number),
                sample=sample,genotype='CTRL' if 'CTRL' in sample else 'LOY',
                timepoint='D0' if sample.startswith('D0') else 'Endpoint',
                replicate=0 if sample.startswith('D0') else int(sample[-1]),
                raw_count=int(row[sample]),cpm=float(cpm.loc[idx,sample]),
                log2_cpm_plus_1=float(log.loc[idx,sample]),all_samples_zero=bool(row.all_samples_zero)))
    profiles=pd.DataFrame(records)
    profiles.to_csv(OUT/'guide_profiles.csv',index=False)
    assert profiles.groupby(['gene','sgrna']).size().eq(6).all()
    # Explicit sequence/count table for the four established gametologs.
    gm=selected.loc[selected.Gene.isin(GAMETOLOGS)].copy()
    gm['sequence_5prime_to_3prime']=gm.sgRNA.str.rsplit('_',n=1).str[-1]
    assert gm.sequence_5prime_to_3prime.str.fullmatch('[ACGT]{20}').all()
    library=pd.read_csv(SOURCES['library'],sep='\t')
    library_pairs=set(zip(library.gene,library.Seq))
    assert all((r.Gene,r.sequence_5prime_to_3prime) in library_pairs for r in gm.itertuples())
    gm['gene_order']=gm.Gene.map({g:i for i,g in enumerate(GAMETOLOGS)})
    gm=gm.sort_values(['gene_order','guide_number'])
    mg_guides=pd.read_csv(SOURCES['mageck_guides'],sep='\t')
    assert mg_guides.sgrna.is_unique
    gm=gm.merge(mg_guides[['sgrna','Gene','LFC','FDR']],left_on=['sgRNA','Gene'],right_on=['sgrna','Gene'],
        how='left',validate='one_to_one').rename(columns={'LFC':'mageck_reported_LFC','FDR':'mageck_guide_FDR'})
    assert len(gm)==16 and gm.mageck_reported_LFC.notna().all()
    gm['endpoint_detected_samples']=gm[SAMPLES[2:]].gt(0).sum(axis=1)
    gm['mageck_LFC_displayed']=~gm.all_samples_zero
    guide_columns=['Gene','guide_number','sgRNA','sequence_5prime_to_3prime']+SAMPLES+[
        'endpoint_detected_samples','all_samples_zero','mageck_reported_LFC','mageck_guide_FDR','mageck_LFC_displayed']
    gm[guide_columns].to_csv(OUT/'gametolog_guide_table.csv',index=False)
    table=['# HT1080 gametolog guide sequences and performance','',
        'All four Brunello guides per gene. Guide numbers match the figures. Sequences are reported 5′–3′ from the count-matrix identifier and exclude the PAM. Counts are unnormalised sequencing reads. The CSV contains all six sample counts and the saved MAGeCK guide statistics.','',
        '| Gene | Guide | Sequence (5′–3′) | D0 CTRL reads | D0 LOY reads | Endpoint samples detected / 4 |',
        '|---|---:|---|---:|---:|---:|']
    for row in gm.itertuples():
        table.append(f'| {row.Gene} | {row.guide_number} | `{row.sequence_5prime_to_3prime}` | {row.D0_CTRL:,} | {row.D0_LOY:,} | {row.endpoint_detected_samples} |')
    table.extend(['','Detection means at least one read, not confirmed cutting activity. DDX3X guide 4 has zero reads in all samples; its software-reported MAGeCK LFC is retained in the CSV for provenance but omitted from the figure. No sequence-specificity analysis or editing-efficiency assay was performed for this display.',''])
    (OUT/'gametolog_guide_table.md').write_text('\n'.join(table),encoding='utf-8')

    # Reuse the saved network as a supporting MAGeCK-selected display; do not
    # mislabel it as a complete PPI analysis of the 48 drugZ candidates.
    ppi_input=pd.read_excel(SOURCES['ppi_input'],sheet_name='Annotation',usecols=['MyList','Gene ID','Gene Symbol'])
    assert set(ppi_input['Gene Symbol'])==set(mg.loc[mg['neg|fdr'].lt(.25),'GENE'])
    ppi_input['drugz_candidate']=ppi_input['Gene Symbol'].isin(depleted.GENE)
    ppi_input.to_csv(OUT/'ppi_input_audit.csv',index=False)
    raw=SOURCES['ppi_network'].read_text(encoding='utf-8')
    networks=json.loads(raw[raw.index('{'):raw.rindex('}')+1])
    net=networks['MyList_PPIColorByCluster']['elements']
    node_rows=[]
    for n in net['nodes']:
        d=n['data']; gene=d['Symbol']
        score=dz.set_index('GENE').loc[gene]
        node_rows.append(dict(gene=gene,entrez_id=d['Gene'],original_node_id=d['id'],
            original_x=n['position']['x'],original_y=n['position']['y'],
            saved_mcode_cluster=d['MCODE_CLUSTER_ID'],normZ=score.normZ,fdr_synth=score.fdr_synth))
    nodes=pd.DataFrame(node_rows)
    assert len(nodes)==7 and nodes.fdr_synth.lt(.25).all()
    edges=[]
    for e in net['edges']:
        d=e['data']; support=json.loads(d['support'])
        edges.append(dict(gene_a=d['Name_A'],gene_b=d['Name_B'],
            string_physical_score=support['String_Score_Physical']/1000,
            string_combined_score=support['String_Score_Combined']/1000,
            string_evidence=support.get('StringEvidence',''),biogrid_type=support.get('BioGrid_Type',''),
            biogrid_pubmed=support.get('BioGrid_PubMed',''),saved_support_json=d['support']))
    edges=pd.DataFrame(edges)
    assert len(edges)==11 and edges.string_physical_score.ge(.7).sum()==10
    nodes.to_csv(OUT/'ppi_nodes.csv',index=False);edges.to_csv(OUT/'ppi_edges.csv',index=False)
    saga=[]
    indexed=dz.set_index('GENE')
    for module,genes in SAGA.items():
        for gene in genes:
            source_gene='CCDC101' if gene=='SGF29' else gene
            scored=source_gene in indexed.index
            record=dict(gene=gene,source_gene=source_gene,module=module,
                reference_20=module!='Alternatives',scored=scored,
                source='Herbst et al. 2021 Fig. 1a; doi:10.1038/s41594-021-00682-7' if module!='Alternatives' else 'Yayli et al. 2023; doi:10.1016/j.celrep.2023.113099')
            for field in ['normZ','fdr_synth','fdr_supp','rank_synth']:
                record[field]=float(indexed.loc[source_gene,field]) if scored else None
            saga.append(record)
    saga=pd.DataFrame(saga)
    assert saga.reference_20.sum()==20
    saga.to_csv(OUT/'saga_reference_scores.csv',index=False)
    summary=dict(n_genes=len(dz),depleted_025=len(depleted),depleted_005=28,enriched_025=len(enriched),enriched_005=2,
        shared_method_genes=len(join),method_rank_spearman=float(rank_r),mageck_support_025=21,mageck_support_005=9,
        guide_log_cpm_reproduction_max_error=delta,guide_profile_genes=len(profiles.gene.unique()),
        guide_profile_guides=len(selected),zero_guides=selected.loc[selected.all_samples_zero,'sgRNA'].tolist(),
        endpoint_ymax=float(profiles.loc[profiles.timepoint.eq('Endpoint'),'log2_cpm_plus_1'].max()),
        d0_ymax=float(profiles.loc[profiles.timepoint.eq('D0'),'log2_cpm_plus_1'].max()),
        saga_reference_scored=int((saga.reference_20 & saga.scored).sum()),saga_additional_alternatives=SAGA['Alternatives'],
        saga_missing=saga.loc[~saga.scored,'gene'].tolist(),
        saga_reference_depleted_005=int((saga.reference_20 & saga.fdr_synth.lt(.05)).sum()),
        ppi_input_genes=len(ppi_input),ppi_input_method='MAGeCK depletion FDR <0.25',
        ppi_input_overlap_drugz=int(ppi_input.drugz_candidate.sum()),ppi_nodes=len(nodes),ppi_edges=len(edges),
        ppi_strong_string_physical_edges=int(edges.string_physical_score.ge(.7).sum()),
        gametolog_guides=len(gm),gametolog_guide_sequences=16,
        statement='Descriptive display of saved analysis; no hit/FDR refitting or D0 adjustment.')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    manifest={k:{'path':str(v),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()} for k,v in SOURCES.items()}
    manifest['preparation_script']={'path':str(Path(__file__).resolve()),'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
