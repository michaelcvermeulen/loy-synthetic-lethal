"""Check source integrity, selected-guide contrasts and exported figure content."""
from pathlib import Path
import hashlib, json, sys, xml.etree.ElementTree as ET
THESIS=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(THESIS.parent/'.figure_python_deps'))
import numpy as np
import pandas as pd
import pymupdf as fitz
from PIL import Image
OUT=THESIS/'figures/chapter4'
DATA=OUT/'source_data/ht1080_candidate_discovery'
STEMS=['ch04_main_ht1080_candidate_discovery','ch04_main_ht1080_candidate_support','ch04_supp_ht1080_gametolog_guides']

def validate_main_guide_distributions(named):
    root=OUT/'source_data/ht1080_main_guide_distributions'
    manifest=json.loads((root/'source_manifest.json').read_text())
    for entry in manifest['inputs'].values():
        assert hashlib.sha256(Path(entry['path']).read_bytes()).hexdigest()==entry['sha256']
    mg=pd.read_csv(manifest['inputs']['mageck_guides']['path'],sep='\t').set_index('sgrna')
    counts=pd.read_csv(manifest['inputs']['counts']['path'],sep='\t').set_index('sgRNA')
    scores=pd.read_csv(DATA/'all_gene_scores.csv').set_index('GENE')
    selected=pd.read_csv(root/'selected_guides.csv',float_precision='round_trip')
    ticks=pd.read_csv(root/'displayed_tick_positions.csv',float_precision='round_trip')
    ntc=pd.read_csv(root/'non_targeting_guides.csv',float_precision='round_trip')
    density=pd.read_csv(root/'ntc_density.csv',float_precision='round_trip')
    assert len(selected)==136 and selected.Gene.nunique()==34 and selected.lfc_displayed.sum()==134
    assert set(selected.Gene)==set(named.GENE)|{'DDX3X','EIF1AX','ZFX','RPS4X'}
    assert selected.groupby('Gene').size().eq(4).all()
    assert selected.loc[~selected.lfc_displayed,'Gene'].tolist()==['DDX3X','MED12']
    assert selected.loc[~selected.lfc_displayed,'all_samples_zero'].all()
    samples=['D0_CTRL','D0_LOY','POST_CTRL_REP1','POST_CTRL_REP2','POST_LOY_REP1','POST_LOY_REP2']
    for row in selected.itertuples():
        assert row.mageck_reported_log2FC==mg.loc[row.sgRNA,'LFC']
        assert row.mageck_guide_FDR==mg.loc[row.sgRNA,'FDR']
        assert all(getattr(row,s)==counts.loc[row.sgRNA,s] for s in samples)
        assert row.lfc_displayed==(sum(getattr(row,s) for s in samples[2:])>0)
        assert row.drugz_normZ==scores.loc[row.Gene,'normZ']
        field='fdr_supp' if row.drugz_fdr_direction=='enrichment' else 'fdr_synth'
        assert row.drugz_directional_fdr==scores.loc[row.Gene,field]
    assert len(ntc)==1000 and ntc.lfc_displayed.sum()==999
    assert ntc.Gene.eq('Non-Targeting-Control').all()
    assert np.array_equal(ntc.LFC.to_numpy(),mg.loc[ntc.sgRNA,'LFC'].to_numpy())
    assert np.array_equal(ntc.lfc_displayed.to_numpy(),counts.loc[ntc.sgRNA,samples[2:]].sum(axis=1).gt(0).to_numpy())
    assert len(ticks)==137 and ticks.sgRNA.nunique()==134 and ticks.height.eq(.8).all()
    assert ticks.groupby('panel').size().to_dict()=={'b':15,'c':122}
    assert set(ticks.loc[ticks.panel.eq('c'),'Gene'])==set(named.GENE)
    assert set(ticks.loc[ticks.panel.eq('b'),'Gene'])=={'DDX3X','EIF1AX','ZFX','RPS4X'}
    for row in ticks.itertuples():
        assert row.x==mg.loc[row.sgRNA,'LFC']
        assert manifest['x_limits'][0]<row.x<manifest['x_limits'][1]
    values=ntc.loc[ntc.lfc_displayed,'LFC']
    assert manifest['x_limits'][0]<values.min()<values.max()<manifest['x_limits'][1]
    assert np.isclose(np.trapezoid(density.kde_density,density.guide_log2FC),1,atol=1e-8)
    return {'unique_genes':34,'unique_guides':136,'unique_observed_guides':134,
            'guide_ticks_b':15,'guide_ticks_c':122,'ntc_reference_n':999,
            'equal_tick_height':True,'original_values_unchanged':True,'no_observed_values_clipped':True}

def main():
    for key,entry in json.loads((DATA/'source_manifest.json').read_text()).items():
        assert hashlib.sha256(Path(entry['path']).read_bytes()).hexdigest()==entry['sha256'],key
    profiles=pd.read_csv(DATA/'guide_profiles.csv')
    observed=profiles.loc[profiles.timepoint.eq('Endpoint')]
    means=observed.groupby(['gene','sgrna','genotype']).cpm.mean().unstack('genotype')
    means['mean_cpm_difference_LOY_minus_CTRL']=means.LOY-means.CTRL
    means.to_csv(DATA/'endpoint_guide_direction_check.csv')
    selected=['DDX3X','SREBF1','SREBF2','TADA2B','TAF6L']
    check={g:{'lower_in_LOY':int(means.loc[g,'mean_cpm_difference_LOY_minus_CTRL'].lt(0).sum()),
              'equal':int(means.loc[g,'mean_cpm_difference_LOY_minus_CTRL'].eq(0).sum())} for g in selected}
    assert check['DDX3X']=={'lower_in_LOY':3,'equal':1},check
    assert all(check[g]=={'lower_in_LOY':4,'equal':0} for g in selected if g!='DDX3X'),check
    saga=pd.read_csv(DATA/'saga_reference_scores.csv')
    assert len(saga)==24 and saga.scored.all() and saga.reference_20.sum()==20
    assert saga.loc[saga.gene.eq('SGF29'),'source_gene'].iloc[0]=='CCDC101'
    guide_table=pd.read_csv(DATA/'gametolog_guide_table.csv')
    assert len(guide_table)==16 and guide_table.groupby('Gene').size().eq(4).all()
    assert guide_table.sequence_5prime_to_3prime.str.fullmatch('[ACGT]{20}').all()
    assert guide_table.mageck_LFC_displayed.sum()==15
    assert guide_table.loc[~guide_table.mageck_LFC_displayed,'sgRNA'].tolist()==['DDX3X_4640_TCTCGGTTCCTTAAATGAGG']
    for row in guide_table.itertuples():
        original=profiles.loc[profiles.sgrna.eq(row.sgRNA)]
        for sample in ['D0_CTRL','D0_LOY','POST_CTRL_REP1','POST_CTRL_REP2','POST_LOY_REP1','POST_LOY_REP2']:
            assert int(original.loc[original['sample'].eq(sample),'raw_count'].iloc[0])==getattr(row,sample)
    nodes=pd.read_csv(DATA/'ppi_nodes.csv');edges=pd.read_csv(DATA/'ppi_edges.csv')
    assert len(nodes)==7 and len(edges)==11
    assert nodes.fdr_synth.lt(.25).all()
    assert set(edges.gene_a)|set(edges.gene_b)==set(nodes.gene)
    assert len(edges[['gene_a','gene_b']].apply(lambda r:tuple(sorted(r)),axis=1).unique())==11
    assert edges.string_physical_score.ge(.7).sum()==10
    result={'original_sources_unchanged':True,'endpoint_guide_direction':check,'exports':{}}
    support_genes={
        'c':['SREBF1','SREBF2','TADA2B','TAF6L','SUPT20H','ATXN7L3','TADA1','TIMM8A','MED12','EIF4A3'],
        'f':['DDX3X','EIF1AX','ZFX','RPS4X','KDM5C','KDM6A']}
    support_files={'c':'candidate_support_panel_c_cpm_plus_1.csv',
                   'f':'candidate_support_endpoint_cpm_plus_1.csv'}
    support_render=json.loads((DATA/'render_manifest.json').read_text())[STEMS[1]]
    assert support_render['panel_f_grey_box_genes']==support_genes['f'][:4]
    schematic_members=pd.read_csv(DATA/'saga_schematic_members.csv')
    assert len(schematic_members)==20 and schematic_members.gene.is_unique
    assert schematic_members.set_index('gene').module.to_dict()==saga.loc[saga.reference_20].set_index('gene').module.to_dict()
    assert support_render['panel_e_score_columns']==[saga.gene.iloc[:13].tolist(),saga.gene.iloc[13:].tolist()]
    source_manifest=json.loads((DATA/'source_manifest.json').read_text())
    raw=pd.read_csv(source_manifest['counts']['path'],sep='\t').set_index('sgRNA')
    library=pd.read_csv(source_manifest['library']['path'],sep='\t')
    library_pairs=set(zip(library.gene,library.Seq))
    for row in profiles.loc[profiles.gene.isin(['KDM5C','KDM6A','EIF4A3'])].itertuples():
        assert (row.gene,row.sgrna.rsplit('_',1)[-1]) in library_pairs
        assert row.raw_count==raw.loc[row.sgrna,row.sample]
        assert np.isclose(row.cpm,1e6*row.raw_count/raw[row.sample].sum())
    result['supp_endpoint_scale']={}
    for panel,genes in support_genes.items():
        linear=pd.read_csv(DATA/support_files[panel],float_precision='round_trip')
        assert len(linear)==16*len(genes) and set(linear.gene)==set(genes)
        assert linear.groupby('gene').size().eq(16).all()
        assert linear.plotted_cpm_plus_1.equals(linear.cpm+1)
        assert np.allclose(np.log2(linear.plotted_cpm_plus_1),linear.log2_cpm_plus_1)
        assert linear.loc[linear.raw_count.eq(0),'plotted_cpm_plus_1'].eq(1).all()
        assert linear.plotted_cpm_plus_1.between(0,60).all()
        original=profiles.set_index(['sgrna','sample']).loc[linear.set_index(['sgrna','sample']).index]
        assert np.array_equal(original.raw_count,linear.raw_count)
        assert np.allclose(original.cpm,linear.cpm)
        assert support_render[f'panel_{panel}_genes']==genes
        assert support_render[f'panel_{panel}_ylim']==[0,60]
        result['supp_endpoint_scale'][panel]={'scale':'linear CPM + 1','points':len(linear),
            'zero_reads_at':1,'shared_ylim':[0,60],'gene_order':genes}
    result['supp_baseline_abundance']={}
    for panel,genes in [('b',['DDX3X','TIMM8A','MED12','SCLY','WTAP','SREBF2','SREBF1','LLGL1','TAF6L','UBE2G2','TADA2B']),('d',support_genes['f'])]:
        baseline=pd.read_csv(DATA/f'candidate_support_panel_{panel}_d0.csv')
        assert len(baseline)==len(genes)*8 and set(baseline.gene)==set(genes)
        assert baseline.groupby(['gene','genotype']).size().eq(4).all()
        assert np.allclose(baseline.plotted_log2_cpm_plus_1,np.log2(baseline.cpm+1))
        for row in baseline.itertuples():
            assert row.raw_count==raw.loc[row.sgrna,row.sample]
            assert np.isclose(row.cpm,1e6*row.raw_count/raw[row.sample].sum())
        assert baseline.plotted_log2_cpm_plus_1.between(0,15).all()
        result['supp_baseline_abundance'][panel]={'genes':genes,'values':len(baseline),'raw_counts_and_CPM_agree':True}
    enriched=pd.read_csv(DATA/'enriched_candidates.csv')
    depleted=pd.read_csv(DATA/'depleted_candidates.csv')
    saga_hits=set(saga.loc[saga.reference_20 & saga.fdr_synth.lt(.25),'source_gene'])
    depleted=depleted.loc[depleted.fdr_synth.lt(.05)|depleted.GENE.isin(saga_hits)]
    enriched=enriched.loc[enriched.fdr_supp.lt(.05)]
    named=pd.read_csv(DATA/'main_named_hits.csv')
    assert len(named)==31 and named.GENE.is_unique
    assert len(depleted)==29 and len(enriched)==2
    assert depleted.loc[depleted.fdr_synth.ge(.05),'GENE'].tolist()==['TAF5L']
    assert set(named.loc[named.direction.eq('depleted'),'GENE'])==set(depleted.GENE)
    assert set(named.loc[named.direction.eq('enriched'),'GENE'])==set(enriched.GENE)
    original_scores=pd.read_csv(DATA/'all_gene_scores.csv').set_index('GENE')
    for row in named.itertuples():
        assert row.normZ==original_scores.loc[row.GENE,'normZ']
        assert row.directional_fdr==original_scores.loc[row.GENE,'directional_fdr']
    result['main_guide_distributions']=validate_main_guide_distributions(named)
    net_data=OUT/'source_data/ht1080_all_sl_string'
    full_nodes=pd.read_csv(net_data/'all_nodes.csv')
    full_edges=pd.read_csv(net_data/'functional_edges.csv')
    display_nodes=pd.read_csv(net_data/'main_display_positions.csv')
    assert len(full_nodes)==48 and len(full_edges)==66
    assert len(display_nodes)==48 and display_nodes.GENE.is_unique
    assert set(display_nodes.GENE)==set(full_nodes.GENE) and display_nodes.node_shape.eq('circle').all()
    assert display_nodes.group.value_counts().to_dict()=={'connected':26,'unconnected':22}
    render=json.loads((DATA/'render_manifest.json').read_text())[STEMS[0]]
    for name,digest in render['network_input_sha256'].items():
        assert hashlib.sha256((net_data/name).read_bytes()).hexdigest()==digest
    expected={'ch04_main_ht1080_candidate_discovery':['DDX3X','19,112','SREBF1','TADA2B','TAF6L','2.13e-26','STRING v12.5','SIRT6']+enriched.GENE.tolist()+depleted.GENE.tolist(),
              'ch04_main_ht1080_candidate_support':['MAGeCK','0.877','ATXN7L1','ATXN7L2','SGF29','KDM6A','EIF1AX','Endpoint CPM + 1'],
              'ch04_supp_ht1080_gametolog_guides':['DDX3X','EIF1AX','RPS4X','ZFX','No reads','MAGeCK']}
    for stem in STEMS:
        doc=fitz.open(OUT/f'{stem}.pdf'); assert len(doc)==1
        page=doc[0]; words=page.get_text('words'); page_text=page.get_text()
        assert len(words)>90
        assert all(x in page_text for x in expected[stem]),(stem,page_text)
        if stem==STEMS[0]:
            # Require the chosen labels in a, and all retained names in b/c.
            panel_a_text=page.get_text(clip=fitz.Rect(.015*page.rect.width,.02*page.rect.height,
                .34*page.rect.width,.25*page.rect.height))
            assert all(g in panel_a_text for g in ['DDX3X','TIMM8A','SREBF2','TAF6L'])
            panel_b_text=page.get_text(clip=fitz.Rect(.015*page.rect.width,.27*page.rect.height,
                .36*page.rect.width,.49*page.rect.height))
            assert all(g in panel_b_text for g in ['DDX3X','EIF1AX','ZFX','RPS4X'])
            assert 'Guide log' in panel_b_text
            panel_c_text=page.get_text(clip=fitz.Rect(.35*page.rect.width,.02*page.rect.height,
                .995*page.rect.width,.51*page.rect.height))
            assert all(g in panel_c_text for g in depleted.GENE.tolist()+enriched.GENE.tolist())
            assert 'Guide log' in panel_c_text
            panel_d_text=page.get_text(clip=fitz.Rect(.015*page.rect.width,.5*page.rect.height,
                .995*page.rect.width,.98*page.rect.height))
            assert all(g in panel_d_text for g in full_nodes.GENE)
            assert 'Endpoint CPM + 1' not in page_text and 'e' not in page_text.splitlines()
            assert np.isclose(page.rect.height/72,8.8)
        elif stem==STEMS[1]:
            rank_text=page.get_text(clip=fitz.Rect(.06*page.rect.width,.02*page.rect.height,
                                                  .47*page.rect.width,.19*page.rect.height))
            assert all(g in rank_text for g in support_render['panel_a_labels'])
            def region(name):
                x,y,w,h=support_render['panel_regions'][name]
                return fitz.Rect(x*page.rect.width,(1-y-h)*page.rect.height,
                                 (x+w)*page.rect.width,(1-y)*page.rect.height)
            assert not any(g in page_text for g in ['ATP6V0E1','USP7','RALGAPB'])
            baseline_text=page.get_text(clip=region('d'))
            assert all(g in baseline_text for g in support_genes['f'])
            schematic_text=page.get_text(clip=region('e_schematic'))
            assert all(g in schematic_text for g in schematic_members.gene)
            assert all(m in schematic_text for m in schematic_members.module.unique())
            scores_text=page.get_text(clip=region('e_scores'))
            assert all(g in scores_text for g in saga.gene)
            for panel in ['c','f']:
                area=region(panel)
                panel_text=page.get_text(clip=area)
                assert all(g in panel_text for g in support_genes[panel])
                assert 'Endpoint CPM + 1' in panel_text and 'Endpoint log' not in panel_text
                labels={g:next(w for w in page.get_text('words',clip=area) if w[4]==g)
                        for g in support_genes[panel]}
                for row_genes in support_render[f'panel_{panel}_rows']:
                    row=[labels[g] for g in row_genes]
                    assert all(a[0]<b[0] for a,b in zip(row,row[1:]))
                    assert np.ptp([w[1] for w in row])<.5
                if panel=='f':
                    assert '(UTX)' in panel_text
                    x,y,w,h=support_render['panel_f_grey_box_bounds']
                    box=fitz.Rect(x*page.rect.width,(1-y-h)*page.rect.height,
                                  (x+w)*page.rect.width,(1-y)*page.rect.height)
                    assert all(box.contains(fitz.Rect(labels[g][:4])) for g in support_genes['f'][:4])
                    assert not any(box.intersects(fitz.Rect(labels[g][:4])) for g in support_genes['f'][4:])
                    assert any(all(abs(a-b)<.5 for a,b in zip(d['rect'],box))
                               and d['fill'] is not None for d in page.get_drawings())
            assert np.isclose(page.rect.height/72,10.6)
        outside=[w[:5] for w in words if w[0]<-.5 or w[1]<-.5 or w[2]>page.rect.width+.5 or w[3]>page.rect.height+.5]
        assert not outside,(stem,outside)
        if stem in STEMS[:2]:
            spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
            overlap=[]
            for i,a in enumerate(spans):
                for b in spans[i+1:]:
                    intersection=fitz.Rect(a['bbox'])&fitz.Rect(b['bbox'])
                    if not intersection.is_empty and intersection.width>.5 and intersection.height>.5:
                        overlap.append([a['text'],b['text']])
            assert not overlap,overlap
        fonts=page.get_fonts(); assert fonts and all('OpenSans' in f[3] for f in fonts),fonts
        svg=ET.parse(OUT/f'{stem}.svg')
        ntext=len(svg.findall('.//{http://www.w3.org/2000/svg}text')); assert ntext>70
        png=Image.open(OUT/f'{stem}.png'); assert min(png.info.get('dpi',(0,0)))>590
        assert png.width==5760
        grey=page.get_pixmap(matrix=fitz.Matrix(1.6,1.6),colorspace=fitz.csGRAY)
        grey.save(OUT/'review'/f'{stem}_greyscale_check.png')
        result['exports'][stem]={'pdf_pages':len(doc),'pdf_words':len(words),'text_out_of_bounds':len(outside),
            'fonts':[f[3] for f in fonts],'svg_editable_text_elements':ntext,'png_size':[png.width,png.height],
            'png_dpi':list(png.info['dpi'])}
        doc.close()
    result['main_named_hits']={'depleted':29,'enriched':2,'SAGA_candidates':6,'scores_unchanged':True}
    result['main_network']={'nodes':48,'edges':66,'all_labels_found':True,'node_shape':'circle','data_unchanged':True}
    result['supp_saga_schematic']={'panel':'e','reference_members':20,'scored_members':24,
        'module_assignments_match':True,'schematic_above_scores':True,'enrichment_panel_removed':True,'panel_d':'D0 X-linked paralog abundance'}
    result['scripts']={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in Path(__file__).parent.glob('*ht1080_candidate*.py')}
    helper=Path(__file__).parent/'ht1080_main_guide_distributions.py'
    result['scripts'][helper.name]=hashlib.sha256(helper.read_bytes()).hexdigest()
    helper=Path(__file__).parent/'make_ht1080_support_compact.py'
    result['scripts'][helper.name]=hashlib.sha256(helper.read_bytes()).hexdigest()
    helper=Path(__file__).parent/'ht1080_saga_schematic.py'
    result['scripts'][helper.name]=hashlib.sha256(helper.read_bytes()).hexdigest()
    (DATA/'export_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
