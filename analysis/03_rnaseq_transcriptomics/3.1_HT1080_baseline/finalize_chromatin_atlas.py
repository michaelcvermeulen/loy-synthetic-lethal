"""Assemble reference features and prespecified gene sets; no RNA enrichment."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'outputs/tables/chromatin_atlas'
PROV = ROOT/'provenance/chromatin_atlas'
MARKS = ['H2AK119ub1','H3K27me3','H3K9me3','H4K20me3','H4K20me1',
         'H3K4me1','H3K27ac','H3K4me3','LaminB1','LaminB2']


def definitions(frame, suffix='__detected', fraction_suffix='__peak_covered_fraction'):
    def b(mark):
        return frame[mark+suffix].astype('boolean')
    result = {f'{m}_PROMOTER':b(m) for m in MARKS}
    result.update({
        'H3K4me3_H3K27me3_COMARKED_PROMOTER':b('H3K4me3') & b('H3K27me3'),
        'H3K27me3_H2AK119ub1_JOINT_PROMOTER':b('H3K27me3') & b('H2AK119ub1'),
        'H3K4me3_H3K27ac_JOINT_PROMOTER':b('H3K4me3') & b('H3K27ac'),
        'H3K9me3_H4K20me3_JOINT_PROMOTER':b('H3K9me3') & b('H4K20me3'),
        'LaminB1_LaminB2_JOINT_PROMOTER':b('LaminB1') & b('LaminB2'),
        'H3K27me3_DETECTED_H3K4me3_NOT_DETECTED':b('H3K27me3') & ~b('H3K4me3'),
        'H3K4me3_DETECTED_H3K27me3_NOT_DETECTED':b('H3K4me3') & ~b('H3K27me3'),
        'H3K9me3_DETECTED_H3K27me3_NOT_DETECTED':b('H3K9me3') & ~b('H3K27me3'),
    })
    return result


def peak_wide(frame, index='gene_id'):
    columns = ['detected','detected_relaxed_q','peak_covered_bp','peak_covered_fraction','n_merged_peak_segments','evaluable']
    result = frame.pivot(index=index,columns='mark',values=columns)
    result.columns = [f'{mark}__{feature}' for feature,mark in result.columns]
    return result


def main():
    config = json.loads((ROOT/'config/chromatin_atlas.json').read_text())
    genes = pd.read_csv(OUT/'genes_with_primary_TSS.csv').set_index('gene_id',drop=False)
    long = pd.read_csv(OUT/'promoter_peak_features_long.csv.gz')
    body = pd.read_csv(OUT/'gene_body_peak_features_long.csv.gz')
    primary_long = long[long.window.eq(config['primary_promoter_half_width'])].copy()
    primary = genes.join(peak_wide(primary_long))
    quality = primary_long.drop_duplicates('gene_id').set_index('gene_id')[
        ['start0','end0','unmasked_bp','unmasked_fraction','usable']].add_prefix('promoter_')
    primary = primary.join(quality)
    primary['reference_status'] = np.where(primary.promoter_usable.astype('boolean').fillna(False),
        'evaluable_under_frozen_rules',np.where(primary.promoter_usable.isna(),
        'chromosome_not_in_reference_tracks','more_than_20_percent_promoter_masked'))
    body_wide = peak_wide(body).add_prefix('body__')
    primary = primary.join(body_wide)
    # Signal is auxiliary; binary sets derive exclusively from reference peaks and frozen annotation.
    for mark in MARKS:
        sig = pd.read_csv(OUT/f'signal_{mark}.csv.gz')
        for region_type in ['promoter','gene_body']:
            cols = ['mean_signal_observed_bases','mean_signal_missing_as_zero','represented_fraction']
            piece = sig[sig.region_type.eq(region_type)].set_index('gene_id')[cols]
            piece.columns = [f'{region_type}__{mark}__{c}' for c in cols]
            primary = primary.join(piece)
    primary.to_csv(OUT/'HT1080_chromatin_features_all_genes.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    ranked = primary[primary.in_RNA_rank_universe].copy()
    assert len(ranked)==25461 and ranked.gene_name.is_unique
    ranked.to_csv(OUT/'HT1080_chromatin_features_RNA_universe.csv',index=False)
    selected_cols = ['gene_id','gene_name','transcript_id','selection_rule','chrom','strand','tss0',
                     'reference_status','promoter_unmasked_fraction']
    selected_cols += [m+'__detected' for m in MARKS]
    selected_cols += [m+'__peak_covered_fraction' for m in MARKS]
    ranked[selected_cols].to_csv(OUT/'HT1080_promoter_marks_RNA_universe.csv',index=False)

    membership = []
    for window, sub in long.groupby('window',sort=True):
        w = peak_wide(sub)
        for scenario,suffix in [(f'window_{window}','__detected')]+(
            [('relaxed_q_primary_window','__detected_relaxed_q')] if window==config['primary_promoter_half_width'] else []):
            for name,flag in definitions(w,suffix).items():
                for gene_id in w.index[flag.fillna(False)]:
                    membership.append((scenario,name,gene_id))
    # Alternative TSS sensitivity preserves co-marking within the SAME window.
    tss = pd.read_csv(OUT/'all_TSS_peak_features.csv.gz')
    for name,flag in definitions(tss).items():
        for gene_id in tss.loc[flag.fillna(False),'gene_id'].unique():
            membership.append(('any_annotated_TSS',name,gene_id))
    # High-coverage gene-body categories are explicit operational summaries, not inferred HMM domains.
    body_sets = {}
    for mark in ['H3K27me3','H3K9me3']:
        name = f'{mark}_GENE_BODY_PEAK_COVERAGE_GE50PCT'
        flag = primary[f'body__{mark}__peak_covered_fraction'].ge(0.5)
        body_sets[name] = flag
        for gene_id in primary.index[flag]:
            membership.append((f"window_{config['primary_promoter_half_width']}",name,gene_id))
    members = pd.DataFrame(membership,columns=['scenario','set_name','gene_id']).drop_duplicates()
    members = members.merge(genes[['gene_name','stable_gene_id','in_RNA_rank_universe']].reset_index(),on='gene_id',validate='m:1')
    members.to_csv(OUT/'gene_set_membership_all_scenarios.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    main_scenario = f"window_{config['primary_promoter_half_width']}"
    main_members = members[members.scenario.eq(main_scenario) & members.in_RNA_rank_universe]
    main_members.to_csv(OUT/'gene_set_membership_RNA_universe.csv',index=False)
    all_names = list(definitions(primary))+list(body_sets)
    summaries = []
    for name in all_names:
        chosen = main_members[main_members.set_name.eq(name)]
        n = len(chosen)
        summaries.append({'set_name':name,'n_RNA_universe':n,
                          'n_all_annotated_genes':len(members[members.scenario.eq(main_scenario)&members.set_name.eq(name)]),
                          'passes_existing_GSEA_size_15_500':15<=n<=500,
                          'size_status': 'below_15' if n<15 else 'above_500' if n>500 else '15_to_500',
                          'enrichment_tested':False})
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT/'gene_set_catalogue.csv',index=False)
    for field,label in [('gene_name','symbols'),('gene_id','gene_ids'),('stable_gene_id','stable_gene_ids')]:
        with (OUT/f'HT1080_GSE261834_measured_{label}.gmt').open('w',encoding='utf-8',newline='\n') as f:
            for name in all_names:
                values = sorted(main_members.loc[main_members.set_name.eq(name),field].unique())
                if values:
                    f.write('\t'.join(['HT1080_NEEDSEQ_'+name,'GSE261834;measured_reference;primary_TSS_2000bp;no_RNA_direction_selection']+values)+'\n')
    sensitivity = members[members.in_RNA_rank_universe].groupby(['scenario','set_name']).size().unstack(0,fill_value=0)
    sensitivity.to_csv(OUT/'gene_set_size_sensitivity.csv')
    # Mark coverage summaries include uncertainty and the two stricter overlap-fraction alternatives.
    coverage = []
    for mark in MARKS:
        flag = ranked[mark+'__detected'].astype('boolean')
        frac = ranked[mark+'__peak_covered_fraction']
        coverage.append({'mark':mark,'n_ranked_genes':len(ranked), 'evaluable':int(flag.notna().sum()),
                         'detected_any_bp':int(flag.sum()),'unknown':int(flag.isna().sum()),
                         'coverage_ge10pct':int(frac.ge(.1).sum()),'coverage_ge50pct':int(frac.ge(.5).sum())})
    pd.DataFrame(coverage).to_csv(OUT/'promoter_mark_summary.csv',index=False)
    candidates = ['CHRNB2','CHRNA7','CHRNA6','CHRNB4','CHRNE','MMP1','MMP3','ZMYND8',
                  'EZH2','EED','SUZ12','HDAC1','HDAC2','KDM5D','UTY','RCOR1']
    primary[primary.gene_name.isin(candidates)][['gene_id','gene_name','in_RNA_rank_universe','reference_status']+
        [m+'__detected' for m in MARKS]+[m+'__peak_covered_fraction' for m in MARKS]].to_csv(
            OUT/'selected_locus_reference_annotations.csv',index=False)
    unknown = ranked[ranked.reference_status.ne('evaluable_under_frozen_rules')]
    unknown[selected_cols].to_csv(OUT/'RNA_genes_with_uncertain_promoter_annotation.csv',index=False)
    result = {'ranked_genes':len(ranked),'all_annotated_genes':len(primary),
              'ranked_promoters_evaluable':int(ranked.promoter_usable.astype('boolean').fillna(False).sum()),
              'ranked_promoters_unknown':len(unknown),'set_definitions':len(summary),
              'nonempty_GMT_sets':int(summary.n_RNA_universe.gt(0).sum()),
              'sets_within_old_15_500_size_limits':int(summary.passes_existing_GSEA_size_15_500.sum()),
              'enrichment_run':False}
    (PROV/'atlas_summary.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)
    print(summary[['set_name','n_RNA_universe','size_status']].to_string(index=False))


if __name__=='__main__':
    main()
