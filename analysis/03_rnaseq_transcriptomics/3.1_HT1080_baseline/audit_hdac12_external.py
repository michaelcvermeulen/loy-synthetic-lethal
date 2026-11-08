"""Define auditable within-source contrasts before external response modeling."""
from pathlib import Path
import sys
import json
import gzip
import csv
import hashlib
import pandas as pd

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / 'resources/hdac12_external/GSE193885'
sys.path.insert(0, str(ROOT / 'tools/hdac12_python'))
from pyxlsb import open_workbook

def main():
    m = pd.read_csv(ROOT / 'outputs/tables/HDACEXT_GEO_sample_metadata.csv', keep_default_na=False)
    m['sample_key'] = m.title.str.replace('HAP1 cell line ', '', regex=False)
    fields = m.sample_key.str.extract(r'^(.*)_([123])_S(47|52)[0-9]+$')
    assert fields.notna().all().all()
    m['source_group'] = fields[0]
    m['replicate_label'] = fields[1].astype(int)
    m['sample_ID_stratum'] = 'S' + fields[2]
    m['model_group'] = m.sample_ID_stratum + '__' + m.source_group
    m['used_in_model'] = ~m.source_group.isin(['HAP_MS1', 'HAP_MS3'])
    m['model_status'] = m.used_in_model.map({True: 'included', False: 'dose_time_label_unresolved'})
    m['replication_scope'] = 'author_reported_biological_triplicates_of_cell_line_or_treatment; independent_clones_not_demonstrated'
    m['stratum_scope'] = 'derived_sample_ID_group_not_verified_experimental_batch'
    m['transgene_species'] = 'not_applicable'
    m.loc[m.source_group.isin(['HD1_141','HD1_wt','HD2_142','HD2_wt','HAP_141','HAP_142']), 'transgene_species'] = 'mouse_cDNA_in_human_cells'
    m.loc[m.source_group.isin(['HD3_135_G10','HD3_wt_A9','HAP_135']), 'transgene_species'] = 'human_cDNA'
    m.loc[m.source_group.isin(['HD8_143','HD8_wt','HAP_143']), 'transgene_species'] = 'not_independently_verified'
    m['drug_metadata_scope'] = 'not_drug_treated'
    m.loc[m.source_group.eq('HAP_MS'), 'drug_metadata_scope'] = 'MS275_label_MS; generic_protocol_3uM_24h; untreated_HAP_control'
    m.loc[m.source_group.eq('HAP_JQ'), 'drug_metadata_scope'] = 'JQ12_generic_protocol_1uM_24h; untreated_HAP_control'
    m.loc[~m.used_in_model, 'drug_metadata_scope'] = 'MS275_label_MS1_or_MS3; distinct_dose_time_unresolved'
    assert m.groupby('model_group').size().eq(3).all()
    m.to_csv(ROOT / 'outputs/tables/HDACEXT_sample_manifest.csv', index=False)
    rows=[]
    def add(name, stratum, numerator, denominator, target, role, tier):
        num=stratum+'__'+numerator; den=stratum+'__'+denominator
        assert m[m.used_in_model & m.model_group.eq(num)].shape[0] == 3
        assert m[m.used_in_model & m.model_group.eq(den)].shape[0] == 3
        rows.append(dict(contrast=name, sample_ID_stratum=stratum, numerator=num, denominator=den,
            target=target, comparison_type=role, display_tier=tier, n_numerator=3, n_denominator=3,
            direction='numerator_minus_denominator', inference_scope='within_established_source_cell_lines_conditional_on_reported_culture_replicates'))
    for target, stratum, knockout, inactive, active, wt_inactive in [
        ('HDAC1','S47','HD1_KO','HD1_141','HD1_wt','HAP_141'),
        ('HDAC2','S47','HD2_KO','HD2_142','HD2_wt','HAP_142'),
        ('HDAC3','S52','HD3_KO','HD3_135_G10','HD3_wt_A9','HAP_135'),
        ('HDAC8','S47','HD8_KO','HD8_143','HD8_wt','HAP_143')]:
        add(target+'_KO_vs_WT',stratum,knockout,'HAP',target,'protein_loss','primary')
        add(target+'_CI_vs_active_in_KO',stratum,inactive,active,target,'inactive_vs_active_transgene','primary')
        add(target+'_CI_in_WT_vs_WT','S52',wt_inactive,'HAP',target,'inactive_transgene_in_WT','primary')
        add(target+'_CI_in_KO_vs_WT',stratum,inactive,'HAP',target,'inactive_transgene_in_KO_vs_WT','supporting')
        add(target+'_CI_in_KO_vs_KO',stratum,inactive,knockout,target,'inactive_protein_presence_in_KO','supporting')
        add(target+'_active_in_KO_vs_WT',stratum,active,'HAP',target,'active_transgene_control_vs_WT','supporting')
    add('MS275_MS_S47_vs_WT','S47','HAP_MS','HAP','MS275','inhibitor_exposure','primary')
    add('MS275_MS_S52_vs_WT','S52','HAP_MS','HAP','MS275','inhibitor_exposure','primary')
    add('JQ12_S52_vs_WT','S52','HAP_JQ','HAP','JQ12','inhibitor_exposure','primary')
    pd.DataFrame(rows).to_csv(ROOT / 'outputs/tables/HDACEXT_contrast_manifest.csv', index=False)
    skipped=[dict(comparison='MS275_'+group+'_S52_vs_WT', status='not_run', reason='Distinct dose/time meaning unresolved in generic GEO metadata and supplement') for group in ['MS1','MS3']]
    skipped.append(dict(comparison='direct_CI_KO_vs_CI_WT_for_HDAC1_HDAC2_HDAC8', status='not_run', reason='Different sample-ID strata; background effect cannot be separated from unverified stratum differences'))
    pd.DataFrame(skipped).to_csv(ROOT / 'outputs/tables/HDACEXT_skipped_comparisons.csv', index=False)
    output=ROOT/'outputs/tables/HDACEXT_published_DESeq2_results.csv.gz'
    # The original worksheet labels and effect statistics are retained for compatibility checks.
    with open_workbook(str(CACHE/'hess_2022_S3_RNAseq.xlsb')) as book, gzip.open(output,'wt',encoding='utf-8',newline='') as handle:
        writer=csv.writer(handle)
        writer.writerow(['source_worksheet','gene_id','gene_version','gene_name','gene_source','gene_biotype','location','baseMean','published_logFC','published_lfcSE','published_P','published_adjusted_P'])
        counts={}
        for name in book.sheets:
            if name=='Description': continue
            count=0
            with book.get_sheet(name) as sheet:
                for i,row in enumerate(sheet.rows()):
                    if i==0: continue
                    values=[cell.v for cell in row]
                    if values and isinstance(values[0],str) and values[0].startswith('ENSG'):
                        writer.writerow([name]+values[:11]);count+=1
            counts[name]=count
            print('Published sheet',name,count,flush=True)
    audit=dict(study='GSE193885',libraries=69,modeled_libraries=int(m.used_in_model.sum()),
        source_groups=int(m.model_group.nunique()),modeled_groups=int(m.loc[m.used_in_model,'model_group'].nunique()),
        contrasts=len(rows),primary_contrasts=sum(r['display_tier']=='primary' for r in rows),
        sample_ID_strata=m.sample_ID_stratum.value_counts().to_dict(),
        verified_source_fields=['GSM_to_count_column','genetic_background','active_or_inactive_transgene','reported_biological_triplicates'],
        unresolved_fields=['independently_generated_clone_replication','exact_culture_harvest_pairing','experimental_batch_dates','MS_MS1_and_MS3_dose_time','HDAC8_transgene_species'],
        supplemental_ENA_audit='Three experiment and three sample XML records repeat generic treatment metadata; they do not resolve MS/MS1/MS3 dose or time.',
        alignment_limit='Human GRCh38 alignment of mouse HDAC1/2 cDNAs does not distinguish endogenous and introduced transcripts; HDAC-family sensitivity is declared.',
        source_versions='GEO lists DESeq2 1.30.0/R4.0.0; paper lists DESeq2 1.34.0/R4.1.3. This analysis uses a declared independent edgeR fit.',
        published_worksheet_rows=counts,pyxlsb_version='1.0.10')
    (ROOT/'provenance/HDACEXT_source_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit,indent=2),flush=True)

if __name__=='__main__': main()
