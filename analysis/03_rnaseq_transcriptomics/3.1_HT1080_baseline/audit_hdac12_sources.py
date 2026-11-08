"""Freeze public source documents for HDAC audit step 1; never refit expression models."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import re
import requests
import pandas as pd
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / 'resources/hdac12_reference_audit'
SETS = ([f'SENESE_{target}_TARGETS_{arm}' for target in ['HDAC1','HDAC2','HDAC1_AND_HDAC2','HDAC3'] for arm in ['UP','DN']]
        + ['ROPERO_HDAC2_TARGETS']
        + [f'HELLER_HDAC_TARGETS{extra}_{arm}' for extra in ['', '_SILENCED_BY_METHYLATION'] for arm in ['UP','DN']]
        + [f'{family}_{arm}' for family in ['PEART_HDAC_PROLIFERATION_CLUSTER','MARKS_HDAC_TARGETS'] for arm in ['UP','DN']]
        + ['GENES_CORRELATED_WITH_HDAC1_DELETION'])

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def fetch(item):
    name, url = item
    path = CACHE / name
    meta = path.with_suffix(path.suffix + '.retrieval.json')
    if path.exists() and meta.exists():
        old = json.loads(meta.read_text(encoding='utf-8'))
        if old.get('status') == 'retrieved':
            assert sha(path) == old['sha256'], f'Cached source changed: {path}'
            return old
    row = dict(file=str(path.relative_to(ROOT)).replace('\\','/'), url=url,
               retrieved_at=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    try:
        response = requests.get(url, timeout=(12, 35), headers={'User-Agent':'HDAC-reference-audit/1.0'}, allow_redirects=True)
        row.update(http_status=response.status_code, final_url=response.url,
                   content_type=response.headers.get('Content-Type',''))
        response.raise_for_status()
        response.encoding = response.apparent_encoding
        text = response.text[:10000] if 'pdf' not in row['content_type'] else ''
        if any(q in text for q in ['Checking your browser', 'reCAPTCHA', 'Cookies must be enabled']):
            row['status'] = 'browser_challenge_no_document'
        else:
            path.write_bytes(response.content)
            row.update(status='retrieved',sha256=sha(path),bytes=path.stat().st_size)
    except requests.RequestException as error:
        row.update(status='unavailable',error=str(error))
    meta.write_text(json.dumps(row,indent=2),encoding='utf-8')
    return row

def parse_definition(name):
    soup = BeautifulSoup((CACHE / (name + '.html')).read_bytes(), 'html.parser')
    count = re.search(r'show (\d+) source identifiers mapped to (\d+) genes', soup.get_text(' ',strip=True))
    assert count, f'Missing source mapping counts: {name}'
    heading = soup.find('th', title='original member source identifier')
    assert heading, f'Missing source identifier table: {name}'
    rows=[]
    for tr in heading.find_parent('table').find_all('tr'):
        cells=tr.find_all('td',recursive=False)
        if len(cells)==4:
            rows.append(dict(pathway=name,source_identifier=cells[0].get_text(' ',strip=True),
                             current_ncbi_id=cells[1].get_text(' ',strip=True),
                             current_symbol=cells[2].get_text(' ',strip=True)))
    frame=pd.DataFrame(rows)
    assert frame.source_identifier.nunique()==int(count[1]), name
    assert frame.loc[frame.current_symbol!='','current_symbol'].nunique()==int(count[2]), name
    namespace_row=next(tr for tr in soup.find_all('tr') if tr.find('th') and
                       'Source platform or' in tr.find('th').get_text(' ',strip=True))
    namespace=namespace_row.find('td').get_text(' ',strip=True)
    return frame, dict(source_identifiers=int(count[1]),live_unique_symbols=int(count[2]),source_namespace=namespace)

def study(name):
    if name.startswith('SENESE'): return 'Senese_2007'
    if name.startswith('HELLER'): return 'Heller_2008'
    if name.startswith('PEART'): return 'Peart_2005'
    if name.startswith('MARKS'): return 'Marks_2007_review'
    if name.startswith('ROPERO'): return 'Ropero_2008'
    return 'MSigDB_C9_dependency'

def build():
    # Check all pre-existing frozen resources before extending the audit cache.
    frozen=pd.read_csv(ROOT/'resources/manifest.csv',keep_default_na=False)
    assert all(sha(ROOT/r.path)==r.sha256 for r in frozen.itertuples()), 'Prior frozen resource changed'
    metadata=pd.read_csv(ROOT/'resources/msigdb_chromatin_descriptions.csv',keep_default_na=False)
    coverage=pd.read_csv(ROOT/'outputs/tables/pathway_resource_coverage.csv',keep_default_na=False).set_index('pathway')
    representatives=pd.read_csv(ROOT/'outputs/tables/symbol_mapping_representatives.csv',keep_default_na=False)
    symbol_ids=representatives.set_index('gene_name').gene_id.to_dict()
    membership=[]
    for chunk in pd.read_csv(ROOT/'resources/msigdb_all_membership.csv.gz',chunksize=150000,keep_default_na=False):
        membership.append(chunk[chunk.pathway.isin(SETS)])
    membership=pd.concat(membership,ignore_index=True).drop_duplicates(['pathway','symbol'])
    membership['represented_in_HT1080']=membership.symbol.isin(symbol_ids)
    membership['gene_id']=membership.symbol.map(symbol_ids).fillna('')
    mappings=[]; manifests=[]; drift=[]
    for name in SETS:
        mapping, counts=parse_definition(name); mappings.append(mapping)
        old=set(membership.loc[membership.pathway==name,'symbol'])
        live=set(mapping.current_symbol)-{''}
        for symbol in sorted(old^live):
            drift.append(dict(pathway=name,symbol=symbol,status='live_only' if symbol in live else 'frozen_only'))
        row=metadata.loc[metadata.gs_name==name].iloc[0].to_dict()
        row.update(pathway=name,study_id=study(name),**counts,
                   frozen_unique_symbols=len(old),represented_genes=len(old&set(symbol_ids)),
                   original_eligible=bool(coverage.loc[name,'eligible']),original_status=coverage.loc[name,'status'],
                   live_matches_frozen=old==live,
                   source_direction=('UP' if name=='ROPERO_HDAC2_TARGETS' else
                                     'UNSIGNED_ASSOCIATION' if name.startswith('GENES_') else name.rsplit('_',1)[1]),
                   definition_url='https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/'+name+'.html')
        row['evidence_class']=({'Senese_2007':'RNAi_expression_response', 'Ropero_2008':'cross_cell_line_mutation_association',
             'Heller_2008':'combined_drug_response' if 'SILENCED_BY_METHYLATION' in name else 'HDAC_inhibitor_response',
             'Peart_2005':'inhibitor_response_proliferation_subset', 'Marks_2007_review':'review_compilation',
             'MSigDB_C9_dependency':'cross_cell_line_dependency_association'})[row['study_id']]
        row['independence_unit']=row['study_id']
        targets=(name.removeprefix('SENESE_').split('_TARGETS')[0].replace('_AND_', ';').split(';')
                 if name.startswith('SENESE_') else ['HDAC2'] if name.startswith('ROPERO_') else [])
        row['perturbed_regulators_in_set']=';'.join(sorted(old & set(targets)))
        assert row['frozen_unique_symbols']==int(coverage.loc[name,'resource_size'])
        assert row['represented_genes']==int(coverage.loc[name,'overlap_size'])
        manifests.append(row)
    manifest=pd.DataFrame(manifests); mapping=pd.concat(mappings,ignore_index=True)
    output=ROOT/'outputs/tables'
    manifest.to_csv(output/'HDAC12_reference_manifest.csv',index=False)
    membership.to_csv(output/'HDAC12_frozen_membership.csv',index=False)
    mapping.to_csv(output/'HDAC12_source_identifier_mapping.csv',index=False)
    pd.DataFrame(drift,columns=['pathway','symbol','status']).to_csv(output/'HDAC12_live_vs_frozen_membership.csv',index=False)
    paper_counts={'HDAC1':987,'HDAC2':283,'HDAC1_AND_HDAC2':601,'HDAC3':1317}
    reconciled=[]
    for target,paper_count in paper_counts.items():
        a=manifest.set_index('pathway').loc['SENESE_'+target+'_TARGETS_UP']
        b=manifest.set_index('pathway').loc['SENESE_'+target+'_TARGETS_DN']
        names=['SENESE_'+target+'_TARGETS_'+arm for arm in ['UP','DN']]
        unique_probes=mapping.loc[mapping.pathway.isin(names),'source_identifier'].nunique()
        reconciled.append(dict(perturbation=target,paper_reported_entries=paper_count,
            MSigDB_UP_source_ids=int(a.source_identifiers),MSigDB_DN_source_ids=int(b.source_identifiers),
            MSigDB_source_id_sum=int(a.source_identifiers+b.source_identifiers),unique_source_ids_both_arms=unique_probes,
            source_sum_minus_paper=int(a.source_identifiers+b.source_identifiers-paper_count),
            status='count_reconciled_not_original_table_verified' if a.source_identifiers+b.source_identifiers==paper_count
                   else 'unresolved_count_discrepancy_original_table_unavailable'))
    pd.DataFrame(reconciled).to_csv(output/'HDAC12_Senese_count_reconciliation.csv',index=False)
    # This DOC stores its table as a contiguous ANSI run with Word cell delimiters.
    raw=(CACHE/'ropero_table1.doc').read_bytes(); start=raw.index(b'Gene accession number')
    ropero=[]
    for item in raw[start:].decode('latin1').split('\x07\x07')[1:]:
        cells=item.split('\x07')
        if len(cells)!=3: break
        accession,symbol,process=[v.strip() for v in cells]
        if not re.fullmatch(r'(?:gb:)?[A-Z]{1,3}_?\d+(?:\.\d+)?',accession): break
        ropero.append(dict(source_accession=accession,source_symbol=symbol,source_process_annotation=process))
    ropero=pd.DataFrame(ropero)
    ropero.to_csv(output/'HDAC12_Ropero_original_table.csv',index=False)
    declared=set(mapping.loc[mapping.pathway=='ROPERO_HDAC2_TARGETS','source_identifier'])
    observed=set(ropero.source_symbol)-{'---',''}
    assert len(ropero)>100, 'Incomplete Ropero table extraction'
    comparison=pd.DataFrame([dict(source_symbol=s,in_original_table=s in observed,in_MSigDB_source_ids=s in declared)
                             for s in sorted(declared|observed)])
    comparison.to_csv(output/'HDAC12_Ropero_table_reconciliation.csv',index=False)
    # Normalize retrieval status: a saved cookie challenge is not a recovered article.
    retrievals=[]
    for p in sorted(CACHE.glob('*.retrieval.json')):
        row=json.loads(p.read_text(encoding='utf-8'))
        source=ROOT/row['file']
        if row.get('status')=='retrieved':
            assert sha(source)==row['sha256'], f'Cached source changed: {source}'
            if source.suffix=='.html' and 'Cookies must be enabled' in source.read_text(errors='ignore'):
                row['status']='browser_challenge_no_document'
        retrievals.append(row)
    pd.DataFrame(retrievals).to_csv(output/'HDAC12_retrieval_manifest.csv',index=False)
    pd.DataFrame(STUDIES).to_csv(output/'HDAC12_source_studies.csv',index=False)
    audit=dict(sets=len(manifest),study_groups=manifest.study_id.nunique(),eligible_sets=int(manifest.original_eligible.sum()),
               mapping_rows=len(mapping),frozen_membership_rows=len(membership),live_frozen_differences=len(drift),
               ropero_original_rows=len(ropero),ropero_original_symbols=len(observed),ropero_curated_source_symbols=len(declared),
               ropero_rows_without_named_symbol=int(ropero.source_symbol.isin(['---','']).sum()),
               ropero_source_symmetric_difference=sorted(observed^declared),
               previous_resource_hashes_verified=len(frozen),new_expression_or_enrichment_tests_run=False)
    assert audit['sets']==18 and audit['study_groups']==6
    assert observed==declared, 'Ropero original named symbols differ from curated source IDs'
    (ROOT/'provenance/HDAC12_audit_checks.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    # Separate namespace records the exact audit inputs without changing original resources.
    paths=[p for p in CACHE.iterdir() if p.is_file() and p.name!='manifest.csv']
    pd.DataFrame([dict(path=str(p.relative_to(ROOT)).replace('\\','/'),sha256=sha(p),bytes=p.stat().st_size)
                  for p in sorted(paths)]).to_csv(CACHE/'manifest.csv',index=False)
    print(json.dumps(audit,indent=2))
    print(manifest[['pathway','source_identifiers','frozen_unique_symbols','represented_genes','original_eligible']].to_string(index=False))

STUDIES=[
 dict(study_id='Senese_2007',pmid='17470557',experiment='HDAC1;HDAC2;HDAC1+2;HDAC3 RNAi',cell_context='U2OS expression profiling',
      control='antiluciferase siRNA',timing='general protocol: second transfection at 24 h; harvest at 72 h; array-specific timing not independently verified',
      platform='Affymetrix HG-U133 Plus; MSigDB source namespace AFFY_HG_U133',
      selection_rule='Table 1S U/D flags; exact thresholds and multiplicity unavailable without supplement',
      replication='two chips per labeled target are technical duplicates; number of independent array RNA preparations unresolved',
      validation='independent RNA/qPCR for selected HDAC1 and double-knockdown genes; no equivalent validation of the full HDAC2 list identified',
      evidence_limit='one source study; protein cross-compensation; HDAC3 also reduced HDAC1/2 proteins; indirect responses included',
      verification='full primary paper and all eight MSigDB source-ID tables; original supplement unavailable',
      source_url='https://pubmed.ncbi.nlm.nih.gov/17470557/'),
 dict(study_id='Ropero_2008',pmid='18264134',experiment='HDAC2-mutant versus wild-type colorectal cell lines',
      cell_context='mutants RKO/Co115 versus HCT116/SW48/HCT15/LoVo',control='four different wild-type cell lines; not isogenic controls',
      timing='baseline mutation-associated profiles; not acute exposure',platform='expression microarray; exact platform not verified',
      selection_rule='Supplementary Table 1: Mann-Whitney P<0.01 and mean fold change>2; common UP genes in both mutant lines; adjustment not specified',
      replication='two mutant and four wild-type cell lines; within-line replicate counts unverified',
      validation='paper abstract reports selected promoter ChIP and wild-type HDAC2 re-expression',
      evidence_limit='cell-line genomic background confounding; targeted follow-up does not make the discovery list isogenic',
      verification='primary abstract and original Supplementary Table 1 recovered; full main text unavailable',
      source_url='https://www.nature.com/articles/onc200831'),
 dict(study_id='Heller_2008',pmid='18172295',experiment='TSA alone; separate Aza-dC plus TSA arms',
      cell_context='MM1;U266;NCI-H929 multiple-myeloma cells',control='untreated corresponding cell line',
      timing='not verified in accessible source methods',platform='Affymetrix HG-U133A',
      selection_rule='union across at least one of three cell lines; exact selection thresholds unverified',
      replication='indexed primary Results describe independently treated replicate RNA; per-condition count unverified',
      validation='selected induced genes investigated by methylation assays; not isoform-specific validation',
      evidence_limit='broad inhibitor; combination arms also perturb DNA methylation; lists are not common responses in all three lines',
      verification='primary abstract and indexed Results plus four MSigDB definitions; full methods blocked',
      source_url='https://aacrjournals.org/cancerres/article-abstract/68/1/44/540556/Genome-Wide-Transcriptional-Response-to-5-Aza-2'),
 dict(study_id='Peart_2005',pmid='15738394',experiment='SAHA and depsipeptide time-course; proliferation-selected subsets',
      cell_context='CEM cells',control='time-zero expression reference in reported clustering',timing='up to 16 h; not a single uniform endpoint',
      platform='microarray; source identifiers are sequence accessions',
      selection_rule='figure 2 requires at least 1.5-fold change at one later time; Table 3S proliferation subset; full statistical cutoff unverified',
      replication='figure 6 uses same microarray RNA for qPCR and reports at least three experiments; array replication structure unverified',
      validation='selected qPCR and APAF1 functional experiments; not genome-wide independent RNA validation',
      evidence_limit='preselected proliferation biology; broad pharmacology; two compounds within one study',
      verification='primary abstract/figure captions plus two MSigDB definitions; full methods unavailable',
      source_url='https://pubmed.ncbi.nlm.nih.gov/15738394/'),
 dict(study_id='Marks_2007_review',pmid='17322921',experiment='review compilation of inhibitor-responsive genes',
      cell_context='multiple literature contexts',control='not one controlled experiment',timing='heterogeneous/not one exposure',
      platform='review Table 2; source gene symbols',selection_rule='Table 2 induced/repressed lists; no single statistical threshold',
      replication='not applicable to a review compilation',validation='underlying studies require separate tracing',
      evidence_limit='not independent perturbation replication; biological drug specificity not inferred',
      verification='PubMed publication type Review plus MSigDB source-table definitions',source_url='https://pubmed.ncbi.nlm.nih.gov/17322921/'),
 dict(study_id='MSigDB_C9_dependency',pmid='',experiment='association of baseline expression with HDAC1 CRISPR dependency',
      cell_context='1106 shared DepMap/CCLE cell lines, release 25Q3',control='cross-cell-line association; not treated/control RNA',
      timing='not applicable',platform='baseline expression plus CRISPR viability dependency',
      selection_rule='information-theoretic ranking; ssGSEA size optimization; 66/34 train/test; held-out association threshold 0.20',
      replication='cell-line train/test split; not replicated HDAC1-deletion transcription',validation='held-out dependency association',
      evidence_limit='unsigned dependency association; cannot label UP or DN after deletion',
      verification='MSigDB C9 full definition',source_url='https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GENES_CORRELATED_WITH_HDAC1_DELETION.html')]

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--fetch',action='store_true')
    args=parser.parse_args()
    CACHE.mkdir(parents=True,exist_ok=True)
    if args.fetch:
        jobs=[(name+'.html','https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/'+name+'.html') for name in SETS]
        jobs += [('senese_2007.pdf','https://air.unimi.it/retrieve/dfa8b98e-e302-748b-e053-3a05fe0a3a96/4784.pdf'),
                 ('ropero_2008.html','https://www.nature.com/articles/onc200831'),
                 ('ropero_table1.doc','https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fonc.2008.31/MediaObjects/41388_2008_BFonc200831_MOESM5_ESM.doc'),
                 ('senese_supplement.html','http://bio.ifom-ieo-campus.it/supplementary/MCB49407/'),
                 ('senese_supplement_https.html','https://bio.ifom-ieo-campus.it/supplementary/MCB49407/')]
        jobs += [(f'pubmed_{pmid}.xml',f'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}%20AND%20SRC:MED&format=json&resultType=core')
                 for pmid in ['17470557','18172295','15738394','17322921','18264134']]
        with ThreadPoolExecutor(max_workers=3) as pool:
            rows=list(pool.map(fetch,jobs))
        (CACHE/'retrieval_manifest.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
        for row in rows: print(row['status'], row['file'], row.get('http_status',''))
    else:
        build()

if __name__ == '__main__':
    main()
