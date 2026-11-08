"""Freeze HGNC ribosomal protein groups and map the exclusion to existing GSEA genes."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import hashlib
import json
import pandas as pd
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'integration/pathway_GSEA/ribosomal_sensitivity'
(OUT/'sources').mkdir(parents=True,exist_ok=True)
GROUPS={728:('cytosolic','S ribosomal proteins'),729:('cytosolic','L ribosomal proteins'),646:('mitochondrial','Mitochondrial ribosomal proteins')}


def fetch(gid):
    path=OUT/f'sources/HGNC_group_{gid}.json'
    url=f'https://rest.genenames.org/fetch/gene_group_id/{gid}'
    if path.exists():data=json.loads(path.read_text())
    else:
        response=requests.get(url,headers={'Accept':'application/json'},timeout=30);response.raise_for_status();data=response.json()
        assert data['response']['numFound']>0
        path.write_text(json.dumps(data,indent=2),encoding='utf-8')
    return gid,data,url,path


records=[];manifest=[]
for gid,data,url,path in ThreadPoolExecutor(3).map(fetch,GROUPS):
    for d in data['response']['docs']:
        assert gid in d['gene_group_id'] and d['status']=='Approved'
        records.append({'group_id':gid,'class':GROUPS[gid][0],'group_name':GROUPS[gid][1],'hgnc_id':d['hgnc_id'],'symbol':d['symbol'],'entrez_id':d.get('entrez_id',''),'gene_name':d['name'],'locus_type':d['locus_type']})
    manifest.append({'url':url,'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'retrieved_UTC':datetime.now(timezone.utc).isoformat()})
catalog=pd.DataFrame(records)
catalog.to_csv(OUT/'sources/HGNC_ribosomal_protein_catalog.csv',index=False)
(OUT/'sources/source_manifest.json').write_text(json.dumps(manifest,indent=2))
ranks=pd.read_csv(ROOT/'integration/pathway_GSEA/gene_ranks.csv',dtype={'entrez_id':str})
mapped=[]
for r in ranks.itertuples():
    match=catalog[(catalog.entrez_id==r.entrez_id)|(catalog.symbol==r.symbol)]
    if match.empty:continue
    assert match['class'].nunique()==1
    exact_id=bool((match.entrez_id==r.entrez_id).any())
    exact_symbol=bool((match.symbol==r.symbol).any())
    mapped.append({'entrez_id':r.entrez_id,'source_symbol':r.symbol,'class':match['class'].iloc[0],
        'hgnc_symbol':';'.join(sorted(set(match.symbol))),'hgnc_entrez_id':';'.join(sorted(set(match.entrez_id))),
        'match_by_id':exact_id,'match_by_symbol':exact_symbol,'mapping_note':'ID and symbol agree' if exact_id and exact_symbol else 'Source symbol matches; Entrez differs' if exact_symbol else 'Entrez matches; source symbol differs'})
mapped=pd.DataFrame(mapped);mapped.to_csv(OUT/'excluded_genes.csv',index=False)
mapped[~(mapped.match_by_id&mapped.match_by_symbol)].to_csv(OUT/'exclusion_mapping_discrepancies.csv',index=False)
prefix=ranks[ranks.symbol.str.match(r'^(RPL|RPS|MRPL|MRPS)',na=False)]
prefix[~prefix.entrez_id.isin(mapped.entrez_id)].to_csv(OUT/'prefix_matches_retained.csv',index=False)
print(mapped.groupby('class').size().to_string())
print(mapped[~(mapped.match_by_id&mapped.match_by_symbol)].to_string(index=False))
