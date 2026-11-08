"""Freeze source hashes and choose RNA-independent GENCODE transcript TSSs."""
from pathlib import Path
from datetime import datetime, timezone
import gzip
import hashlib
import json
import re
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
OUT = ROOT/'outputs/tables/chromatin_atlas'
PROV = ROOT/'provenance/chromatin_atlas'


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(2**20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    PROV.mkdir(parents=True, exist_ok=True)
    config_path = ROOT/'config/chromatin_atlas.json'
    config = json.loads(config_path.read_text())
    inputs = [config_path, REPO/config['gtf'], ROOT/config['gene_annotation'],
              ROOT/config['rank_universe'], Path(config['blacklist']),
              ROOT/'provenance/chipseq_resources/geo_local_sample_metadata.csv']
    marks = []
    for mark in config['marks']:
        row = dict(mark)
        row['peak_path'] = str(Path(config['data_dir'])/(mark['stem']+mark['peak_suffix']))
        row['bigwig_path'] = str(Path(config['data_dir'])/(mark['stem']+'.bw'))
        row['min_minuslog10q'] = config[f"{mark['format']}_min_minuslog10q"]
        row['sensitivity_min_minuslog10q'] = config[f"sensitivity_{mark['format']}_min_minuslog10q"]
        marks.append(row)
        inputs.extend([Path(row['peak_path']),Path(row['bigwig_path'])])
    manifest = pd.DataFrame([{'path':str(p),'bytes':p.stat().st_size,'sha256':sha256(p)} for p in inputs])
    frozen = PROV/'input_manifest.csv'
    if frozen.exists():
        old = pd.read_csv(frozen)
        assert old.equals(manifest), 'Frozen inputs changed; review before replacing the atlas manifest.'
    else:
        manifest.to_csv(frozen,index=False)
        (PROV/'frozen_decisions.json').write_text(json.dumps({
            'frozen_utc':datetime.now(timezone.utc).isoformat(),
            'config_sha256':sha256(config_path), 'RNA_effect_values_used_for_selection':False,
            'config':config},indent=2),encoding='utf-8')
    pd.DataFrame(marks).to_csv(OUT/'mark_inputs.csv',index=False)
    genes = pd.read_csv(ROOT/config['gene_annotation'])
    rank = pd.read_csv(ROOT/config['rank_universe'],usecols=['gene_id','gene_name'])
    assert len(rank)==25461 and rank.gene_id.is_unique and rank.gene_name.is_unique
    assert rank.gene_id.isin(genes.gene_id).all()
    genes['in_RNA_rank_universe'] = genes.gene_id.isin(rank.gene_id)
    expected = dict(zip(genes.gene_id,genes.gene_name))
    assert all(expected[g]==n for g,n in rank.itertuples(index=False,name=None))
    print('Frozen inputs; extracting GENCODE transcripts.',flush=True)
    rows = []
    pattern = re.compile(r'(\w+) "([^"]*)";')
    with gzip.open(REPO/config['gtf'],'rt',encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or '\ttranscript\t' not in line:
                continue
            fields = line.rstrip('\n').split('\t')
            attrs = pattern.findall(fields[8])
            a = dict(attrs)
            tags = {v for k,v in attrs if k=='tag'}
            start,end = int(fields[3]),int(fields[4])
            # GTF is 1-based closed, so the minus-strand TSS base is end-1 in BED.
            tss = start-1 if fields[6]=='+' else end-1
            priority = next((i for i,t in enumerate(config['transcript_priority'][:3]) if t in tags),3)
            rows.append({'gene_id':a['gene_id'],'transcript_id':a['transcript_id'],
                         'chrom':fields[0],'strand':fields[6],'tx_start0':start-1,'tx_end0':end,
                         'tss0':tss,'span_bp':end-start+1,'priority':priority,
                         'tags':';'.join(sorted(tags))})
    tx = pd.DataFrame(rows)
    tx = tx[tx.gene_id.isin(genes.gene_id)].copy()
    assert tx.transcript_id.is_unique
    tx = tx.sort_values(['gene_id','priority','span_bp','transcript_id'],ascending=[True,True,False,True])
    selected = tx.drop_duplicates('gene_id').copy()
    selected['selection_rule'] = selected.priority.map(dict(enumerate(config['transcript_priority'][:4])))
    primary = genes.merge(selected[['gene_id','transcript_id','tss0','selection_rule']],on='gene_id',how='left',validate='1:1')
    primary['tss_missing'] = primary.tss0.isna()
    assert not primary.loc[primary.in_RNA_rank_universe,'tss_missing'].any()
    primary.to_csv(OUT/'genes_with_primary_TSS.csv',index=False)
    tx['is_selected_transcript'] = tx.transcript_id.isin(selected.transcript_id)
    tx.to_csv(OUT/'transcripts.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    tss = tx.groupby(['gene_id','chrom','strand','tss0'],sort=False).agg(
        transcript_ids=('transcript_id',lambda x:';'.join(x)),
        transcript_count=('transcript_id','size'),primary_TSS=('is_selected_transcript','any')).reset_index()
    tss['tss_id'] = tss.gene_id+':'+tss.tss0.astype(str)
    assert tss.tss_id.is_unique
    tss.to_csv(OUT/'distinct_TSS.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    summary = {'genes':len(genes),'ranked_genes':len(rank),'transcripts':len(tx),'distinct_TSS':len(tss),
               'missing_TSS':int(primary.tss_missing.sum()),
               'selection_rules':primary.selection_rule.value_counts().to_dict()}
    (PROV/'annotation_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary),flush=True)


if __name__=='__main__':
    main()
