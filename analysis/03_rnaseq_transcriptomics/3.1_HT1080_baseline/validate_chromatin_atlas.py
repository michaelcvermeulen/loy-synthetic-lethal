"""Scientific validation against source intervals and frozen resource hashes."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs/tables/chromatin_atlas'
PROV=ROOT/'provenance/chromatin_atlas'


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):
            h.update(b)
    return h.hexdigest()


def union(intervals):
    result=[]
    for a,b in sorted((int(a),int(b)) for a,b in intervals):
        if not result or a>result[-1][1]:
            result.append([int(a),int(b)])
        else:
            result[-1][1]=max(result[-1][1],int(b))
    return result


def bp_without_mask(intervals,mask):
    total=0
    for a,b in union(intervals):
        excluded=sum(max(0,min(b,y)-max(a,x)) for x,y in mask)
        total+=b-a-excluded
    return total


def main():
    config=json.loads((ROOT/'config/chromatin_atlas.json').read_text())
    checks=[]
    for r in pd.read_csv(PROV/'input_manifest.csv').itertuples():
        assert sha256(r.path)==r.sha256,r.path
    checks.append('All frozen source hashes unchanged')
    full=pd.read_csv(OUT/'HT1080_chromatin_features_all_genes.csv.gz',low_memory=False)
    ranked=pd.read_csv(OUT/'HT1080_chromatin_features_RNA_universe.csv')
    original=pd.read_csv(ROOT/config['rank_universe'],usecols=['gene_id','gene_name'])
    assert len(full)==78899 and full.gene_id.is_unique
    assert len(ranked)==25461 and ranked.gene_name.is_unique
    assert set(ranked.gene_id)==set(original.gene_id)
    assert dict(zip(ranked.gene_id,ranked.gene_name))==dict(zip(original.gene_id,original.gene_name))
    assert not any(x in ranked.columns for x in ['combined','mean_difference','AUG_2025','MAY_2026','pvalue','FDR'])
    checks.append('Full annotation and fixed ranked universe preserved; no RNA effects in atlas')
    tx=pd.read_csv(OUT/'transcripts.csv.gz')
    sel=tx[tx.is_selected_transcript].set_index('gene_id')
    assert len(sel)==len(full)
    expected=np.where(sel.strand.eq('+'),sel.tx_start0,sel.tx_end0-1)
    assert np.array_equal(sel.tss0,expected)
    assert dict(zip(full.gene_id,full.tss0))==sel.tss0.to_dict()
    assert (sel.tx_start0<=sel.tss0).all() and (sel.tss0<sel.tx_end0).all()
    checks.append('Plus/minus TSS coordinates verified against selected GTF transcript bounds')
    raw_mask=pd.read_csv(config['blacklist'],sep='\t',header=None,usecols=[0,1,2])
    masks={chrom:union(sub[[1,2]].to_numpy()) for chrom,sub in raw_mask.groupby(0)}
    peaks=pd.read_csv(OUT/'promoter_peak_features_long.csv.gz')
    checked=0
    for m in pd.read_csv(OUT/'mark_inputs.csv').itertuples():
        source=pd.read_csv(m.peak_path,sep='\t',header=None)
        source=source[source[8]>=m.min_minuslog10q]
        rows=peaks[peaks.mark.eq(m.mark)&peaks.window.eq(config['primary_promoter_half_width'])]
        # Deliberately include positive, negative, masked and both-strand intervals.
        subsets=[rows[rows.detected.eq(True)],rows[rows.detected.eq(False)],
                 rows[~rows.evaluable],rows[rows.strand.eq('-')]]
        sample=pd.concat([s.sample(min(12,len(s)),random_state=42) for s in subsets]).drop_duplicates('gene_id')
        for r in sample.itertuples():
            mask=masks.get(r.chrom,[])
            in_window=source[(source[0]==r.chrom)&(source[1]<r.end0)&(source[2]>r.start0)]
            clipped=[(max(r.start0,a),min(r.end0,b)) for a,b in in_window[[1,2]].itertuples(index=False,name=None)]
            bp=bp_without_mask(clipped,mask)
            denom=bp_without_mask([(r.start0,r.end0)],mask)
            assert bp==r.peak_covered_bp,(m.mark,r.gene_id,bp,r.peak_covered_bp)
            assert denom==r.unmasked_bp
            if r.evaluable:
                assert np.isclose(bp/denom,r.peak_covered_fraction)
                assert (bp>0)==bool(r.detected)
            else:
                assert pd.isna(r.detected) and pd.isna(r.peak_covered_fraction)
            checked+=1
    checks.append(f'{checked} peak/mask overlaps reproduced independently from raw peak files')
    columns=[c for c in ranked.columns if c.endswith('__peak_covered_fraction')]
    for c in columns:
        x=pd.to_numeric(ranked[c],errors='coerce').dropna()
        assert x.between(0,1).all(),c
    members=pd.read_csv(OUT/'gene_set_membership_RNA_universe.csv')
    catalogue=pd.read_csv(OUT/'gene_set_catalogue.csv')
    assert not members.duplicated(['set_name','gene_id']).any()
    for row in catalogue.itertuples():
        assert len(members[members.set_name.eq(row.set_name)])==row.n_RNA_universe
    for field,label in [('gene_name','symbols'),('gene_id','gene_ids'),('stable_gene_id','stable_gene_ids')]:
        for line in (OUT/f'HT1080_GSE261834_measured_{label}.gmt').read_text().splitlines():
            name,description,*genes=line.split('\t')
            assert len(genes)==len(set(genes))
            expected=set(members.loc[members.set_name.eq(name.removeprefix('HT1080_NEEDSEQ_')),field])
            assert set(genes)==expected
    checks.append('Coverage fractions valid; GMT membership matches gene IDs, symbols and catalogue')
    tss=pd.read_csv(OUT/'all_TSS_peak_features.csv.gz',usecols=['gene_id','H3K4me3__detected','H3K27me3__detected'])
    same_window=tss['H3K4me3__detected'].eq(True)&tss['H3K27me3__detected'].eq(True)
    expected=set(tss.loc[same_window,'gene_id'])
    all_members=pd.read_csv(OUT/'gene_set_membership_all_scenarios.csv.gz')
    observed=set(all_members.loc[all_members.scenario.eq('any_annotated_TSS')&
        all_members.set_name.eq('H3K4me3_H3K27me3_COMARKED_PROMOTER'),'gene_id'])
    assert observed==expected
    checks.append('Alternative-TSS co-marking requires both marks at the same TSS window')
    (PROV/'validation.json').write_text(json.dumps({'status':'passed','checks':checks},indent=2))
    print(json.dumps({'status':'passed','checks':checks},indent=2))


if __name__=='__main__':
    main()
