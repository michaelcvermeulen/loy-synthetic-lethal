"""Read local bigWigs under WSL; exact interval integrals, no imputation."""
from pathlib import Path
import argparse
import gzip
import json
import sys
import numpy as np
import pandas as pd
import pyBigWig

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'outputs/tables/chromatin_atlas'
PROV = ROOT/'provenance/chromatin_atlas'


def local_path(text):
    text = str(text).replace('\\','/')
    if len(text)>2 and text[1]==':':
        return '/mnt/'+text[0].lower()+text[2:]
    return text


def integral_at(points, starts, ends, vals):
    lengths = ends-starts
    area = np.r_[0.,np.cumsum(lengths*vals)]
    coverage = np.r_[0.,np.cumsum(lengths)]
    n = np.searchsorted(ends,points,side='right')
    current = np.minimum(n,len(starts)-1)
    extra = np.clip(points-starts[current],0,lengths[current])
    extra[n==len(starts)] = 0
    return area[n]+extra*vals[current], coverage[n]+extra


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--metadata-only',action='store_true')
    parser.add_argument('--validate-only',action='store_true')
    args = parser.parse_args()
    marks = pd.read_csv(OUT/'mark_inputs.csv')
    headers = {}; chroms = {}
    for row in marks.itertuples():
        with pyBigWig.open(local_path(row.bigwig_path)) as bw:
            assert bw.isBigWig()
            sizes = bw.chroms()
            assert sizes['chr1']==248956422 and sizes['chrX']==156040895, 'Not expected GRCh38 lengths'
            for chrom,length in sizes.items():
                assert chrom not in chroms or chroms[chrom]==length, 'Conflicting sequence dictionaries'
                chroms[chrom]=length
            headers[row.mark] = {'chromosomes':sizes,'header':bw.header(),'pyBigWig':pyBigWig.__version__}
    (PROV/'bigwig_headers.json').write_text(json.dumps(headers,indent=2))
    pd.DataFrame(sorted(chroms.items()),columns=['chrom','length']).to_csv(OUT/'reference_chromosomes.csv',index=False)
    if args.metadata_only:
        print('Verified bigWig headers:',len(headers),flush=True)
        return
    regions = pd.read_csv(OUT/'signal_regions.csv')
    pieces = pd.read_csv(OUT/'signal_region_unmasked_pieces.csv.gz')
    if args.validate_only:
        # Independent basewise extraction includes interval boundaries, gaps and blacklist masks.
        possible = regions[(regions.end0-regions.start0<80000)&regions.unmasked_bp.gt(0)]
        sample = pd.concat([s.sample(min(15,len(s)),random_state=42)
                            for _,s in possible.groupby('region_type')])
        checks = []
        for row in marks.itertuples():
            saved = pd.read_csv(OUT/f'signal_{row.mark}.csv.gz').set_index('region_id')
            with pyBigWig.open(local_path(row.bigwig_path)) as bw:
                for r in sample.itertuples():
                    if r.chrom not in bw.chroms():
                        continue
                    vals = np.array(bw.values(r.chrom,int(r.start0),int(r.end0),numpy=True),dtype=float)
                    keep = np.zeros(len(vals),dtype=bool)
                    for piece in pieces[pieces.region_id.eq(r.region_id)].itertuples():
                        keep[int(piece.start0-r.start0):int(piece.end0-r.start0)] = True
                    valid = keep & np.isfinite(vals)
                    total = float(vals[valid].sum())
                    assert valid.sum()==saved.loc[r.region_id,'represented_bp']
                    assert np.isclose(total,saved.loc[r.region_id,'signal_integral'],rtol=1e-7,atol=1e-4)
                    checks.append({'mark':row.mark,'region_id':r.region_id,'basewise_integral':total,
                                   'stored_integral':saved.loc[r.region_id,'signal_integral']})
        pd.DataFrame(checks).to_csv(OUT/'signal_independent_validation.csv',index=False)
        (PROV/'signal_validation.json').write_text(json.dumps({'status':'passed',
            'regions_compared':len(checks),'method':'Independent basewise bigWig.values over unmasked intervals'},indent=2))
        print('Independent basewise signal checks passed:',len(checks),flush=True)
        return
    ids = {r:i for i,r in enumerate(regions.region_id)}
    pieces['row_index'] = pieces.region_id.map(ids)
    assert pieces.row_index.notna().all()
    record_counts=[]
    for row in marks.itertuples():
        sums = np.zeros(len(regions)); covered = np.zeros(len(regions))
        represented = np.zeros(len(regions),dtype=bool)
        with pyBigWig.open(local_path(row.bigwig_path)) as bw:
            for chrom,sub in pieces.groupby('chrom',sort=True):
                if chrom not in bw.chroms():
                    continue
                represented[regions.chrom.eq(chrom).to_numpy()] = True
                intervals = bw.intervals(chrom)
                if not intervals:
                    continue
                a = np.asarray(intervals,dtype=np.float64)
                starts,ends,values = a.T
                assert np.isfinite(values).all() and (ends>starts).all()
                assert (starts[1:]>=ends[:-1]).all()
                lo,c_lo = integral_at(sub.start0.to_numpy(),starts,ends,values)
                hi,c_hi = integral_at(sub.end0.to_numpy(),starts,ends,values)
                np.add.at(sums,sub.row_index.to_numpy(),hi-lo)
                np.add.at(covered,sub.row_index.to_numpy(),c_hi-c_lo)
                record_counts.append({'mark':row.mark,'chrom':chrom,'intervals':len(a),
                                      'median_interval_width':float(np.median(ends-starts)),
                                      'minimum_signal':float(values.min()),'maximum_signal':float(values.max())})
        denom = regions.unmasked_bp.to_numpy(dtype=float)
        result = regions[['region_id','gene_id','region_type']].copy()
        result['mark'] = row.mark
        result['chromosome_in_track'] = represented
        result['unmasked_bp'] = denom
        result['represented_bp'] = covered
        result['represented_fraction'] = np.divide(covered,denom,out=np.full_like(denom,np.nan),where=denom>0)
        result['mean_signal_observed_bases'] = np.divide(sums,covered,out=np.full_like(sums,np.nan),where=covered>0)
        result['mean_signal_missing_as_zero'] = np.divide(sums,denom,out=np.full_like(sums,np.nan),where=(denom>0)&represented)
        result['signal_integral'] = np.where(represented,sums,np.nan)
        assert (covered<=denom+1e-6).all()
        result.to_csv(OUT/f'signal_{row.mark}.csv.gz',index=False,compression={'method':'gzip','mtime':0})
        print(f'{row.mark}: {len(result)} signal regions complete',flush=True)
    pd.DataFrame(record_counts).to_csv(OUT/'bigwig_interval_qc.csv',index=False)


if __name__=='__main__':
    main()
