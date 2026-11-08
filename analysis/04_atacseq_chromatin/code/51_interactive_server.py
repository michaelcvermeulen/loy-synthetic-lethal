"""Read-only localhost viewer with on-demand coverage from the accepted ATAC fragments."""
from pathlib import Path
from collections import OrderedDict
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs, unquote
import argparse
import json
import math
import threading
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = ['CTRL_11','CTRL_12','CTRL_21','LOY_11','LOY_12','LOY_21']

class LocusData:
    def __init__(self):
        self.genes=pd.read_csv(ROOT/'integration/gene_atlas/data/gene_models.csv',dtype={'gene_id':str}).set_index('gene_id')
        self.exons=pd.read_csv(ROOT/'integration/gene_atlas/data/merged_gene_exons.csv',dtype={'gene_id':str})
        self.peaks=pd.read_csv(ROOT/'annotation/annotated_differential_peaks.csv').set_index('peak_id')
        links=pd.read_csv(ROOT/'integration/feature_overlaps/peak_gene_link_audit.csv',dtype={'entrez_id':str})
        links=links[['feature','entrez_id','peak_id']].drop_duplicates()
        self.links={(g,f):list(v.peak_id) for (g,f),v in links.groupby(['entrez_id','feature'])}
        self.counts=pd.read_csv(ROOT/'counts/peak_counts.tsv',sep='\t').set_index('peak_id')
        norm=pd.read_csv(ROOT/'differential/normalization_factors.csv')
        norm=norm[norm.arm=='primary_TMM'].set_index('sample').loc[SAMPLES,'effective_size']
        self.cpm=self.counts[SAMPLES].div(norm,axis=1)*1e6
        self.npz={s:np.load(ROOT/f'counts/{s}_fragments.npz') for s in SAMPLES}
        self.denominator={s:json.loads((ROOT/f'qc/{s}_preparation.json').read_text())['nonY_fragments'] for s in SAMPLES}
        self.sizes={k:int(v) for k,v in (line.split() for line in (ROOT/'config/bam.chrom.sizes').read_text().splitlines())}
        self.cache=OrderedDict();self.responses=OrderedDict();self.lock=threading.RLock()

    def coverage(self,chrom,left,right):
        width=max(25,math.ceil((right-left)/1500/25)*25)
        starts=np.arange(left,right,width);ends=np.minimum(starts+width,right)
        with self.lock:
            if chrom not in self.cache:
                self.cache[chrom]={s:(self.npz[s][chrom+'_start'],self.npz[s][chrom+'_end']) for s in SAMPLES}
                while len(self.cache)>2:self.cache.popitem(last=False)
            self.cache.move_to_end(chrom)
            matrix=[]
            for sample in SAMPLES:
                a,b=self.cache[chrom][sample]
                stop=np.searchsorted(a,right,side='left')
                ix=np.flatnonzero(b[:stop]>left)
                a=np.maximum(a[ix].astype(np.int64),left)-left
                b=np.minimum(b[ix].astype(np.int64),right)-left
                n=len(starts);sb=a//width;eb=b//width
                delta=np.bincount(sb,minlength=n+1)-np.bincount(eb,minlength=n+1)
                integral=np.cumsum(delta)[:n].astype(float)*width
                integral-=np.bincount(sb,weights=a%width,minlength=n+1)[:n]
                integral+=np.bincount(eb,weights=b%width,minlength=n+1)[:n]
                assert np.isclose(integral.sum(),(b-a).sum()) and np.all(integral>=0)
                matrix.append((integral/(ends-starts)*1e6/self.denominator[sample]).round(7).tolist())
        return starts.tolist(),ends.tolist(),matrix,width

    def locus(self,gene_id,feature,view):
        if (gene_id,feature) not in self.links or gene_id not in self.genes.index:
            raise ValueError('Gene is not linked to tested peaks in this feature.')
        key=(gene_id,feature,view)
        with self.lock:
            if key in self.responses:return self.responses[key]
        gene=self.genes.loc[gene_id];chrom=gene.chrom
        linked=self.peaks.loc[self.links[key[:2]]].sort_values(['FDR','start'])
        if view=='overview':
            left=min(int(gene.start)-1,int(linked.start.min()))-10000
            right=max(int(gene.end),int(linked.end.max()))+10000
            selected_peak=None
        else:
            selected_peak=linked.index[0] if view=='best' else view
            if selected_peak not in linked.index:raise ValueError('Peak is not linked to this gene and feature.')
            p=linked.loc[selected_peak];left=int(p.start)-5000;right=int(p.end)+5000
        left=max(0,left//25*25);right=min(self.sizes[chrom],math.ceil(right/25)*25)
        starts,ends,coverage,width=self.coverage(chrom,left,right)
        local=self.peaks[(self.peaks.chrom==chrom)&(self.peaks.start<right)&(self.peaks.end>left)]
        models=self.genes[(self.genes.chrom==chrom)&(self.genes.start-1<right)&(self.genes.end>left)]
        exon=self.exons[(self.exons.chrom==chrom)&(self.exons.start0<right)&(self.exons.end>left)]
        model_rows=[]
        for gid,row in models.iterrows():
            e=exon[exon.gene_id==gid]
            model_rows.append({'id':gid,'symbol':row.symbol,'start':int(row.start)-1,'end':int(row.end),'strand':row.strand,
                               'exons':e[['start0','end']].values.tolist()})
        linked_rows=[]
        for pid,p in linked.iterrows():
            linked_rows.append({'id':pid,'start':int(p.start),'end':int(p.end),'logFC':p.logFC,'FDR':p.FDR,
                                'raw':[int(x) for x in self.counts.loc[pid,SAMPLES]],'cpm':self.cpm.loc[pid,SAMPLES].round(7).tolist()})
        linked_set=set(linked.index)
        result={'gene':gene_id,'feature':feature,'chrom':chrom,'start':left,'end':right,'bin':width,'samples':SAMPLES,
                'starts':starts,'ends':ends,'coverage':coverage,'selected_peak':selected_peak or linked.index[0],
                'linked_peaks':linked_rows,'genes':model_rows,
                'peaks':[{'id':pid,'start':int(p.start),'end':int(p.end),'logFC':p.logFC,'FDR':p.FDR,'linked':pid in linked_set} for pid,p in local.iterrows()]}
        with self.lock:
            self.responses[key]=result
            while len(self.responses)>60:self.responses.popitem(last=False)
        return result

class Handler(SimpleHTTPRequestHandler):
    def __init__(self,*args,**kwargs):super().__init__(*args,directory=str(ROOT),**kwargs)
    def log_message(self,fmt,*args):pass
    def send_json(self,obj,status=200):
        encoded=json.dumps(obj,separators=(',',':'),allow_nan=False).encode()
        self.send_response(status);self.send_header('Content-Type','application/json');self.send_header('Cache-Control','no-store')
        self.send_header('Content-Length',str(len(encoded)))
        if self.headers.get('Origin')=='null':self.send_header('Access-Control-Allow-Origin','null')
        self.end_headers();self.wfile.write(encoded)
    def do_GET(self):
        host=self.headers.get('Host','')
        if host not in (f'127.0.0.1:{self.server.server_port}',f'localhost:{self.server.server_port}'):
            self.send_error(403);return
        origin=self.headers.get('Origin')
        if origin and origin not in ('null',f'http://127.0.0.1:{self.server.server_port}',f'http://localhost:{self.server.server_port}'):
            self.send_error(403);return
        path=urlparse(self.path)
        if path.path=='/api/health':self.send_json({'status':'ready','application':'HT1080_ATAC_RNA'});return
        if path.path=='/api/locus':
            q=parse_qs(path.query)
            try:self.send_json(self.server.data.locus(q.get('gene',[''])[0],q.get('feature',['promoter'])[0],q.get('view',['best'])[0]))
            except ValueError as e:self.send_json({'error':str(e)},400)
            except Exception as e:
                print('Locus request failed:',repr(e),flush=True);self.send_json({'error':'The locus could not be loaded. See the viewer log.'},500)
            return
        if path.path=='/':
            self.send_response(302);self.send_header('Location','/report/interactive_ATAC_RNA.html');self.end_headers();return
        target=(ROOT/unquote(path.path).lstrip('/')).resolve()
        if not target.is_relative_to(ROOT) or not target.is_file() or target.parts[len(ROOT.parts)] not in ('report','integration'):
            self.send_error(404);return
        super().do_GET()

def main():
    p=argparse.ArgumentParser();p.add_argument('--port',type=int,default=8765);args=p.parse_args()
    server=ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    server.data=LocusData()
    print(f'Viewer ready: http://127.0.0.1:{args.port}/report/interactive_ATAC_RNA.html',flush=True)
    server.serve_forever()

if __name__=='__main__':main()
