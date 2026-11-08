"""Render co-dependency heatmaps and complete top-partner atlases."""
from pathlib import Path
import sys, os, json, math, argparse, hashlib

THESIS=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(THESIS.parent/'.figure_python_deps'))
OUT=THESIS/'figures/chapter4'
DATA=OUT/'source_data/ht1080_depmap_codependency'
os.environ.setdefault('MPLCONFIGDIR',str(DATA/'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
FONT=apply_style(base_size=9)
REVIEW=OUT/'review'
REVIEW.mkdir(exist_ok=True)
ATLAS=OUT/'depmap_codependency_partner_atlas'
ATLAS.mkdir(exist_ok=True)
A=pd.read_csv(DATA/'seed_annotations.csv').fillna('').set_index('gene')
R=pd.read_csv(DATA/'seed_pairwise_pearson.csv',index_col=0)
ADJ=pd.read_csv(DATA/'seed_pairwise_lineage_sex_residual_pearson.csv',index_col=0)
TOP=pd.read_csv(DATA/'top_partners.csv').fillna('')
MANIFEST=json.loads((DATA/'analysis_manifest.json').read_text())
FIGURE_ROLES={stem:'supplementary draft' for stem in [
    'ch04_supp_ht1080_depmap_hit_codependencies_draft',
    'ch04_supp_ht1080_depmap_saga_codependencies_draft',
    'ch04_supp_ht1080_depmap_hit_partner_atlas',
    'ch04_supp_ht1080_depmap_saga_partner_atlas']}
DIR={'Depleted':'#795082','Enriched':'#00877C','Not nominated':'#737B84'}
MOD={'Core':'#B49C28','HAT':'#3CA4AE','DUB':'#8A4D99','TRRAP':'#737373','Splicing':'#7896AA','Alternatives':'#CF8B4F'}
CMAP=plt.get_cmap('RdBu_r').copy(); CMAP.set_bad('#D8DCE0')


def average_link_order(corr):
    """UPGMA: minimum average 1-r between clusters, deterministic tie handling."""
    genes=list(corr.index); n=len(genes)
    distances={(i,j):float(1-corr.iloc[i,j]) for i in range(n) for j in range(i+1,n)}
    active={i:[i] for i in range(n)}; sizes={i:1 for i in range(n)}
    merges=[]
    while len(active)>1:
        a,b=min(((i,j) for i in active for j in active if i<j),key=lambda ij:(distances[ij],ij))
        k=n+len(merges); distance=distances[(a,b)]
        for c in active:
            if c in (a,b): continue
            da=distances[tuple(sorted((a,c)))]; db=distances[tuple(sorted((b,c)))]
            distances[tuple(sorted((c,k)))]=(sizes[a]*da+sizes[b]*db)/(sizes[a]+sizes[b])
        active[k]=active.pop(a)+active.pop(b)
        sizes[k]=sizes[a]+sizes[b]
        merges.append([a,b,distance,sizes[k]])
    order=[genes[i] for i in next(iter(active.values()))]
    assert len(set(order))==n and np.all(np.diff([m[2] for m in merges])>=-1e-10)
    return order,merges


def matrix(ax,table,genes,size=7.5,annotate=False):
    vals=table.reindex(index=genes,columns=genes).to_numpy()
    im=ax.imshow(vals,cmap=CMAP,vmin=-1,vmax=1,interpolation='none',aspect='equal')
    ax.set_xticks(range(len(genes)),genes,rotation=90,fontsize=size)
    ax.set_yticks(range(len(genes)),genes,fontsize=size)
    ax.tick_params(length=0,pad=3)
    for spine in ax.spines.values(): spine.set_visible(False)
    if annotate:
        for i in range(len(genes)):
            for j in range(i):
                v=vals[i,j]
                if np.isfinite(v) and abs(v)>=.4:
                    ax.text(j,i,f'{v:.2f}'.replace('0.','.'),ha='center',va='center',fontsize=6,
                            color='white' if abs(v)>.62 else '#222222')
    return im


def save(fig,stem):
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{stem}.{ext}',dpi=600)
    fig.savefig(REVIEW/f'{stem}_preview.png',dpi=160)
    plt.close(fig)


def draw_hit_dendrograms(top_ax,left_ax,genes,order,merges):
    """Draw the same linkage tree on both axes, aligned to heatmap cell centres."""
    nodes={i:(order.index(g),0.) for i,g in enumerate(genes)}
    for i,(left,right,height,size) in enumerate(merges):
        x1,h1=nodes[int(left)]; x2,h2=nodes[int(right)]
        leaf_position=[x1,x1,x2,x2]
        merge_height=[h1,height,height,h2]
        top_ax.plot(leaf_position,merge_height,color='#636B75',lw=.65)
        left_ax.plot(merge_height,leaf_position,color='#636B75',lw=.65)
        nodes[len(genes)+i]=((x1+x2)/2,height)
    limit=merges[-1][2]*1.06
    top_ax.set(xlim=(-.5,len(order)-.5),ylim=(0,limit),ylabel='1 − r')
    top_ax.set_xticks([]); top_ax.set_yticks([0,.5,1])
    top_ax.tick_params(axis='y',labelsize=7,length=2,pad=2)
    top_ax.yaxis.label.set_size(7.5)
    for spine in top_ax.spines.values(): spine.set_visible(False)
    left_ax.set(xlim=(limit,0),ylim=(len(order)-.5,-.5))
    left_ax.axis('off')


def hits_heatmap():
    hits=A.index[A.screen_hit.eq(True)].tolist()
    present=[g for g in hits if A.loc[g,'depmap_available']]
    order,merges=average_link_order(R.loc[present,present])
    missing=[g for g in hits if g not in present]
    order+=missing
    pd.DataFrame({'gene':order,'row':range(1,len(order)+1)}).to_csv(DATA/'hit_heatmap_order.csv',index=False)
    pd.DataFrame(merges,columns=['left_cluster','right_cluster','distance_1_minus_r','size']).to_csv(DATA/'hit_heatmap_linkage.csv',index=False)
    pd.DataFrame({'leaf_id':range(len(present)),'gene':present}).to_csv(DATA/'hit_heatmap_leaf_ids.csv',index=False)
    clustering={'method':'Average linkage (UPGMA)', 'distance':'1 - Pearson r',
                'clustered_genes':len(present),'input_order':present,'display_order':order,
                'missing_excluded_from_clustering':missing,'same_tree_for_rows_and_columns':True,
                'cluster_cut':'No discrete cluster cut or significance threshold imposed.',
                'source_sha256':hashlib.sha256((DATA/'seed_pairwise_pearson.csv').read_bytes()).hexdigest()}
    (DATA/'hit_clustering_manifest.json').write_text(json.dumps(clustering,indent=2),encoding='utf-8')
    fig=plt.figure(figsize=(10,10.7))
    x,y,w=.125,.16,.76
    h=w*10/10.7
    ax=fig.add_axes([x,y,w,h])
    im=matrix(ax,R,order,size=7.3)
    ax.yaxis.tick_right()
    ax.tick_params(axis='y',pad=10,labelright=True,labelleft=False)
    top=fig.add_axes([x,y+h+.008,w,.075])
    left=fig.add_axes([.022,y,.095,h])
    draw_hit_dendrograms(top,left,present,order,merges)
    for t in ax.get_xticklabels()+ax.get_yticklabels():
        g=t.get_text()
        t.set_color(DIR[A.loc[g,'screen_direction']])
        if A.loc[g,'screen_fdr']<.05: t.set_fontweight('semibold')
        if g in missing: t.set_color('#777777')
    for i,g in enumerate(order):
        if bool(A.loc[g,'saga_reference_20']):
            ax.add_patch(Rectangle((len(order)-.25,i-.37),.22,.74,facecolor='#BAA04B',edgecolor='none',clip_on=False))
    if missing:
        ax.axvline(len(present)-.5,color='white',lw=1.5)
        ax.axhline(len(present)-.5,color='white',lw=1.5)
    fig.text(.125,.978,'56 / 57 screen hits in DepMap · 1,186 models',fontsize=9)
    cb=fig.colorbar(im,cax=fig.add_axes([.67,.064,.25,.013]),orientation='horizontal',ticks=[-1,0,1])
    cb.set_label('Pearson r',fontsize=8,labelpad=1); cb.ax.tick_params(labelsize=7.5,length=2)
    handles=[Patch(facecolor=DIR['Depleted'],label='Depleted in LOY (48)'),
             Patch(facecolor=DIR['Enriched'],label='Enriched in LOY (9)'),
             Patch(facecolor='#BAA04B',label='SAGA reference member'),
             Patch(facecolor='#D8DCE0',label='Unavailable: ZCCHC16')]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.12,.008),ncol=2,fontsize=8,frameon=False)
    fig.text(.67,.009,'Semibold: screen FDR <0.05',fontsize=7.5,va='bottom')
    save(fig,'ch04_supp_ht1080_depmap_hit_codependencies_draft')
    return order


def saga_heatmaps():
    genes=A.index[A.saga_module.ne('')].tolist()
    # Keep the frozen structural grouping, then the four separately labelled alternatives.
    genes=[g for module in MOD for g in genes if A.loc[g,'saga_module']==module]
    pd.DataFrame({'gene':genes,'module':[A.loc[g,'saga_module'] for g in genes]}).to_csv(DATA/'saga_heatmap_order.csv',index=False)
    fig=plt.figure(figsize=(9.6,5.7))
    for x,table,letter,heading in [(.105,R,'a','Unadjusted'),(.60,ADJ,'b','Lineage + recorded sex residuals')]:
        ax=fig.add_axes([x,.265,.35,.35*9.6/5.7])
        im=matrix(ax,table,genes,size=7.2,annotate=True)
        for t in ax.get_xticklabels()+ax.get_yticklabels():
            g=t.get_text(); t.set_color(MOD[A.loc[g,'saga_module']])
            if bool(A.loc[g,'screen_hit']): t.set_fontweight('semibold')
        bounds=[]
        for i,g in enumerate(genes[1:],1):
            if A.loc[g,'saga_module']!=A.loc[genes[i-1],'saga_module']: bounds.append(i-.5)
        for v in bounds:
            ax.axhline(v,color='white',lw=1); ax.axvline(v,color='white',lw=1)
        ax.text(0,1.03,heading,transform=ax.transAxes,fontsize=8.5)
        fig.text(x-.075,.905,letter,fontsize=13,fontweight='semibold')
    fig.legend(handles=[Patch(facecolor=c,label=m if m!='Alternatives' else 'Alternative members') for m,c in MOD.items()],
               loc='lower left',bbox_to_anchor=(.075,.055),ncol=3,fontsize=7.7,frameon=False,columnspacing=1.5)
    cb=fig.colorbar(ScalarMappable(norm=Normalize(-1,1),cmap=CMAP),cax=fig.add_axes([.655,.1,.295,.018]),orientation='horizontal',ticks=[-1,0,1])
    cb.set_label('Pearson r',fontsize=8,labelpad=2); cb.ax.tick_params(labelsize=7.5,length=2)
    fig.text(.105,.017,'Semibold: six screen hits (FDR <0.25). Lower-triangle values shown for |r| ≥0.40.',fontsize=7.6)
    save(fig,'ch04_supp_ht1080_depmap_saga_codependencies_draft')


def partner_atlas(genes,stem):
    pages=math.ceil(len(genes)/8)
    with PdfPages(OUT/f'{stem}.pdf') as pdf:
        for page in range(pages):
            fig,axes=plt.subplots(4,2,figsize=(9.6,10.2))
            fig.subplots_adjust(left=.16,right=.965,bottom=.115,top=.90,hspace=.72,wspace=.70)
            for ax,g in zip(axes.flat,genes[page*8:(page+1)*8]):
                row=A.loc[g]
                if not row.depmap_available:
                    ax.axis('off'); ax.text(0,.85,g,fontweight='semibold',fontsize=10)
                    ax.text(0,.55,'Not present in the DepMap\ngene-effect matrix.',fontsize=9)
                    continue
                t=TOP.loc[TOP.seed.eq(g)&TOP.selection.eq('Top positive')].head(8)
                y=np.arange(len(t))
                colours=[DIR[v] for v in t.partner_screen_direction]
                ax.hlines(y,0,t.pearson_r,color='#C8CDD3',lw=1)
                ax.scatter(t.pearson_r,y,c=colours,s=22,zorder=3,edgecolor='white',linewidth=.35)
                ax.set_yticks(y,t.partner,fontsize=8)
                ax.set_xlim(0,1); ax.set_ylim(len(t)-.5,-.65)
                ax.set_xticks([0,.5,1]); ax.tick_params(axis='y',length=0,pad=4)
                ax.tick_params(axis='x',labelsize=7)
                for label in ax.get_yticklabels():
                    if label.get_text() in A.index and bool(A.loc[label.get_text(),'saga_reference_20']): label.set_fontweight('semibold')
                for yy,rr in zip(y,t.pearson_r): ax.text(rr+.025 if rr<.88 else rr-.025,yy,f'{rr:.2f}',va='center',ha='left' if rr<.88 else 'right',fontsize=7)
                suffix=f" · {row.saga_module}" if row.saga_module else ''
                ax.set_title(g+suffix,loc='left',fontweight='semibold',fontsize=9,pad=5)
                ax.set_xlabel('Pearson r',fontsize=8,labelpad=2)
            for ax in list(axes.flat)[len(genes[page*8:(page+1)*8]):]: ax.axis('off')
            fig.text(.16,.97,'Genome-wide positive co-dependency partners',fontsize=10,fontweight='semibold')
            fig.text(.16,.945,'DepMap 25Q3 · ≥200 shared models · top 8 per seed, without an r cutoff',fontsize=8)
            fig.legend(handles=[Line2D([0],[0],marker='o',ls='',color=c,label=d) for d,c in DIR.items()],
                       loc='lower left',bbox_to_anchor=(.15,.018),ncol=3,frameon=False,fontsize=8)
            fig.text(.965,.035,f'{page+1} / {pages}',ha='right',fontsize=8)
            fig.text(.16,.012,'Point colour: screen-hit direction. Semibold partner labels: SAGA reference members.',fontsize=7)
            pdf.savefig(fig)
            fig.savefig(ATLAS/f'{stem}_page{page+1:02}.svg')
            fig.savefig(ATLAS/f'{stem}_page{page+1:02}.png',dpi=160)
            plt.close(fig)
    return pages


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--hits-only',action='store_true',help='Redraw only the clustered all-hit heatmap.')
    args=parser.parse_args()
    if args.hits_only:
        hits_heatmap()
        path=DATA/'render_manifest.json'
        manifest=json.loads(path.read_text())
        manifest.update(hit_dendrograms='Identical average-linkage tree above and left of the heatmap; 56 observed genes only.',
                        hit_size_inches=[10,10.7],hit_row_labels='Right',
                        hit_clustering_manifest='hit_clustering_manifest.json',figure_roles=FIGURE_ROLES)
        path.write_text(json.dumps(manifest,indent=2),encoding='utf-8')
        print(json.dumps(manifest,indent=2))
        return
    hits_heatmap(); saga_heatmaps()
    pages={}
    pages['hits']=partner_atlas(A.index[A.screen_hit.eq(True)].tolist(),'ch04_supp_ht1080_depmap_hit_partner_atlas')
    pages['saga']=partner_atlas(A.index[A.saga_module.ne('')].tolist(),'ch04_supp_ht1080_depmap_saga_partner_atlas')
    manifest={'font':FONT,'heatmap_scale':[-1,1],'missing_value_colour':'#D8DCE0','figure_roles':FIGURE_ROLES,
              'hit_order':'Average linkage on 1 - Pearson r; missing ZCCHC16 appended and grey.',
              'hit_dendrograms':'Identical average-linkage tree above and left of the heatmap; 56 observed genes only.',
              'hit_size_inches':[10,10.7],'hit_row_labels':'Right',
              'hit_clustering_manifest':'hit_clustering_manifest.json',
              'saga_order':'Frozen structural module order, then four alternative members.',
              'partner_atlas_pages':pages,'atlas_display':'Top 8 positive partners per seed, without a minimum r cutoff; same 0-1 r axis across panels.',
              'heatmap_exports':'PDF, editable SVG, 600 dpi PNG; 160 dpi review PNG.',
              'atlas_exports':'Multipage vector PDF and per-page editable SVG / 160 dpi review PNG.'}
    (DATA/'render_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(manifest,indent=2))


if __name__=='__main__': main()
