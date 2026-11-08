"""Presentation-only revision of saved TF results; no analysis is refitted."""
from make_tf_modules_figures import *

plt.rcParams.update({'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6,
                     'ytick.labelsize': 6, 'legend.fontsize': 6})

def mark(ax, s):
    ax.text(-.13, 1.04, s, transform=ax.transAxes, fontsize=12, weight='semibold')

def line_compact(line):
    R = regulator_rows(line)
    mem = tab(line, 'TF_module_membership').sort_values('order')
    S = tab(line, 'TF_module_similarity_weighted').set_index('source')
    mods = tab(line, 'TF_modules').set_index('module')
    drv = tab(line, 'TF_module_driver_targets')
    n = len(R)
    H = max(145 if line=='DU145' else 100, n * 2.45 + 38)
    fig = plt.figure(figsize=(183*MM, H*MM))
    # Checklist spans the left column; overlap and contributing genes on the right.
    bottom, height = (.57,.26) if line=='A549' else (.17,.66)
    a = fig.add_axes([.10, bottom, .14, height])
    q = fig.add_axes([.25, bottom, .28, height])
    y = np.arange(n)[::-1]
    a.barh(y, R.ulm_score, color=np.where(R.ulm_score > 0, RED, BLUE), height=.7)
    a.axvline(0, color='#999999', lw=.5)
    a.set_yticks(y, R.source, fontsize=5.5 if n>30 else 6)
    for t, tier in zip(a.get_yticklabels(), R.tier):
        t.set_weight('semibold' if tier=='robust' else 'normal')
    a.set_xlabel('ULM score')
    a.set_xlim(min(-1,R.ulm_score.min()-1), max(1,R.ulm_score.max()+1))
    for ax in [a,q]: ax.set_ylim(-.7,n-.3)
    for i, rule in enumerate(RULES):
        ok=R[rule].astype(bool).to_numpy()
        q.scatter(np.full(n,i),y,s=10,facecolors=np.where(ok,DARK,'white'),edgecolors=DARK,lw=.5)
    for yy, p, m in zip(y,R.ulm_FDR,R.module):
        q.text(6,yy,pshort(p),fontsize=5.3,va='center')
        if pd.notna(m): q.text(8.1,yy,f'{int(m)}',fontsize=6,va='center',color=DARK)
    for x,label in enumerate(RULES): q.text(x,n-.1,label,ha='center',fontsize=6)
    q.text(6,n-.1,'ULM FDR',fontsize=6);q.text(8.1,n-.1,'Module',fontsize=6)
    q.set_xlim(-.5,9.6);q.axis('off');mark(a,'a')
    fig.text(.1,.07,'Bold: all six checks passed; filled/open circles: passed/failed',fontsize=6)
    b=fig.add_axes([.65,.59,.25,.33])
    order=mem.source.tolist()
    im=b.imshow(S.loc[order,order],cmap='magma_r',vmin=0,vmax=1)
    b.set_xticks(range(len(order)),order,rotation=90,fontsize=4.4 if n>30 else 6)
    b.set_yticks(range(len(order)),order,fontsize=4.4 if n>30 else 6)
    b.tick_params(length=1,pad=1)
    for m in mods.index:
        idx=np.flatnonzero(mem.module.to_numpy()==m)
        b.add_patch(Rectangle((idx.min()-.5,idx.min()-.5),len(idx),len(idx),fill=False,
                             ec=MOD_COLOUR[int(m)],lw=1.1,ls='-' if mods.loc[m,'stable'] else '--'))
    cb=fig.colorbar(im,cax=fig.add_axes([.93,.64,.012,.22]))
    cb.set_label('Weighted target overlap',fontsize=6);mark(b,'b')
    c=fig.add_axes([.70,.12 if line=='DU145' else .17,.23,.37 if line=='DU145' else .27])
    rows=[]
    for m in mods.index:
        ds=drv[drv.module==m].sort_values('summed_contribution',ascending=False).head(6 if len(mods)>1 else 8)
        for _,r in ds.iterrows(): rows.append((m,r.target,r.t))
    yy=np.arange(len(rows))[::-1]
    c.barh(yy,[r[2] for r in rows],color=[RED if r[2]>0 else BLUE for r in rows],height=.7)
    c.axvline(0,color='#999999',lw=.5)
    c.set_yticks(yy,[r[1] for r in rows],fontsize=5.5)
    for m in mods.index:
        ids=[i for i,r in enumerate(rows) if r[0]==m]
        cy=np.mean(yy[ids]); label=f'M{int(m)}'+(' *' if not mods.loc[m,'stable'] else '')
        c.text(-.55,cy,label,transform=c.get_yaxis_transform(),va='center',fontsize=6)
        if ids[-1]<len(rows)-1:c.axhline(yy[ids[-1]]-.5,color='#cccccc',lw=.5)
    c.set_xlabel('Target-gene moderated t');mark(c,'c')
    if not mods.stable.all():fig.text(.65,.035,'* Below module stability threshold',fontsize=6)
    export(fig,f'ch03_supp_{line.lower()}_tf_regulator_modules',
           {'regulators':R,'modules':mods.reset_index(),'drivers':drv,'similarity':S.reset_index()},
           {'inputs':set(USED),'figure_size_mm':[183,H],'display_note':'Dashed outlines mark modules below the stability threshold.'})

def crossline_compact():
    rs={l:tab(l,'TF_robustness').set_index('source') for l in LINES}
    keep=set.union(*[set(r.index[r.tier=='robust']) for r in rs.values()])|{'FOXA2','GRHL2'}
    def key(s):
        flags=[s in rs[l].index and rs[l].loc[s,'tier']=='robust' for l in LINES]
        return (-sum(flags),next((i for i,f in enumerate(flags) if f),9),s)
    order=['STAT1','MYC']+sorted(keep-{'STAT1','MYC'},key=key)
    score=pd.DataFrame({l:rs[l].ulm_score.reindex(order) for l in LINES})
    tiers=pd.DataFrame({l:rs[l].tier.reindex(order) for l in LINES})
    fig=plt.figure(figsize=(183*MM,158*MM))
    a=fig.add_axes([.12,.15,.18,.78])
    im=a.imshow(score,cmap=SIGNED,vmin=-8,vmax=8,aspect='auto')
    for i,s in enumerate(order):
        for j,l in enumerate(LINES):
            t=tiers.loc[s,l]
            if t in ['robust','supported']:
                a.scatter(j,i,s=10,facecolors=DARK if t=='robust' else 'none',edgecolors=DARK,lw=.6)
    a.set_yticks(range(len(order)),order,fontsize=5.5)
    a.set_xticks(range(3),LINES,fontsize=6);a.xaxis.tick_top();a.tick_params(length=0)
    a.axhline(1.5,color='white',lw=1.5)
    mark(a,'a')
    cb=fig.colorbar(im,cax=fig.add_axes([.12,.07,.18,.012]),orientation='horizontal',extend='both')
    cb.set_label('ULM score (clipped at ±8)',fontsize=6)
    a.legend(handles=[Line2D([],[],marker='o',ls='',color=DARK,ms=3,label='Robust'),
                      Line2D([],[],marker='o',ls='',mfc='none',mec=DARK,ms=3,label='Supported')],
             bbox_to_anchor=(0,-.015),loc='upper left',ncol=2,fontsize=5.5)
    b=fig.add_axes([.49,.72,.43,.19]); counts=[]
    names=['robust','supported','consistent_only','overall_only'];colors=['#333333','#888888','#bbbbbb','#e5e5e5']
    for i,l in enumerate(LINES):
        left=0
        for t,col in zip(names,colors):
            k=int((rs[l].tier==t).sum());counts.append({'line':l,'tier':t,'n':k})
            b.barh(2-i,k,left=left,color=col,height=.6,edgecolor='white',lw=.4);left+=k
        b.text(left+1,2-i,f'{left}/{len(rs[l])}',va='center',fontsize=6)
    b.set_yticks([2,1,0],LINES);b.set_xlim(0,70);b.set_xlabel('Regulators with ULM FDR <0.05')
    b.legend(handles=[Rectangle((0,0),1,1,color=c,label=t) for c,t in zip(colors,['Robust','Supported','R1 + R2 only','R1 only'])],
             ncol=2,loc='lower left',bbox_to_anchor=(0,1.03),fontsize=6);mark(b,'b')
    rows=[]
    for l in LINES:
        m=tab(l,'TF_modules').copy();u=tab(l,'TF_module_target_union_tests')
        m['line']=l;rows.append(m.merge(u[['module','fry_p']],on='module'))
    D=pd.concat(rows,ignore_index=True)
    c=fig.add_axes([.49,.16,.23,.39]);yy=np.arange(len(D))[::-1]
    c.axvspan(.75,1.02,color='#f0f0f0');c.axvline(.75,color='#aaaaaa',lw=.5,ls='--')
    for y,(_,r) in zip(yy,D.iterrows()):
        col=LINE_COLOUR[r.line]
        c.plot([r.bootstrap_jaccard,r.edge_drop_jaccard],[y,y],color=col,lw=.7)
        c.scatter(r.bootstrap_jaccard,y,s=20,color=col)
        c.scatter(r.edge_drop_jaccard,y,s=20,marker='s',facecolors='white',edgecolors=col)
        c.text(1.07,y,pshort(r.fry_p),va='center',fontsize=6,transform=c.get_yaxis_transform())
        members=r.members.split(';');label=', '.join(members[:3])+(f' +{len(members)-3}' if len(members)>3 else '')
        c.text(1.40,y,label+'\n'+r.direction+(' *' if not r.stable else ''),va='center',fontsize=5.7,transform=c.get_yaxis_transform())
    c.set_yticks(yy,[f'{r.line} M{r.module}' for _,r in D.iterrows()],fontsize=6)
    c.set_xlim(.55,1.02);c.set_ylim(-.6,len(D)-.3);c.set_xlabel('Mean module Jaccard')
    c.text(1.07,len(D)-.3,'fry P',fontsize=6,transform=c.get_yaxis_transform())
    c.text(1.40,len(D)-.3,'Members; direction',fontsize=6,transform=c.get_yaxis_transform())
    c.legend(handles=[Line2D([],[],marker='o',ls='',color='#777777',label='Gene bootstrap'),
                      Line2D([],[],marker='s',ls='',mfc='white',mec='#777777',label='Edge drop')],
             loc='lower left',bbox_to_anchor=(0,1.1),fontsize=6)
    fig.text(.49,.075,'Shading: stability ≥0.75 in both tests\n* Below stability threshold; fry P values are unadjusted',fontsize=6)
    mark(c,'c')
    export(fig,'ch03_main_crossline_tf_regulators',{'scores':score.reset_index(),'tiers':tiers.reset_index(),
            'counts':pd.DataFrame(counts),'modules':D},{'inputs':set(USED),'display_note':'Module fry P values are unadjusted saved results.'})

def sensitivity_figure():
    fig=plt.figure(figsize=(183*MM,173*MM))
    a=fig.add_axes([.29,.75,.64,.20]);tests=list(NULL_LABEL);nr=[]
    for j,l in enumerate(LINES):
        df=tab(l,'TF_module_nulls').set_index('test')
        for i,t in enumerate(tests):
            r=df.loc[t];a.scatter(r.empirical_p,3-i+(1-j)*.22,s=17,color=LINE_COLOUR[l])
            nr.append({'line':l,'test':t,**r.to_dict()})
    a.set_xscale('log');a.axvline(.05,color='#999999',lw=.6,ls='--');a.set_xlim(.0005,2)
    a.set_yticks(range(4),[NULL_LABEL[t] for t in tests[::-1]],fontsize=6)
    a.set_xlabel('Empirical P versus null');mark(a,'a')
    a.legend(handles=[Line2D([],[],marker='o',ls='',color=LINE_COLOUR[l],label=l) for l in LINES],
             ncol=3,loc='lower left',bbox_to_anchor=(0,1.02),fontsize=6)
    allrows=[];uts=[]
    for l in LINES:
        co=tab(l,'TF_module_confounders').copy();co['line']=l;allrows.append(co)
        u=tab(l,'TF_module_target_union_tests').copy();u['line']=l;uts.append(u)
    U=pd.concat(uts);C=pd.concat(allrows)
    for j,(_,r) in enumerate(U.iterrows()):
        ax=fig.add_axes([.28+j*.138,.22,.12,.34]);cc=C[(C.line==r.line)&(C.module==r.module)].set_index('programme').reindex(PROGRAMME)
        yy=np.arange(len(cc))[::-1];fail=cc.collapses.astype(bool)
        ax.axvline(r.union_ulm_score,color=LINE_COLOUR[r.line],lw=.8)
        ax.scatter(cc.removed_score,yy,s=16,facecolors=np.where(fail,'white',DARK),edgecolors=np.where(fail,RED,DARK),lw=.7)
        ax.set_yticks(yy,list(PROGRAMME.values()) if j==0 else [],fontsize=6)
        ax.set_ylim(-.6,len(cc)-.4);ax.tick_params(axis='x',labelsize=5.5)
        ax.set_xlabel(f'{r.line}\nM{r.module}',fontsize=6)
        ax.text(.5,1.04,f'perm. P\n{pshort(r.perm_p)}\nfry P\n{pshort(r.fry_p)}',transform=ax.transAxes,ha='center',fontsize=6)
        if j==0:mark(ax,'b')
    fig.text(.43,.135,'Signed target-union score after gene-set exclusion',fontsize=7)
    fig.legend(handles=[Line2D([],[],marker='o',ls='',color=DARK,label='Same direction, P <0.05'),
                        Line2D([],[],marker='o',ls='',mfc='white',mec=RED,label='No longer significant'),
                        Line2D([],[],color='#777777',lw=.8,label='All targets')],
               loc='lower center',bbox_to_anchor=(.57,.075),ncol=3,fontsize=6)
    export(fig,'ch03_supp_crossline_tf_module_sensitivity',{'nulls':pd.DataFrame(nr),'exclusions':C,'union_tests':U},
           {'inputs':set(USED),'display_note':'Positive signed scores denote agreement with inferred regulator directions, not necessarily increased expression.'})

if __name__=='__main__':
    for l in LINES:
        USED.clear();line_compact(l)
    USED.clear();crossline_compact()
    USED.clear();sensitivity_figure()
