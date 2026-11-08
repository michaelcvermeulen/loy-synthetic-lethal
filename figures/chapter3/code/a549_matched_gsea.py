"""HT1080-matched GSEA appearance, preserving A549-specific analyses."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm

LABELS={
'HALLMARK_PI3K_AKT_MTOR_SIGNALING':'[H] PI3K–AKT–mTOR signalling',
'HALLMARK_DNA_REPAIR':'[H] DNA repair',
'HALLMARK_UV_RESPONSE_UP':'[H] UV response up',
'HALLMARK_UV_RESPONSE_DN':'[H] UV response down',
'HALLMARK_IL6_JAK_STAT3_SIGNALING':'[H] IL-6/JAK/STAT3 signalling',
'HALLMARK_IL2_STAT5_SIGNALING':'[H] IL-2/STAT5 signalling',
'HALLMARK_TGF_BETA_SIGNALING':'[H] TGF-beta signalling',
'HALLMARK_WNT_BETA_CATENIN_SIGNALING':'[H] WNT/beta-catenin signalling',
'HALLMARK_KRAS_SIGNALING_UP':'[H] KRAS signalling up',
'HALLMARK_KRAS_SIGNALING_DN':'[H] KRAS signalling down',
'HALLMARK_MYC_TARGETS_V1':'[H] MYC targets V1',
'HALLMARK_MYC_TARGETS_V2':'[H] MYC targets V2',
'HALLMARK_E2F_TARGETS':'[H] E2F targets',
'HALLMARK_G2M_CHECKPOINT':'[H] G2M checkpoint',
'HALLMARK_MTORC1_SIGNALING':'[H] mTORC1 signalling',
'HALLMARK_P53_PATHWAY':'[H] p53 pathway',
'HALLMARK_TNFA_SIGNALING_VIA_NFKB':'[H] TNF-alpha signalling via NF-κB',
'HALLMARK_OXIDATIVE_PHOSPHORYLATION':'[H] Oxidative phosphorylation',
'REACTOME_CELL_CYCLE_CHECKPOINTS':'[R] Cell-cycle checkpoints',
'REACTOME_DNA_REPLICATION':'[R] DNA replication',
'REACTOME_DNA_REPAIR':'[R] DNA repair',
'GOBP_RIBOSOME_BIOGENESIS':'[GO] Ribosome biogenesis',
'GOBP_MITOCHONDRIAL_GENE_EXPRESSION':'[GO] Mitochondrial gene expression',
'GOBP_TRNA_METABOLIC_PROCESS':'[GO] tRNA metabolism',
'REACTOME_FORMATION_OF_FIBRIN_CLOT_CLOTTING_CASCADE':'[R] Fibrin-clot formation',
'REACTOME_GPCR_LIGAND_BINDING':'[R] GPCR ligand binding',
'GOBP_PROTEIN_ACTIVATION_CASCADE':'[GO] Protein activation cascade'}
THEMES=[
('Cell cycle and DNA maintenance',['HALLMARK_MYC_TARGETS_V1','HALLMARK_MYC_TARGETS_V2','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','REACTOME_CELL_CYCLE_CHECKPOINTS','REACTOME_DNA_REPLICATION','REACTOME_DNA_REPAIR','HALLMARK_P53_PATHWAY']),
('RNA processing and metabolism',['GOBP_RIBOSOME_BIOGENESIS','GOBP_TRNA_METABOLIC_PROCESS','GOBP_MITOCHONDRIAL_GENE_EXPRESSION','HALLMARK_MTORC1_SIGNALING','HALLMARK_OXIDATIVE_PHOSPHORYLATION']),
('Inflammatory signalling',['HALLMARK_TNFA_SIGNALING_VIA_NFKB']),
('Coagulation and ligand-associated sets',['REACTOME_FORMATION_OF_FIBRIN_CLOT_CLOTTING_CASCADE','GOBP_PROTEIN_ACTIVATION_CASCADE','REACTOME_GPCR_LIGAND_BINDING'])]

def draw(data,fixed,broad=False):
    assert data.NES.abs().max()<=2.6, 'Do not silently clip enrichment scores'
    height=6.4 if broad else 10.0
    fig=plt.figure(figsize=(183/25.4,height))
    themes=THEMES if broad else [('Negatively enriched sets',list(data.index[data.NES<0])),('Positively enriched sets',list(data.index[data.NES>=0]))]
    rows=[];heads=[];values=[]
    for title,keys in themes:
        keys=[x for x in keys if x in data.index]
        if not keys:continue
        heads.append((len(values),title,float(data.loc[keys,'NES'].mean())>0));values.append([np.nan])
        for key in keys:
            rows.append((len(values),key));values.append([float(data.loc[key,'NES'])])
    assert len(rows)==len(data)
    ax=fig.add_axes([.51,.16,.12,.66 if broad else .70])
    cmap=LinearSegmentedColormap.from_list('NES',['#3C5488','#FAFAFA','#E64B35'])
    im=ax.imshow(values,aspect='auto',cmap=cmap,norm=TwoSlopeNorm(vmin=-2.6,vcenter=0,vmax=2.6),interpolation='nearest')
    labels=[LABELS.get(key,'[H] '+key.removeprefix('HALLMARK_').replace('_',' ').capitalize().replace('signaling','signalling')) for _,key in rows]
    ax.set_yticks([i for i,_ in rows],labels,fontsize=6.4 if broad else 6.1)
    ax.set_xticks([0],['Overall'],fontsize=7);ax.tick_params(length=0);ax.spines[:].set_visible(False)
    for i,title,positive in heads:
        ax.text(-.56,i,title,ha='right',va='center',fontsize=6.8 if broad else 6.5,weight='semibold',color='#B5422B' if positive else '#3C5488',clip_on=False)
    stats=fig.add_axes([.66,.16,.31,.66 if broad else .70])
    stats.set(xlim=(-.5,1.5),ylim=(len(values)-.5,-.5));stats.axis('off')
    for i,key in rows:
        nes=float(data.loc[key,'NES']);fdr=float(data.loc[key,'BH_across_collections'])
        if fdr<.05:ax.text(0,i,'*',ha='center',va='center',fontsize=9,color='white' if abs(nes)>2.6*.55 else '#242424')
        vals=[fdr,float(fixed.loc[key,'BH_across_collections'])]
        for j,val in enumerate(vals):
            stats.text(j,i,f'{val:.1e}',ha='center',va='center',fontsize=6.2 if broad else 5.9,weight='semibold' if val<.05 else 'normal',color='#444444')
    for j,label in enumerate(['fgsea','CAMERA\ndefault (0.01)']):stats.text(j,-.85,label,ha='center',va='bottom',fontsize=6.5)
    fig.text(.035,.945,'a',fontsize=12,weight='semibold')
    fig.text(.08,.945,'Gene-set enrichment' if broad else 'Complete Hallmark enrichment',fontsize=10,weight='semibold')
    fig.text(.66,.945,'b',fontsize=12,weight='semibold')
    fig.text(.70,.945,'Global FDR',fontsize=9,weight='semibold')
    fig.text(.08,.907,'All collections | 23,457 sets tested',fontsize=7)
    if broad:fig.text(.08,.878,'1,366 significant: 1,300 lower, 66 higher | selected sets shown',fontsize=6.8,color='#555555')
    cax=fig.add_axes([.34,.080,.35,.017]);cb=fig.colorbar(im,cax=cax,orientation='horizontal')
    cb.ax.tick_params(labelsize=7);cb.ax.set_title('Normalised enrichment score (NES)',fontsize=7,pad=6)
    fig.text(.08,.038,'* fgsea FDR < 0.05. Bold values: FDR < 0.05 in the indicated analysis.',fontsize=6.5)
    fig.text(.08,.018,'H, Hallmark; R, Reactome; GO, GO Biological Process.',fontsize=6.5)
    display=data.loc[[key for _,key in rows]].copy();display['display_label']=labels
    return fig,display
