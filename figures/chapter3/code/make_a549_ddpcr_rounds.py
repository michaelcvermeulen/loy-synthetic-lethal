"""A549 ddPCR guide comparison using instrument concentration ratios and explicit round labels."""
from pathlib import Path
import argparse, sys, json, hashlib, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
THESIS=Path(__file__).resolve().parents[3]
CHAPTER=THESIS/'figures/chapter3'
STYLE=THESIS/'skills/thesis-figure-style'
sys.path.insert(0,str(STYLE/'scripts'))
from thesis_style import apply_style
MAIN='ch03_main_a549_ddpcr_nucleofection_rounds'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=CHAPTER)
    parser.add_argument('--raw-root',type=Path,default=Path('A:/Dropbox/LOY_paralogs/data/ddPCR'))
    parser.add_argument('--dpi',type=int,default=600)
    args=parser.parse_args();out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    family=apply_style(8.5)
    font=Path(font_manager.findfont(font_manager.FontProperties(family=family)))
    assert font.resolve()==(STYLE/'assets/fonts/OpenSans-Regular.ttf').resolve()
    target=CHAPTER/'source_data/a549_ddpcr_nucleofection_rounds';target.mkdir(exist_ok=True)
    snapshots=target/'raw_exports';snapshots.mkdir(exist_ok=True)
    master_path=args.raw_root/'summary_table_ddPCR.csv'
    master=pd.read_csv(master_path)
    metadata=master[master.cell_line.eq('A549')].copy()
    metadata.to_csv(target/'a549_master_metadata.csv',index=False)
    names={'CP':'CP','Q53':'Q53','CP+Q53':'CP + Q53','CTRL2':'CTRL2','CTRL':'CTRL2','MIX/MIX':'MIX/MIX'}
    mapping={}
    for sid,g in metadata.groupby('Sample'):
        labels={names[v] for v in g['name'] if v in names}
        if labels:
            assert len(labels)==1,(sid,labels)
            label=next(iter(labels));mapping[sid]=(label,2 if label=='MIX/MIX' else 1)
    # January exports use condition labels; the master links these to DNA IDs E1/E2/E3/E5.
    jan_alias={'A549_CP':'E1','A549_Q53':'E3','A549_CP+Q53':'E5','A549_Q53+CP':'E5','A549_CTRL2':'E2'}
    runs=[('2025-01-26','AMELX Craig_JAN_26_2025_2025-01-26-21-39'),
          ('2025-01-27','AMELX Craig_JAN_27_2025_2025-01-27-16-45'),
          ('2025-02-24','AMELX Craig_FEB_24_2025_2025-02-25-10-57'),
          ('2025-02-27','AMELX Craig_FEB_27_2025_2025-02-27-17-17'),
          ('2025-03-09','AMELX Craig_MARCH_9_2025_2025-03-09-22-48'),
          ('2025-03-25','AMELX Craig_MARCH_25_2025_2025-03-25-14-12')]
    records=[];manifest=[]
    for date,folder in runs:
        path=args.raw_root/folder/'data_table.csv';d=pd.read_csv(path)
        conc=next(c for c in d if c.startswith('Conc(copies/'))
        for _,y in d[d.Target.isin(['AMELY','DDX3Y'])].iterrows():
            source_label=y['Sample description 1']
            sid=jan_alias.get(source_label,source_label) if date.startswith('2025-01') else source_label
            if sid not in mapping:continue
            condition,round_n=mapping[sid]
            ref='AMELX' if y.Target=='AMELY' else 'RNaseP'
            refs=d[d.Well.eq(y.Well)&d.Target.str.upper().eq(ref.upper())];assert len(refs)==1
            x=refs.iloc[0];yc=float(y[conc]);xc=float(x[conc]);assert xc>0
            ratio=yc/xc;assert np.isfinite(ratio)
            assert np.isclose(ratio,float(y['Ratio']))
            assert set(master.loc[master.Sample.eq(sid),'cell_line'])=={'A549'}
            excluded=sid=='I1' and y.Target=='AMELY'
            records.append(dict(date=date,condition=condition,sample_id=sid,source_label=source_label,
                export_well=y.Well,merged_wells=y.get('MergedWells',''),target=y.Target,reference=ref,
                assay=y.Target+'/'+ref,target_concentration=yc,reference_concentration=xc,
                concentration_ratio=ratio,instrument_ratio=y['Ratio'],accepted_droplets=y['Accepted Droplets'],
                target_status=y['Status'],target_status_reason=y['Status Reason'],reference_status=x['Status'],
                source_path=str(path),nucleofection_count=round_n,included=not excluded,
                exclusion_reason='Historical A549 script excludes I1 AMELY; retained without re-diagnosing assay discordance' if excluded else ''))
        snap=date+'_data_table.csv';shutil.copy2(path,snapshots/snap)
        manifest.append(dict(source_path=str(path),snapshot='raw_exports/'+snap,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    manifest.append(dict(source_path=str(master_path),sha256=hashlib.sha256(master_path.read_bytes()).hexdigest()))
    all_records=pd.DataFrame(records)
    all_records.to_csv(target/'all_assessed_measurements.csv',index=False)
    all_records[~all_records.included].to_csv(target/'excluded_measurements.csv',index=False)
    p=all_records[all_records.included].copy()
    assert set(p[p.nucleofection_count.eq(2)].sample_id)=={'L3','L6'}
    assert len(all_records[~all_records.included])==2
    assert not p.duplicated(['date','export_well','target']).any()
    reference_copies={'AMELY/AMELX': 2, 'DDX3Y/RNaseP': 3}
    p['reference_copy_number']=p.assay.map(reference_copies)
    p['estimated_target_copies_per_cell']=p.concentration_ratio*p.reference_copy_number
    p.to_csv(target/'main_measurements.csv',index=False)
    pd.DataFrame([dict(sample_id=k,condition=v[0],nucleofection_count=v[1],basis='Saved A549 plotting script: MIX/MIX is CP+Q53 twice; Version is not round count') for k,v in mapping.items()]).to_csv(target/'sample_round_mapping.csv',index=False)
    (target/'raw_source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    for scaled in [False,True]:
        stem=MAIN if scaled else 'ch03_supp_a549_ddpcr_raw_ratios'
        value_column='estimated_target_copies_per_cell' if scaled else 'concentration_ratio'
        groups=['CP','Q53','CP + Q53','MIX/MIX','CTRL2'];assays=['AMELY/AMELX','DDX3Y/RNaseP'];summaries=[]
        fig,axes=plt.subplots(1,2,figsize=(183/25.4,108/25.4))
        fig.subplots_adjust(left=.093,right=.987,bottom=.23,top=.85,wspace=.25)
        for ax,assay,letter in zip(axes,assays,['a','b']):
            ymax=1.35 if scaled else (.70 if assay==assays[0] else .44)
            for i,guide in enumerate(groups):
                g=p[p.assay.eq(assay)&p.condition.eq(guide)];vals=g[value_column].to_numpy()
                if len(vals):
                    med=float(np.median(vals));q1,q3=np.quantile(vals,[.25,.75]);color='#3C5488' if guide=='CTRL2' else '#E64B35'
                    ax.bar(i,med,width=.67,color=color,edgecolor=color,linewidth=.7,hatch='///' if guide=='MIX/MIX' else None,zorder=2)
                    if guide=='MIX/MIX':
                        ax.bar(i,med,width=.67,facecolor='none',edgecolor='#87301E',linewidth=.7,hatch='///',zorder=3)
                    ax.errorbar(i,med,yerr=[[med-q1],[q3-med]],fmt='none',ecolor='#333333',elinewidth=.8,capsize=2,capthick=.8,zorder=5)
                    offsets=np.linspace(-.19,.19,len(vals)) if len(vals)>1 else np.array([0.])
                    ax.scatter(i+offsets,vals,s=14,facecolor='white',edgecolor='#222222',linewidth=.65,zorder=6)
                    y=max(vals.max(),q3)+ymax*.044
                    if scaled and .93 < y < 1.05:y=1.05
                    ax.text(i,y,f'{med:.3f}',ha='center',va='bottom',fontsize=7.4,fontweight='semibold')
                    ax.text(i,y+ymax*.047,f'n={len(vals)}',ha='center',va='bottom',fontsize=6.8,color='#666666')
                else:
                    med=q1=q3=None
                    ax.text(i,.045,'NA',transform=ax.get_xaxis_transform(),ha='center',fontsize=8,color='#777777')
                summaries.append(dict(assay=assay,reference_copy_number=reference_copies[assay],unit='estimated target copies per cell' if scaled else 'target/reference concentration ratio',condition=guide,nucleofection_count=2 if guide=='MIX/MIX' else 1,n_measurements=len(vals),n_sample_ids=g.sample_id.nunique(),median=med,q1=q1,q3=q3,status='recorded' if len(vals) else 'not available'))
            if scaled:ax.axhline(1,color='#888888',lw=.8,ls=(0,(4,3)),zorder=1)
            ax.axvline(3.5,color='#DDDDDD',lw=.65)
            ax.set_xlim(-.55,4.55);ax.set_ylim(0,ymax)
            ax.set_xticks(range(5),['CP\n1\u00d7','Q53\n1\u00d7','CP + Q53\n1\u00d7','CP + Q53\n2\u00d7','sgCTRL'])
            ax.tick_params(axis='x',length=0,pad=8,labelsize=7.5)
            ax.set_yticks(np.arange(0,1.21,.2) if scaled else (np.arange(0,.61,.1) if assay==assays[0] else np.arange(0,.41,.1)))
            ax.set_title(assay.split('/')[0] if scaled else assay.replace('/',' / '),pad=12,fontsize=10)
            ax.text(-.14,1.08,letter,transform=ax.transAxes,fontsize=12,fontweight='semibold')
            ax.grid(axis='y',color='#ECECEC',lw=.45,zorder=0);ax.set_axisbelow(True)
        axes[0].set_ylabel('Estimated Y-linked target copies per cell' if scaled else 'Target/reference concentration ratio',labelpad=9)
        fig.text(.54,.065,'1\u00d7 = one nucleofection; 2\u00d7 = two successive nucleofections',ha='center',fontsize=8,color='#444444')
        for ext in ['pdf','svg','png']:fig.savefig(out/f'{stem}.{ext}',dpi=args.dpi,facecolor='white')
        fig.savefig(out/f'{stem}_preview.png',dpi=170,facecolor='white');plt.close(fig)
        pd.DataFrame(summaries).to_csv(target/('main_group_summary.csv' if scaled else 'raw_ratio_group_summary.csv'),index=False)
    meta=dict(font=str(font),dimensions_mm=[183,108],dpi=args.dpi,main_measurements=len(p),
        excluded_measurements=2,double_nucleofection_sample_ids=['L3','L6'],version_used_as_round=False,
        summary='median and interquartile range; all individual assay measurements; no inferential tests',
        source_concentration_recheck='passed',normalization='target/reference concentration ratio; no copy-number scaling or control normalization')
    meta.update(reference_copy_numbers=reference_copies,normalization='concentration ratio multiplied by reference-locus copy number; no control normalization',main_axis_limits=[0,1.35],reference_line=1)
    (target/'figure_validation.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print(pd.DataFrame(summaries).to_string(index=False))
if __name__=='__main__':main()
