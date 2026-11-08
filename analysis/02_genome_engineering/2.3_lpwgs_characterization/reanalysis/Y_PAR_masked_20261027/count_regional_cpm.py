"""Count primary read ends in complete regions and exact plotted-bin unions."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import csv,hashlib,json,subprocess

ROOT=Path(__file__).resolve().parent
THESIS=ROOT.parents[2]/'SEPT_11_2026_THESIS'
OUT=ROOT/'regional_cpm_20261028'
OUT.mkdir(exist_ok=True)
(OUT/'beds').mkdir(exist_ok=True)
SAMTOOLS=str(ROOT/'runtime/env/bin/samtools')
PAIRS={'A549':['A549-CTRL2','A549-LOY'],'HT1080':['HT1080-sgCTRL-D0','HT1080-LOY-D0'],
       'DU145':['LNCaP-sgCTRL','LNCaP-LOY-2X']}
REGIONS={'MSY_Yp':('chrY',2458320,10565750),'PAR1_X':('chrX',0,2394410),
         'euchromatic_Yq':('chrY',10883085,27449937),'PAR2_X':('chrX',153925834,154259566),
         'Yq_including_Yq12_excluding_PAR2':('chrY',10883085,62122809),
         'Yq12':('chrY',27449937,62122809),'distal_MSY_Yp':('chrY',2458320,4800000),
         'proximal_MSY_Yp':('chrY',4800000,10565750)}
PRIMARY=['MSY_Yp','PAR1_X','euchromatic_Yq','PAR2_X']
sources={};bedsets={};figure_medians={}

def rows(path):
    with path.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))

def save(path,data):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)

def write_bed(path,intervals):
    merged=[]
    for chrom,start,end in sorted(intervals):
        if merged and merged[-1][0]==chrom and start<=merged[-1][2]:
            merged[-1]=(chrom,merged[-1][1],max(end,merged[-1][2]))
        else:merged.append((chrom,start,end))
    assert merged
    path.write_text(''.join(f'{c}\t{s}\t{e}\n' for c,s,e in merged))
    return {'path':str(path),'length_bp':sum(e-s for c,s,e in merged),
            'n_intervals':len(merged),'n_input_bins':len(intervals)}

for region,interval in REGIONS.items():
    bedsets[('whole_region','all',region)]=write_bed(OUT/'beds'/f'whole_{region}.bed',[interval])

for cell in PAIRS:
    slug=cell.lower()
    base=THESIS/'figures/chapter3/drafts'
    base=base/'a549_dense_bins_20261027/100kb' if cell=='A549' else base/'ht1080_du145_ypar_20261027'/slug
    source=base/'source_data/regional_y_atlases_normalised'/slug
    sources[cell]=source
    y=[r for r in rows(source/'paired_Y_bins_and_exclusions.csv') if r['included'].lower()=='true']
    par=[r for r in rows(source/'X_PAR_bin_ratios.csv') if r['included'].lower()=='true']
    labels={'MSY_Yp':['MSY Yp','distal MSY Yp','proximal MSY Yp'],'euchromatic_Yq':['euchromatic Yq']}
    if cell=='DU145':labels.update({'distal_MSY_Yp':['distal MSY Yp'],'proximal_MSY_Yp':['proximal MSY Yp']})
    for region,names in labels.items():
        intervals=[('chrY',int(float(r['start']))-1,int(float(r['end']))) for r in y if r['region'] in names]
        bedsets[('plotted_bins',cell,region)]=write_bed(OUT/'beds'/f'{cell}_{region}_plotted.bed',intervals)
    for region in ['PAR1_X','PAR2_X']:
        intervals=[('chrX',int(r['start0']),int(r['end'])) for r in par if r['track']==region]
        bedsets[('plotted_bins',cell,region)]=write_bed(OUT/'beds'/f'{cell}_{region}_plotted.bed',intervals)
    names={'MSY Yp':'MSY_Yp','euchromatic Yq':'euchromatic_Yq','distal MSY Yp':'distal_MSY_Yp','proximal MSY Yp':'proximal_MSY_Yp'}
    for r in rows(source/'regional_median_paired_ratios.csv'):
        figure_medians[(cell,names[r['region']])]=float(r['median_percent_of_sgCTRL'])
    for r in rows(source/'X_PAR_regional_medians.csv'):
        figure_medians[(cell,r.get('track',r.get('region')))]=float(r['median_percent_of_sgCTRL'])

def count(bam,bed=None,region=None):
    cmd=[SAMTOOLS,'view','-c','-q','37','-F','0xF04']
    if bed:cmd+=['-M','-L',bed]
    elif not region:cmd+=['-@','2']
    cmd.append(str(bam))
    if region:cmd.append(region)
    return int(subprocess.run(cmd,check=True,capture_output=True,text=True).stdout.strip())

def process(item):
    cell,condition,sample=item
    bam=ROOT/'alignment'/f'{sample}.YPARmasked.bam'
    marker=ROOT/'qc'/('REALIGNMENT_COMPLETE.txt' if sample=='A549-CTRL2' else f'{sample}.REALIGNMENT_COMPLETE.txt')
    assert marker.exists() and bam.exists() and Path(str(bam)+'.bai').exists()
    stat=(bam.stat().st_size,bam.stat().st_mtime_ns)
    subprocess.run([SAMTOOLS,'quickcheck',str(bam)],check=True)
    print('Counting library denominator:',sample,flush=True)
    denominator=count(bam)
    assert denominator>5_000_000
    result=[]
    for (scope,line,region),bed in bedsets.items():
        if line not in ['all',cell]:continue
        if cell!='DU145' and region in ['distal_MSY_Yp','proximal_MSY_Yp']:continue
        n=count(bam,bed=bed['path'])
        result.append({'cell':cell,'condition':condition,'sample':sample,'scope':scope,'region':region,
            'primary_read_ends':n,'eligible_library_read_ends':denominator,'CPM':n/denominator*1e6,
            'region_length_bp':bed['length_bp'],'n_intervals':bed['n_intervals'],
            'n_plotted_bins':bed['n_input_bins'] if scope=='plotted_bins' else ''})
    for region in ['chrY:1-2458320','chrY:62122810-62460029']:
        assert count(bam,region=region)==0,'Masked Y-PAR unexpectedly has high-quality mapped reads'
    prior=(ROOT/'a549_pair_comparison/regional_coverage.csv' if cell=='A549' else ROOT/'atlas_pair_coverage'/cell/'masked/regional_PAR_coverage.csv')
    for r in rows(prior):
        if r['sample']!=sample or r['chromosome']!='chrX' or r['alignment'] not in ['masked','Y_PAR_masked']:continue
        observed=next(z['primary_read_ends'] for z in result if z['scope']=='whole_region' and z['region']==r['region']+'_X')
        assert observed==int(r['read_count']),(sample,r,observed)
    for z in result:
        if z['scope']=='plotted_bins':
            full=next(r for r in result if r['scope']=='whole_region' and r['region']==z['region'])
            assert z['primary_read_ends']<=full['primary_read_ends']
    assert stat==(bam.stat().st_size,bam.stat().st_mtime_ns),'Input BAM changed'
    (OUT/f'{sample}.counts.json').write_text(json.dumps({'bam':str(bam),'size':stat[0],'mtime_ns':stat[1],'rows':result},indent=2))
    print('Complete:',sample,denominator,'eligible read ends',flush=True)
    return result

work=[(cell,condition,sample) for cell,pair in PAIRS.items() for condition,sample in zip(['sgCTRL','Y-depleted'],pair)]
with ThreadPoolExecutor(max_workers=3) as pool:
    allrows=[r for group in pool.map(process,work) for r in group]
save(OUT/'regional_reads_and_cpm_all_definitions.csv',allrows)
save(OUT/'whole_region_reads_and_cpm.csv',[r for r in allrows if r['scope']=='whole_region' and r['region'] in PRIMARY])
save(OUT/'plotted_bins_reads_and_cpm.csv',[r for r in allrows if r['scope']=='plotted_bins'])
comparisons=[]
for cell in PAIRS:
    for scope in ['whole_region','plotted_bins']:
        subset=[r for r in allrows if r['cell']==cell and r['scope']==scope]
        for ctrl in [r for r in subset if r['condition']=='sgCTRL']:
            loy=next(r for r in subset if r['condition']=='Y-depleted' and r['region']==ctrl['region'])
            comparisons.append({'cell':cell,'scope':scope,'region':ctrl['region'],
                'control_read_ends':ctrl['primary_read_ends'],'LOY_read_ends':loy['primary_read_ends'],
                'control_CPM':ctrl['CPM'],'LOY_CPM':loy['CPM'],
                'LOY_CPM_percent_of_control':100*loy['CPM']/ctrl['CPM'],
                'figure_median_percent_of_control':figure_medians.get((cell,ctrl['region']),''),
                'n_plotted_bins':ctrl['n_plotted_bins']})
save(OUT/'CPM_ratios_vs_figure_medians.csv',comparisons)
source_files=[p/f for p in sources.values() for f in ['paired_Y_bins_and_exclusions.csv','X_PAR_bin_ratios.csv','regional_median_paired_ratios.csv','X_PAR_regional_medians.csv']]
(OUT/'provenance.json').write_text(json.dumps({'minimum_MAPQ':37,'excluded_flags':'0xF04',
    'excluded_categories':['unmapped','secondary','supplementary','QC-failed','marked duplicates'],
    'count_unit':'primary read end; R1 and R2 counted separately, not read pairs',
    'CPM_denominator':'all mapped primary read ends passing identical MAPQ/flag filters across all reference contigs, including chrM',
    'interval_rule':'samtools multi-region iterator counts each alignment once per region union; do not sum per-bin overlap counts',
    'PAR_interpretation':'X-mapped combined X/Y signal after Y-PAR masking',
    'main_Yq_definition':'euchromatic chrY:10883086-27449937, matching atlas anatomical definition',
    'full_Yq_including_Yq12_also_reported':True,'regions_0based_halfopen':REGIONS,
    'samtools_version':subprocess.run([SAMTOOLS,'--version'],check=True,capture_output=True,text=True).stdout.splitlines()[0],
    'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files},
    'prior_PAR_counts_match':True,'input_BAMs_unchanged_by_size_mtime':True},indent=2))
print('Saved',OUT,flush=True)
