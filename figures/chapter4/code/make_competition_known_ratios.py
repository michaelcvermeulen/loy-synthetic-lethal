"""Known-ratio recovery from saved September measurements; no FCS reanalysis."""
from pathlib import Path
import hashlib
import json
import os
import sys
import xml.etree.ElementTree as ET

THESIS = Path(__file__).resolve().parents[3]
PROJECTS = THESIS.parent
sys.path.insert(0, str(PROJECTS / '.figure_python_deps'))
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/competition_known_ratios'
REVIEW = OUT / 'review'
DATA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pymupdf
from PIL import Image
from thesis_style import apply_style

STEM = 'ch04_main_competition_known_ratios_draft'
SOURCE = PROJECTS / 'comp_assay/analysis/spectral_competition_20260901/results.json'
CHECKS = THESIS / 'research/chapter4_competition_qc_20261019/numerical_checks.json'
ORIENTATIONS = {
    'LOY_GFP': ('LOY–GFP / CTRL–mCherry', '#32965B', 'o', '-'),
    'CTRL_GFP': ('LOY–mCherry / CTRL–GFP', '#D34F4F', 's', '--'),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare():
    source = json.loads(SOURCE.read_text())
    previous = json.loads(CHECKS.read_text())
    expected = next(e['sha256'] for e in previous['sources'] if Path(e['path']).name == 'results.json')
    assert sha(SOURCE) == expected, 'Saved analysis changed since the existing thesis check.'
    wells = pd.DataFrame(source['ratio_ladder'])
    assert len(wells) == 20 and wells['sample'].is_unique
    assert set(wells.orientation) == set(ORIENTATIONS)
    assert np.allclose(wells.observed_loy, np.where(wells.orientation.eq('LOY_GFP'),
                                                  wells.observed_gfp, 1-wells.observed_gfp))
    wells['nominal_percent'] = wells.nominal_loy * 100
    wells['observed_percent'] = wells.observed_loy * 100
    wells['deviation_pp'] = wells.observed_percent - wells.nominal_percent
    wells['orientation_label'] = wells.orientation.map({k:v[0] for k,v in ORIENTATIONS.items()})
    summary = wells.groupby(['orientation','nominal_percent']).observed_percent.agg(
        mean_percent='mean', minimum_percent='min', maximum_percent='max', technical_n='size').reset_index()
    assert len(summary) == 10 and summary.technical_n.eq(2).all()
    summary['deviation_pp'] = summary.mean_percent-summary.nominal_percent
    summary['technical_difference_pp'] = summary.maximum_percent-summary.minimum_percent
    for old in previous['september_ladder']:
        row = summary.loc[summary.orientation.eq(old['orientation']) &
                          summary.nominal_percent.eq(old['nominal_percent'])].iloc[0]
        assert np.isclose(row.mean_percent, old['mean_percent'], rtol=0, atol=1e-12)
        assert np.isclose(row.technical_difference_pp, old['technical_difference_pp'], rtol=0, atol=1e-12)
    wells.to_csv(DATA/'measurements.csv', index=False)
    summary.to_csv(DATA/'summary.csv', index=False)
    sources = [SOURCE, CHECKS,
               PROJECTS/'comp_assay/spectral/a8_pipeline/competition_report.py',
               PROJECTS/'comp_assay/analysis/review_20260907/review.md']
    manifest = {
        'acquisition_date':'2026-09-01', 'instrument':'BD FACSDiscover A8',
        'sources':{p.name:{'path':str(p), 'sha256':sha(p)} for p in sources},
        'denominator':'All primary singlets after run-specific Scatter -> SSC Singlets -> FSC Singlets gates.',
        'classifier':source['classifier']['name'], 'unmixing_matrix':source['unmixing_matrix'],
        'nominal_percent':[10,25,50,75,90], 'measurements':20, 'mixture_groups':10,
        'technical_measurements_per_group':2, 'summary':'Arithmetic mean of two fractions, not event-count pooling.',
        'intervals':'Observed minimum to maximum of the two technical measurements; not confidence intervals.',
        'maximum_technical_difference_pp':float(summary.technical_difference_pp.max()),
        'maximum_absolute_mean_deviation_pp':float(summary.deviation_pp.abs().max()),
        'scope':'Saved September classification retained, including low-reporter and Saturated-flagged singlets. No October background exclusion, refitting, inferential tests or calibration correction.',
    }
    (DATA/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    return summary, manifest


def render(summary, manifest):
    apply_style(8.5)
    fig = plt.figure(figsize=(183/25.4, 86/25.4))
    a = fig.add_axes([.105,.185,.345,.655])
    b = fig.add_axes([.625,.185,.345,.655])
    a.plot([0,100],[0,100],color='#A7ADB5',lw=.8,ls=(0,(3,3)),zorder=0)
    b.axhline(0,color='#A7ADB5',lw=.8,ls=(0,(3,3)),zorder=0)
    handles = []
    for orientation,(label,colour,marker,style) in ORIENTATIONS.items():
        d = summary.loc[summary.orientation.eq(orientation)].sort_values('nominal_percent')
        err = np.vstack([d.mean_percent-d.minimum_percent, d.maximum_percent-d.mean_percent])
        for ax,field in [(a,'mean_percent'),(b,'deviation_pp')]:
            ax.errorbar(d.nominal_percent,d[field],yerr=err,color=colour,
                        marker=marker,ms=4,mfc='white',mec=colour,mew=.9,
                        lw=1.1,ls=style,elinewidth=.8,capsize=2.2,capthick=.8,zorder=3)
        handles.append(Line2D([],[],color=colour,marker=marker,mfc='white',ms=4,
                              lw=1.1,ls=style,label=label))
    for ax in (a,b):
        ax.set_xlim(0,100)
        ax.set_xticks([10,25,50,75,90])
        ax.set_xlabel('Nominal LOY fraction (%)',labelpad=5)
        ax.tick_params(length=3,pad=3)
    a.set(ylim=(0,100),yticks=[0,25,50,75,100],ylabel='Measured LOY fraction (%)')
    a.text(.035,.965,'Dashed grey line: identity',transform=a.transAxes,
           ha='left',va='top',fontsize=7.3,color='#6C737C')
    b.set(ylim=(-17,3),yticks=[-15,-10,-5,0],ylabel='Measured − nominal\n(percentage points)')
    b.yaxis.labelpad=6
    worst = summary.loc[summary.deviation_pp.idxmin()]
    b.annotate(f'{worst.deviation_pp:.2f} pp', (worst.nominal_percent,worst.deviation_pp),
               xytext=(-7,-11),textcoords='offset points',ha='right',va='top',
               color=ORIENTATIONS[worst.orientation][1],fontsize=8)
    b.text(.035,.14,f'Max. technical difference\n{manifest["maximum_technical_difference_pp"]:.2f} percentage points',
           transform=b.transAxes,fontsize=7.5,va='bottom',linespacing=1.2)
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.54,.992),ncol=2,
               fontsize=8.2,handlelength=2.2,columnspacing=2.4,borderaxespad=0)
    for x,label in [(.015,'a'),(.515,'b')]:
        fig.text(x,.855,label,fontweight='semibold',fontsize=12,va='top')
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{STEM}.{ext}',dpi=600)
    fig.savefig(REVIEW/f'{STEM}_preview.png',dpi=200)
    plt.close(fig)


def validate(manifest):
    for entry in manifest['sources'].values():
        assert sha(Path(entry['path'])) == entry['sha256']
    pdf = pymupdf.open(OUT/f'{STEM}.pdf')
    assert len(pdf)==1
    page=pdf[0]
    words=page.get_text('words')
    outside=[w for w in words if w[0]<-.5 or w[1]<-.5 or w[2]>page.rect.width+.5 or w[3]>page.rect.height+.5]
    assert not outside,outside
    assert all('OpenSans' in f[3] for f in page.get_fonts())
    page.get_pixmap(matrix=pymupdf.Matrix(2,2),colorspace=pymupdf.csGRAY).save(REVIEW/f'{STEM}_greyscale_check.png')
    svg=ET.parse(OUT/f'{STEM}.svg')
    ntext=len(svg.findall('.//{http://www.w3.org/2000/svg}text'))
    assert ntext>20
    png=Image.open(OUT/f'{STEM}.png')
    assert min(png.info['dpi'])>590
    result={'source_hashes_unchanged':True,'all_20_measurements_preserved':True,
            'ten_means_and_technical_differences_match_existing_thesis_check':True,
            'pdf_text_outside_page':0,'fonts':[f[3] for f in page.get_fonts()],
            'svg_editable_text_elements':ntext,'png_pixels':list(png.size),'png_dpi':list(png.info['dpi']),
            'renderer_sha256':sha(Path(__file__))}
    (DATA/'export_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    pdf.close()


if __name__=='__main__':
    summary,manifest=prepare()
    render(summary,manifest)
    validate(manifest)
