"""Package and plot the supplied arm-loss results; no DE/GSEA is rerun.

Run from any directory: python path/to/arm_loss_comparison.py
Inputs remain unchanged. Figure captions live in Chapter3_figure_captions.md.
"""
from pathlib import Path
import hashlib
import io
import json
import shutil
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

THESIS = Path(__file__).resolve().parents[3]
BUNDLE = THESIS.parent / 'RNAseq/analysis/arm_loss_bundle'
FIG = THESIS / 'figures/chapter3'
SRC = FIG / 'source_data/arm_loss_comparison'
TABLE = THESIS / 'tables/chapter3/arm_loss_comparison'
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style

ORDER = ['A2058_1q', 'A2780_1q_vsWT', 'A2780_1q_vsTK', 'AGS_1q_vsTK',
         'MCF10A_1q', 'MCF10Ahr_1q', 'RPTEC_3p_early', 'RPTEC_3p_late',
         'A2058_7p', 'A2058_8q', 'HCT116_8q', 'RPE1_mono1018',
         'RPE1_mono13', 'RPE1_mono19p']
LABELS = ['1q | A2058 [1]', '1q | A2780 vs WT [2]', '1q | A2780 vs TK [2]',
          '1q | AGS vs TK', '1q | MCF10A', '1q | MCF10A-HRAS',
          '3p | RPTEC day 18 [3]', '3p | RPTEC day 90 [3]',
          '7p | A2058 [1]', '8q | A2058 [1]', '8q | HCT116',
          '10+18 | RPE1 p53 KO [4]', '13 | RPE1 p53 KO [4]', '19p | RPE1 p53 KD']
PW = {
    'HALLMARK_P53_PATHWAY': 'p53',
    'HALLMARK_E2F_TARGETS': 'E2F',
    'HALLMARK_G2M_CHECKPOINT': 'G2/M',
    'HALLMARK_MYC_TARGETS_V1': 'MYC V1',
    'HALLMARK_OXIDATIVE_PHOSPHORYLATION': 'OXPHOS',
    'RIBOSOMAL_PROTEIN': 'Ribosomal',
    'HALLMARK_INTERFERON_ALPHA_RESPONSE': 'IFN-α',
    'HALLMARK_INTERFERON_GAMMA_RESPONSE': 'IFN-γ',
    'HALLMARK_TNFA_SIGNALING_VIA_NFKB': 'TNFα–NF-κB',
    'HALLMARK_INFLAMMATORY_RESPONSE': 'Inflammatory',
    'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION': 'EMT',
}
TIME_PW = ['HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
           'HALLMARK_E2F_TARGETS', 'HALLMARK_G2M_CHECKPOINT', 'HALLMARK_MYC_TARGETS_V1',
           'HALLMARK_P53_PATHWAY', 'RIBOSOMAL_PROTEIN']
COL = {'neither': '#C8C8C8', '3p_only': '#00A087', 'Yq_only': '#F39B7F',
       'shared': '#3C5488', 'opposite': '#E64B35'}
MARK = {'neither': '.', '3p_only': '^', 'Yq_only': 's', 'shared': 'o', 'opposite': 'D'}
inputs = {}


def read(rel):
    path = BUNDLE / rel
    inputs[str(rel)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path, sep='\t')


def source(data, name):
    data.to_csv(SRC / name, sep='\t', index=False)


def export(fig, stem):
    for ext in ['pdf', 'svg', 'png']:
        buffer = io.BytesIO()
        fig.savefig(buffer, format=ext, dpi=300)
        (FIG / f'{stem}.{ext}').write_bytes(buffer.getvalue())
    plt.close(fig)


def title(ax, letter, text):
    ax.set_title(text, loc='left', pad=6)
    ax.text(-.13, 1.035, letter, transform=ax.transAxes, weight='semibold', fontsize=11)


def scatter(ax, d, letter, heading, rho):
    ax.axhline(0, color='.84', lw=.5)
    ax.axvline(0, color='.84', lw=.5)
    ax.plot([-3.5, 3.5], [-3.5, 3.5], '--', color='.75', lw=.5)
    for cls in COL:
        q = d[d['class'] == cls]
        ax.scatter(q.NES_3p, q.NES_Yq, c=COL[cls], marker=MARK[cls],
                   s=12 if len(d) < 100 else (15 if cls == 'opposite' else 5),
                   linewidths=0, alpha=.65 if cls == 'neither' else .95)
    ax.set(xlim=(-3.5, 3.5), ylim=(-3.5, 3.5), xticks=[-3, 0, 3], yticks=[-3, 0, 3],
           xlabel='NES, 3p loss vs intact', ylabel='NES, Yq loss vs intact')
    ax.set_aspect('equal')
    ax.text(.98, .02, f'ρ = {rho:.2f}; {len(d):,} sets', ha='right', va='bottom',
            transform=ax.transAxes, fontsize=7.5)
    title(ax, letter, heading)


def annotate(ax, d, name, label, xytext):
    r = d.set_index('pathway').loc[name]
    ax.annotate(label, (r.NES_3p, r.NES_Yq), xytext=xytext,
                fontsize=7.3, ha='left', va='center',
                arrowprops={'arrowstyle': '-', 'lw': .45, 'color': '.4'},
                bbox={'facecolor': 'white', 'alpha': .9, 'edgecolor': 'none', 'pad': .8})


def prepare():
    SRC.mkdir(parents=True, exist_ok=True)
    TABLE.mkdir(parents=True, exist_ok=True)
    a = read('figures/source_data/main_gse336056_3p_vs_yq_nes_a.tsv')
    b = read('figures/source_data/main_gse336056_3p_vs_yq_nes_b.tsv')
    c = read('figures/source_data/main_arm_loss_meta_a.tsv')
    d = read('figures/source_data/main_arm_loss_meta_d.tsv')
    qc = read('tables/arm_loss_meta/contrast_qc.tsv').query("type == 'arm_loss'")
    qc = qc.set_index('contrast').loc[ORDER].reset_index()
    assert len(a) == 49 and len(b) == 3178 and len(c) == 154 and len(qc) == 14
    assert set(c.pathway) == set(PW) and set(c.contrast) == set(ORDER)
    # Check plotted values against the complete numerical results.
    for arm, prefix in [('Yq', 'sgYq_YqNull'), ('3p', 'sg3p_3pLoss')]:
        full = read(f'tables/gse336056/GSE336056_{prefix}_vs_intact_GSEA.tsv.gz').set_index('pathway')
        for panel in [a, b]:
            np.testing.assert_allclose(panel[f'NES_{arm}'], full.loc[panel.pathway, 'NES'])
            np.testing.assert_allclose(panel[f'padj_{arm}'], full.loc[panel.pathway, 'padj'])
    for contrast in ORDER:
        full = read(f'tables/arm_loss_meta/gsea/{contrast}_GSEA.tsv.gz').set_index('pathway')
        for panel in [c, d]:
            sub = panel[panel.contrast == contrast]
            np.testing.assert_allclose(sub.NES, full.loc[sub.pathway, 'NES'])
            np.testing.assert_allclose(sub.padj, full.loc[sub.pathway, 'padj'])
        for folder, suffix in [('de', 'DE'), ('gsea', 'GSEA')]:
            rel = f'tables/arm_loss_meta/{folder}/{contrast}_{suffix}.tsv.gz'
            path = BUNDLE / rel
            inputs[rel] = hashlib.sha256(path.read_bytes()).hexdigest()
            shutil.copyfile(path, TABLE / path.name)
    for path in (BUNDLE / 'tables/gse336056').glob('*_vs_intact_*.tsv.gz'):
        inputs[str(path.relative_to(BUNDLE))] = hashlib.sha256(path.read_bytes()).hexdigest()
        shutil.copyfile(path, TABLE / path.name)
    for data, name in [(a, 'main_a_hallmark.tsv'), (b, 'main_b_go.tsv'),
                       (c, 'main_c_programmes.tsv'), (d, 'main_d_time.tsv'), (qc, 'supp_d_bulk_qc.tsv')]:
        source(data, name)
    qc.to_csv(TABLE / 'contrast_design_qc.tsv', sep='\t', index=False)
    return a, b, c, d, qc


def main_figure(a, b, c, d):
    fig = plt.figure(figsize=(7.1, 8.6))
    axa = fig.add_axes([.10, .735, .36, .235])
    axb = fig.add_axes([.59, .735, .36, .235])
    scatter(axa, a, 'a', 'Early response: Hallmark', .71390)
    scatter(axb, b, 'b', 'Early response: GO BP', .49870)
    annotate(axa, a, 'HALLMARK_P53_PATHWAY', 'p53', (-.1, 2.9))
    annotate(axa, a, 'HALLMARK_E2F_TARGETS', 'E2F', (-3.25, -1.7))
    annotate(axa, a, 'HALLMARK_INTERFERON_ALPHA_RESPONSE', 'IFN-α', (.1, .5))
    interferon = b.set_index('pathway').loc['GOBP_RESPONSE_TO_TYPE_I_INTERFERON']
    axb.scatter(interferon.NES_3p, interferon.NES_Yq, marker='D', s=27,
                facecolor=COL['opposite'], edgecolor='#222222', linewidth=.65, zorder=6)
    axb.annotate('Response to\ntype I interferon',
                 xy=(interferon.NES_3p, interferon.NES_Yq), xytext=(-3.15, .12),
                 fontsize=7.5, ha='left', va='center', zorder=7,
                 arrowprops={'arrowstyle': '-|>', 'lw': .85, 'color': '#222222',
                             'shrinkA': 4, 'shrinkB': 4, 'mutation_scale': 8},
                 bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 1.3})
    handles = [Line2D([], [], color=COL[k], marker=MARK[k], lw=0, markersize=4,
                       label=label) for k, label in [
                           ('shared', 'Both, same direction'), ('opposite', 'Both, opposite'),
                           ('3p_only', '3p only'), ('Yq_only', 'Yq only'), ('neither', 'Neither')]]
    fig.legend(handles=handles, loc='center', bbox_to_anchor=(.62, .668),
               ncol=3, fontsize=7.2,
               columnspacing=1.2, handletextpad=.3)
    fig.text(.10, .666, 'FDR < 0.05 in:', fontsize=7.5, va='center')
    axc = fig.add_axes([.355, .335, .565, .285])
    mat = c.pivot(index='contrast', columns='pathway', values='NES').loc[ORDER, list(PW)]
    fdr = c.pivot(index='contrast', columns='pathway', values='padj').loc[ORDER, list(PW)]
    cmap = LinearSegmentedColormap.from_list('signed', ['#3C5488', '#FFFFFF', '#E64B35'])
    assert np.isfinite(mat.to_numpy()).all() and mat.abs().max().max() <= 4
    im = axc.imshow(mat, cmap=cmap, vmin=-4, vmax=4, aspect='auto', interpolation='none')
    yy, xx = np.where(fdr.to_numpy() < .05)
    axc.scatter(xx, yy, c='#222222', s=5)
    axc.set_yticks(range(14), LABELS, fontsize=7.6)
    axc.set_xticks(range(11), list(PW.values()), rotation=48, ha='right', fontsize=7.5)
    axc.tick_params(length=0)
    for sp in axc.spines.values(): sp.set_visible(False)
    title(axc, 'c', 'Established models: selected programmes')
    cbax = fig.add_axes([.936, .405, .012, .14])
    cb = fig.colorbar(im, cax=cbax, ticks=[-4, 0, 4]); cb.ax.tick_params(labelsize=7)
    cb.ax.set_title('NES', fontsize=7)
    axd = fig.add_axes([.355, .058, .585, .165])
    early = d.query("contrast == 'RPTEC_3p_early'").set_index('pathway').loc[TIME_PW]
    late = d.query("contrast == 'RPTEC_3p_late'").set_index('pathway').loc[TIME_PW]
    for i, pw in enumerate(TIME_PW):
        axd.plot([early.loc[pw, 'NES'], late.loc[pw, 'NES']], [i, i], color='.65', lw=.7)
        for frame, colour, marker in [(early, '#00A087', 'o'), (late, '#7E6148', 's')]:
            row = frame.loc[pw]
            axd.scatter(row.NES, i, marker=marker, s=24,
                        facecolor=colour if row.padj < .05 else 'white', edgecolor=colour, lw=.8, zorder=3)
    axd.axvline(0, color='.8', lw=.5)
    axd.set(yticks=range(7), yticklabels=[PW[p] for p in TIME_PW], ylim=(6.6, -.6),
            xlim=(-3.2, 3.2), xticks=[-3, -2, -1, 0, 1, 2, 3], xlabel='NES relative to unedited RPTEC')
    axd.tick_params(axis='y', length=0)
    title(axd, 'd', 'RPTEC 3p-loss clones: day 18 and day 90')
    fig.legend(handles=[Line2D([], [], marker='o', color='#00A087', lw=0, label='Day 18'),
                        Line2D([], [], marker='s', color='#7E6148', lw=0, label='Day 90')],
               loc='center', bbox_to_anchor=(.16, .13), ncol=1, fontsize=8, markerscale=.75)
    fig.text(.355, .005, 'c, dot: FDR < 0.05. d, filled: FDR < 0.05; open: FDR ≥ 0.05.', fontsize=7)
    export(fig, 'ch03_main_arm_loss_comparison')


def supplementary(qc):
    y = read('figures/source_data/main_gse336056_yq_breakpoint_b.tsv')
    y = y.loc[~y.gene.isin(['TTTY14', 'TTTY10'])].reset_index(drop=True)
    cells = read('tables/gse336056/GSE336056_chr3p_dosage_per_cell.tsv.gz').query("arm == 'sg3p'").copy()
    bins = read('tables/gse336056/GSE336056_chr3_breakpoint_bins.tsv')
    # Preserve exact calls and denominators; plotted classes come from supplied tables.
    call = (cells.chr3p_l2r < np.log2(.8)) & (cells.chr3q_l2r.abs() < .25)
    assert np.array_equal(call, cells.loss3p)
    assert len(cells) == 9583 and int(call.sum()) == 649
    for data, name in [(y, 'supp_a_y_regional_detection.tsv'), (cells, 'supp_b_3p_calls.tsv.gz'),
                       (bins, 'supp_c_chr3_bins.tsv')]: source(data, name)
    fig = plt.figure(figsize=(7.1, 8.6))
    axa = fig.add_axes([.135, .625, .34, .33])
    positions = np.arange(len(y))
    axa.barh(positions - .19, y.det_int, height=.36, color='#3C5488', label='Yq-intact')
    axa.barh(positions + .19, y.det_null, height=.36, color='#E64B35', label='Yq-null')
    labels = [r.gene + (' *' if r.selection else '') for r in y.itertuples()]
    axa.set(yticks=range(len(y)), yticklabels=labels, ylim=(len(y)-.4, -.6), xlim=(0, 1.05),
            xticks=[0, .5, 1], xlabel='Fraction of cells with detected RNA')
    axa.tick_params(axis='y', length=0, labelsize=7.5)
    axa.axhline(4.5, color='.8', lw=.5)
    title(axa, 'a', '')
    axa.legend(handles=[Patch(facecolor='#3C5488', label='Yq-intact'),
                        Patch(facecolor='#E64B35', label='Yq-null')],
               loc='lower center', bbox_to_anchor=(.5, -.23), ncol=2, fontsize=7, markerscale=.7)
    axb = fig.add_axes([.64, .555, .315, .165])
    for status, colour in [(False, '#B8B8B8'), (True, '#E64B35')]:
        sub = cells[cells.loss3p == status]
        axb.scatter(sub.chr3p_l2r, sub.chr3q_l2r, s=2, color=colour, alpha=.45,
                    rasterized=True)
    axb.axvline(np.log2(.8), ls='--', color='.3', lw=.6)
    for v in [-.25, .25]: axb.axhline(v, ls='--', color='.3', lw=.6)
    axb.set(xlabel='3p RNA log₂ ratio', ylabel='3q RNA log₂ ratio')
    title(axb, 'c', '')
    axb.text(.98, .04, '649 loss / 8,934 other cells', transform=axb.transAxes, ha='right', fontsize=7)
    axc = fig.add_axes([.135, .415, .34, .12])
    axc.plot(bins.bin_mb+2.5, bins.l2r_loss_vs_ntc, marker='o', ms=2, lw=.6, color='#3C5488')
    axc.axvline(90.5, color='.5', ls='--', lw=.6)
    axc.axhline(0, color='.8', lw=.5)
    axc.set(xlabel='Chromosome 3 position (Mb)', ylabel='RNA log₂ ratio', xticks=[0, 100, 200])
    title(axc, 'd', '')
    axd = fig.add_axes([.355, .09, .58, .22])
    for i, row in qc.iterrows():
        axd.plot([row.median_lfc_cis, row.median_lfc_trans], [i, i], color='.75', lw=.6)
        axd.scatter(row.median_lfc_cis, i, c='#3C5488', marker='s', s=18)
        axd.scatter(row.median_lfc_trans, i, c='#00A087', marker='o', s=18)
    axd.axvline(0, color='.75', lw=.5)
    axd.set(yticks=range(14), yticklabels=LABELS, ylim=(13.7, -.7),
            xlabel='Median gene log₂ fold change (loss vs control)', xlim=(-1, .2))
    axd.tick_params(axis='y', length=0, labelsize=7.5)
    title(axd, 'f', '')
    axd.legend(handles=[Line2D([], [], marker='s', color='#3C5488', lw=0, label='Within lost region'),
                        Line2D([], [], marker='o', color='#00A087', lw=0, label='Outside lost region')],
               loc='upper center', bbox_to_anchor=(.5, -.27), ncol=2, fontsize=7, markerscale=.7)
    # Compact Y classification, reconstructed raw panel UMI counts from saved values.
    yc = read('figures/source_data/main_gse336056_yq_breakpoint_c.tsv')
    yc = yc.loc[yc.arm_grp.eq('sgYq')].copy()
    yc['Yp_UMI'] = np.rint(yc.yp_per10k * yc.nUMI / 10000).astype(int)
    yc['Yq_UMI'] = np.rint(yc.yq_per10k * yc.nUMI / 10000).astype(int)
    expected = np.where((yc.Yp_UMI >= 5) & (yc.Yq_UMI == 0), 'Yq_null',
                        np.where((yc.Yp_UMI >= 5) & (yc.Yq_UMI >= 2), 'Yq_intact', 'ambiguous'))
    assert np.array_equal(expected, yc.status)
    assert len(yc) == 7068 and yc.status.eq('Yq_null').sum() == 194
    source(yc, 'supp_b_y_classification.tsv.gz')
    ay = fig.add_axes([.64, .795, .315, .16])
    for status, colour in [('Yq_intact', '#3C5488'), ('ambiguous', '#BDBDBD'), ('Yq_null', '#E64B35')]:
        q = yc.loc[yc.status.eq(status)]
        ay.scatter(np.log10(q.Yp_UMI + 1), np.log10(q.Yq_UMI + 1),
                   s=3, color=colour, alpha=.45, linewidths=0, rasterized=True)
    ay.axvline(np.log10(6), color='.35', ls='--', lw=.6)
    ay.axhline(np.log10(3), color='.35', ls='--', lw=.6)
    ay.set(xlabel='Yp-panel UMI count', ylabel='Yq-panel UMI count',
           xlim=(-.04, 2.1), ylim=(-.08, 2.4))
    ay.set_xticks(np.log10(np.array([0,5,20,100])+1), ['0','5','20','100'])
    ay.set_yticks(np.log10(np.array([0,2,20,200])+1), ['0','2','20','200'])
    title(ay, 'b', '')
    # Denominators include all depth-filtered cells within the corresponding sgRNA group.
    counts = pd.DataFrame([
        {'population':'sgYq', 'loss_class':'Yq-null', 'loss_cells':194, 'total_cells':7068},
        {'population':'sg3p', 'loss_class':'3p loss', 'loss_cells':649, 'total_cells':9583}])
    counts['percent'] = 100 * counts.loss_cells / counts.total_cells
    source(counts, 'supp_e_population_counts.tsv')
    ac = fig.add_axes([.68, .415, .275, .075])
    ac.barh([1,0], counts.percent, color='#E64B35', height=.52)
    for yy, row in zip([1,0], counts.itertuples()):
        ac.text(row.percent+.15, yy, f'{row.loss_cells:,}/{row.total_cells:,} ({row.percent:.1f}%)',
                va='center', fontsize=6.5)
    ac.set(yticks=[1,0], yticklabels=['Yq-null','3p loss'], xlim=(0,14), xticks=[0,5,10],
           xlabel='Cells in respective sgRNA group (%)', ylim=(-.6,1.6))
    ac.tick_params(axis='y',length=0)
    title(ac, 'e', '')
    export(fig, 'ch03_supp_arm_loss_rna_checks')


def design_tables(qc):
    cd = read('tables/arm_loss_meta/coldata_all.tsv')
    rows = []
    specifications = {}
    for contrast in ORDER:
        if contrast.startswith('RPTEC'):
            sub = cd[cd.dataset == 'GSE336576']
            case = sub[(sub.condition == 'loss') & (sub.passage == ('Early' if contrast.endswith('early') else 'Late'))]
            ctrl = sub[sub.condition == 'WT']
            unit = 'Three clones; repeat libraries summed within clone/time'
            control_unit = 'Three libraries from unedited control'
        elif contrast.startswith('RPE1'):
            base, parent = {'RPE1_mono1018': ('RM_1018', 'RPE1_p53KO'),
                            'RPE1_mono13': ('RM_13', 'RPE1_p53KO'),
                            'RPE1_mono19p': ('RM_19p', 'RPE1_p53_KD')}[contrast]
            sub = cd[cd.dataset == 'GSE150686']
            case = sub[sub.base == base]; ctrl = sub[sub.base == parent]
            unit = 'Three samples from one established monosomic line'
            control_unit = 'Three samples from matched p53-deficient parental line'
        else:
            line, seg = contrast.split('_')[:2]
            if line == 'MCF10Ahr': line = 'MCF10A_hRas'
            seg = seg.split('_')[0]
            sub = cd[(cd.dataset == 'GSE222379') & (cd.cell_line == line)]
            case = sub[(sub.condition == 'loss') & (sub.arm_lost == seg)]
            ctrl = sub[sub.condition == ('TK' if contrast.endswith('vsTK') else 'WT')]
            unit = 'Loss clones; repeat libraries summed within clone'
            control_unit = 'Two TK control clones' if line == 'AGS' else 'Separate control libraries retained'
        ncase, nctrl = case.unit.nunique(), ctrl.unit.nunique()
        qr = qc.set_index('contrast').loc[contrast]
        assert (ncase, nctrl) == (qr.n_case, qr.n_ctrl), contrast
        specifications[contrast] = {'case': list(case['sample']), 'control': list(ctrl['sample'])}
        for group, frame in [('case', case), ('control', ctrl)]:
            for _, r in frame.iterrows():
                rows.append({'contrast': contrast, 'group': group, **r.to_dict()})
        qc.loc[qc.contrast == contrast, 'case_unit_description'] = unit
        qc.loc[qc.contrast == contrast, 'control_unit_description'] = control_unit
        qc.loc[qc.contrast == contrast, 'case_libraries'] = len(case)
        qc.loc[qc.contrast == contrast, 'control_libraries'] = len(ctrl)
    pd.DataFrame(rows).to_csv(TABLE / 'sample_contrast_membership.tsv', sep='\t', index=False)
    qc.to_csv(TABLE / 'contrast_design_qc.tsv', sep='\t', index=False)
    md = ['# Chromosome-loss comparison: data and analysis tables', '',
          'Four public datasets; two single-cell contrasts and 14 bulk contrasts. These are additional to the 18 external LOY contrasts in Section 3.4.6.4.', '',
          '| Accession | Included comparison | Sampling |',
          '|---|---|---|',
          '| [GSE336056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336056) | Two within-guide, RNA-defined Yq/3p contrasts | Days 5 and 10; four HTO samples per guide; one RPTEC background |',
          '| [GSE336576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336576) | Two 3p-loss-versus-unedited contrasts | Three clones at days 18 and 90; 21 libraries |',
          '| [GSE222379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222379) | Nine ReDACT extra-arm-removal contrasts | Six backgrounds; 61 included libraries |',
          '| [GSE150686](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150686) | Three RPE1 monosomy contrasts | p53-deficient loss/control lines; 15 included libraries |', '',
          '| Bulk contrast | Case units | Control units | Case unit definition | Control unit definition |',
          '|---|---:|---:|---|---|']
    for r in qc.itertuples():
        md.append(f'| {r.contrast} | {r.n_case} | {r.n_ctrl} | {r.case_unit_description} | {r.control_unit_description} |')
    md += ['', 'Shared samples: A2058 contrasts reuse WT controls; A2780 contrasts reuse loss clones; RPTEC contrasts reuse controls and sample the same clones at two times; RPE1 10+18 and 13 reuse p53-KO controls.', '',
           '[Dataset provenance and Methods checks](../../../research/molecular_characterisation/arm_loss_integration_notes_2026-09-28.md).', '',
           '[Sample-to-contrast membership](sample_contrast_membership.tsv) records every included library and its aggregation unit. [Design and regional RNA QC](contrast_design_qc.tsv) gives library counts, model-unit counts and regional expression summaries.', '',
           'Complete result tables (unchanged copies from the supplied bundle):', '']
    md += [f'- [{p.name}]({p.name})' for p in sorted(TABLE.glob('*.tsv.gz'))]
    (TABLE / 'README.md').write_text('\n'.join(md)+'\n', encoding='utf-8')


def run():
    font = apply_style(8)
    a, b, c, d, qc = prepare()
    main_figure(a, b, c, d)
    supplementary(qc)
    design_tables(qc)
    manifest = {'source_bundle': str(BUNDLE), 'font': font, 'figure_size_inches': [7.1, 8.6],
                'input_sha256': inputs, 'analysis_rerun': False,
                'checks': ['49 Hallmark / 3178 GO paired sets', '14 bulk contrasts / 11 selected programmes',
                           'All main-panel NES and FDR values match complete GSEA tables',
                           '649/9583 sg3p RNA calls reproduce saved rule',
                           'All bulk model-unit counts agree with sample aggregation']}
    (SRC / 'validation.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'figures': 2, 'formats_each': ['pdf', 'svg', 'png'],
                      'complete_result_tables': len(list(TABLE.glob('*.tsv.gz'))),
                      'validated_bulk_contrasts': len(qc), 'font': font}))


if __name__ == '__main__':
    # Historical packaging functions above are retained for provenance.
    # The default entry point renders only the current concise main figure.
    import runpy
    runpy.run_path(str(Path(__file__).with_name('make_arm_loss_concise.py')), run_name='__main__')
