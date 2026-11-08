"""TF regulator and shared-target module figures for HT1080, A549 and DU145, plus a cross-line summary.

Reads only the saved tf_modules_v1 outputs of each RNA reanalysis (rules in each analysis root's
provenance/DECISIONS.md). No statistic is recomputed here; plotted tables are written to source_data.

    ch03_main_ht1080_tf_regulator_modules    (Results 3.4.6.1)
    ch03_supp_a549_tf_regulator_modules      (Results 3.4.6.2)
    ch03_main_du145_tf_regulator_modules     (Results 3.4.6.3)
    ch03_main_crossline_tf_regulators        (Results 3.4.6.4)
"""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from du145_rna_common import FIG, PROJECTS, BLUE, RED, GREY, DARK, MM, panel  # noqa: F401  (applies thesis style)
from crossline_common import LINE_COLOUR

ROOTS = {'HT1080': PROJECTS / 'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE',
         'A549': PROJECTS / 'RNAseq/analysis/A549/SEPT_18_2026_REANALYZE',
         'DU145': PROJECTS / 'RNAseq/analysis/LNCAP/SEPT_23_2026_REANALYZE'}
STEMS = {'HT1080': 'ch03_main_ht1080_tf_regulator_modules', 'A549': 'ch03_supp_a549_tf_regulator_modules',
         'DU145': 'ch03_main_du145_tf_regulator_modules'}
LINES = ['HT1080', 'A549', 'DU145']
# Line-specific wording for the two substituted rules (see each DECISIONS.md).
R2_LABEL = {'HT1080': 'Both periods', 'A549': '8/8 omissions', 'DU145': '3/3 pairs'}
R5_LABEL = {'HT1080': '6/6 pair refits', 'A549': '8/8 sample refits', 'DU145': '3/3 pair refits'}
RULES = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']
RULE_LABEL = {'R1': 'FDR < 0.05', 'R3': '≥3/4 methods', 'R4': 'Alternative priors', 'R6': 'Top targets removed'}
# Module identity: an NPG subset kept distinct from the signed-effect blue/vermilion and the line colours.
MOD_COLOUR = {1: '#4DBBD5', 2: '#91D1C2', 3: '#B09C85', 4: '#F39B7F', 5: '#3C5488'}
PROGRAMME = {'HALLMARK_INTERFERON_ALPHA_RESPONSE': 'IFN-α response', 'HALLMARK_INTERFERON_GAMMA_RESPONSE': 'IFN-γ response',
             'HALLMARK_E2F_TARGETS': 'E2F targets', 'HALLMARK_G2M_CHECKPOINT': 'G2/M checkpoint',
             'HALLMARK_TNFA_SIGNALING_VIA_NFKB': 'TNFα via NF-κB', 'HALLMARK_P53_PATHWAY': 'p53 pathway',
             'HALLMARK_UNFOLDED_PROTEIN_RESPONSE': 'Unfolded protein resp.', 'HALLMARK_HYPOXIA': 'Hypoxia',
             'IEG_vanDenBrink2017_Fig1': 'Dissociation IEGs'}
NULL_LABEL = {'rewired_R1_R2_count': 'Regulators passing R1 and R2', 'rewired_best_silhouette': 'Best mean silhouette',
              'size_matched_mean_overlap': 'Mean target overlap', 'coactivity_within_minus_between': 'Co-activity, within − between'}
PROG_SHORT = {'HALLMARK_INTERFERON_ALPHA_RESPONSE': 'IFN-α', 'HALLMARK_INTERFERON_GAMMA_RESPONSE': 'IFN-γ', 'HALLMARK_E2F_TARGETS': 'E2F',
              'HALLMARK_G2M_CHECKPOINT': 'G2/M', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB': 'TNFα/NF-κB', 'HALLMARK_P53_PATHWAY': 'p53',
              'HALLMARK_UNFOLDED_PROTEIN_RESPONSE': 'UPR', 'HALLMARK_HYPOXIA': 'hypoxia', 'IEG_vanDenBrink2017_Fig1': 'dissociation IEG'}
NULL_REF = {'rewired_R1_R2_count': 'rewired', 'rewired_best_silhouette': 'rewired',
            'size_matched_mean_overlap': 'size-matched', 'coactivity_within_minus_between': 'label-permuted'}
SIGNED = LinearSegmentedColormap.from_list('signed', [BLUE, '#F7F7F7', RED])
FS = 5.6
USED = set()


def tab(line, name):
    p = ROOTS[line] / 'outputs/tf_modules_v1/tables' / f'{name}.csv'
    USED.add(str(p))
    return pd.read_csv(p) if p.exists() else None


def pfmt(p, label='P'):
    if p is None or not np.isfinite(p):
        return f'{label} n.a.'
    if p < 1e-3:
        m, e = f'{p:.1e}'.split('e')
        return f'{label} = {m}×10$^{{{int(e)}}}$'
    return f'{label} = {p:.2g}' if p < .1 else f'{label} = {p:.2f}'


def pshort(p):
    if not np.isfinite(p):
        return 'n.a.'
    return f'{p:.1e}'.replace('e-0', 'e-') if p < 1e-3 else f'{p:.3f}'


def export(fig, stem, tables, notes):
    """Stage then replace PDF/SVG/PNG; write plotted tables and hashed provenance."""
    src = FIG / 'source_data' / stem.removeprefix('ch03_').replace('main_', '').replace('supp_', '')
    src.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        table.to_csv(src / f'{name}.csv', index=False)
    records = []
    for ext in ['pdf', 'svg', 'png']:
        path = FIG / f'{stem}.{ext}'
        stage = FIG / f'{stem}.staging.{ext}'
        fig.savefig(stage, dpi=600)
        stage.replace(path)
        records.append({'file': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    (FIG / 'previews').mkdir(exist_ok=True)
    fig.savefig(FIG / 'previews' / f'{stem}_preview.png', dpi=160)
    plt.close(fig)
    inputs = sorted(notes.pop('inputs'))
    (src / 'provenance.json').write_text(json.dumps({
        'figure': stem, 'font': 'Open Sans', 'script': Path(__file__).name,
        'inputs': [{'path': p, 'sha256': hashlib.sha256(Path(p).read_bytes()).hexdigest()} for p in inputs],
        **notes, 'exports': records}, indent=2), encoding='utf-8')
    print(stem, flush=True)


def regulator_rows(line):
    """Robust and supported regulators: modules in dendrogram order, then unclustered nodes, then non-node regulators."""
    R = tab(line, 'TF_robustness');mem = tab(line, 'TF_module_membership')
    R = R[R.tier.isin(['robust', 'supported'])].copy()
    R = R.merge(mem[['source', 'module', 'order']], on='source', how='left')
    R['node'] = R.source.isin(mem.source)
    R['k1'] = np.where(R.module.notna(), R.module, np.where(R.node, 90, 99))
    R['k2'] = np.where(R.order.notna(), R.order, -R.consensus_z.abs())
    return R.sort_values(['k1', 'k2']).reset_index(drop=True)


def mm_axes(fig, W, H, x, top, w, h):
    """Axes from millimetre coordinates measured from the top-left corner."""
    return fig.add_axes([x / W, 1 - (top + h) / H, w / W, h / H])


def letter(fig, W, H, x, top, s):
    fig.text(x / W, 1 - top / H, s, fontsize=12, weight='semibold', va='top')


def line_figure(line):
    R = regulator_rows(line)
    mem = tab(line, 'TF_module_membership');S = tab(line, 'TF_module_similarity_weighted')
    drv = tab(line, 'TF_module_driver_targets');nul = tab(line, 'TF_module_nulls')
    conf = tab(line, 'TF_module_confounders');uni = tab(line, 'TF_module_target_union_tests')
    scored = len(tab(line, 'TF_robustness'))
    mods = sorted(mem.module.dropna().astype(int).unique())
    n = len(R);small = n <= 10
    # ---- layout (mm): panel heights follow row counts; short regulator lists move panel e to the left column.
    W = 183;pitch = 2.2 if n > 30 else 3.0;ha = max(n, 6) * pitch
    top_a = 30;a_bot = top_a + ha;d_top = a_bot + 22;d_end = d_top + 50
    b_top = 12;b_size = 55
    per_mod = 8 if len(mods) == 1 else 6
    c_rows = sum(min(len(drv[drv.module == m]), per_mod) for m in mods)
    c_top = b_top + b_size + 21;hc = max(c_rows * 2.8, 22);c_bot = c_top + hc
    e_h = len(PROGRAMME) * 3.2
    e_top = d_end + 14 if small else c_bot + 24
    H = max(e_top + e_h + 17, c_bot + 12) if small else max(d_end, e_top + e_h + 17)
    H = round(H + 2)
    fig = plt.figure(figsize=(W * MM, H * MM))

    # ---- a: regulator scores, adjusted P and the rule checks -----------------------------------------------
    y = (n - 1) - np.arange(n);off = max(n, 6) - n
    a1 = mm_axes(fig, W, H, 19, top_a, 25, ha);a2 = mm_axes(fig, W, H, 46, top_a, 9, ha);a3 = mm_axes(fig, W, H, 56, top_a, 33, ha)
    a1.barh(y, R.ulm_score, color=np.where(R.ulm_score > 0, RED, BLUE), height=.68)
    a1.axvline(0, lw=.6, color='#888888')
    a1.set_yticks(y, R.source, fontsize=FS if n <= 30 else 5.0);a1.tick_params(axis='y', length=1.5, pad=1)
    for t_, tier in zip(a1.get_yticklabels(), R.tier):
        t_.set_fontweight('semibold' if tier == 'robust' else 'normal')
    lim = np.ceil(R.ulm_score.abs().max() + .5);a1.set_xlim(-lim, lim);a1.tick_params(axis='x', labelsize=6)
    a1.set_xlabel('ULM score, Y-depleted vs sgCTRL', fontsize=6.0)
    for ax in (a1, a2, a3):
        ax.set_ylim(-.6 - off, n - .4)
    for yy, p in zip(y, R.ulm_FDR):
        a2.text(0, yy, pshort(p), va='center', ha='left', fontsize=5.0, color=DARK)
    a2.axis('off');a2.text(0, n - .1, 'adj. P', fontsize=5.6, weight='semibold', va='bottom')
    xs = {c: i for i, c in enumerate(RULES)};xs['module'] = 6.3
    for c in RULES:
        ok = R[c].astype(bool).to_numpy()
        a3.scatter(np.full(ok.sum(), xs[c]), y[ok], s=13, color=DARK, lw=0)
        a3.scatter(np.full((~ok).sum(), xs[c]), y[~ok], s=13, facecolors='white', edgecolors='#AAAAAA', lw=.6)
    for yy, m, node in zip(y, R.module, R.node):
        fc = MOD_COLOUR[int(m)] if pd.notna(m) else ('#E4E4E4' if node else 'white')
        a3.add_patch(Rectangle((xs['module'] - .35, yy - .38), .7, .76, fc=fc, ec=fc if pd.notna(m) else '#BBBBBB', lw=.5))
        if pd.notna(m):
            a3.text(xs['module'], yy, f'{int(m)}', ha='center', va='center', fontsize=4.8, color=DARK)
    labels = {**RULE_LABEL, 'R2': R2_LABEL[line], 'R5': R5_LABEL[line]}
    for c in RULES:
        a3.text(xs[c], n - .1, f'{c}  {labels[c]}', rotation=90, ha='center', va='bottom', fontsize=5.4)
    a3.text(xs['module'], n - .1, 'Module', rotation=90, ha='center', va='bottom', fontsize=5.4)
    a3.set_xlim(-.6, 6.9);a3.axis('off')
    fig.legend(handles=[Line2D([], [], marker='s', ls='', color=RED, ms=4, label='Higher in Y-depleted'),
                        Line2D([], [], marker='s', ls='', color=BLUE, ms=4, label='Lower in Y-depleted'),
                        Line2D([], [], marker='o', ls='', color=DARK, ms=3.5, label='Check passed'),
                        Line2D([], [], marker='o', ls='', mfc='white', mec='#AAAAAA', ms=3.5, label='Check failed'),
                        Rectangle((0, 0), 1, 1, fc='#E4E4E4', ec='#BBBBBB', lw=.5, label='Node, no module'),
                        Rectangle((0, 0), 1, 1, fc='white', ec='#BBBBBB', lw=.5, label='Not a node')],
               loc='upper left', bbox_to_anchor=(6 / W, 1 - (a_bot + 10) / H), ncol=3, fontsize=5.3, handletextpad=.3,
               columnspacing=.9, labelspacing=.3, borderaxespad=0, handlelength=1)
    tiers = R.tier.value_counts()
    fig.text(11 / W, 1 - 4.6 / H, f'{tiers.get("robust", 0)} robust (bold) and {tiers.get("supported", 0)} supported of '
             f'{scored} regulators scored', fontsize=5.8, color=DARK, va='top')
    letter(fig, W, H, 3, 3, 'a')

    # ---- b: shared-target similarity among module nodes ------------------------------------------------
    order = mem.sort_values('order').source.tolist()
    M = S.set_index('source').loc[order, order].to_numpy()
    b = mm_axes(fig, W, H, 118, b_top, b_size, b_size)
    im = b.imshow(M, cmap='magma_r', vmin=0, vmax=1, interpolation='nearest')
    fsb = 5.4 if len(order) <= 12 else 4.2
    b.set_xticks(range(len(order)), order, rotation=90, fontsize=fsb);b.set_yticks(range(len(order)), order, fontsize=fsb)
    b.tick_params(length=1.5, pad=1)
    mm_ = mem.set_index('source').loc[order, 'module'].to_numpy()
    for m in mods:
        idx = np.where(mm_ == m)[0]
        b.add_patch(Rectangle((idx.min() - .5, idx.min() - .5), len(idx), len(idx), fill=False, ec=MOD_COLOUR[m], lw=1.3))
    for s_ in b.spines.values():
        s_.set_visible(False)
    cb = mm_axes(fig, W, H, 176, b_top + 18, 1.6, 22);c_ = fig.colorbar(im, cax=cb);c_.ax.tick_params(labelsize=5, length=1.5)
    c_.set_label('Weighted overlap\nof leading targets', fontsize=5.4)
    letter(fig, W, H, 100, 3, 'b')
    node_rule = 'robust only' if len(order) == tiers.get('robust', 0) else 'robust plus supported'
    fig.text(107 / W, 1 - 4.6 / H, f'{len(order)} module nodes ({node_rule})', fontsize=5.8, color=DARK, va='top')

    # ---- c: module driver targets -------------------------------------------------------------------------
    c = mm_axes(fig, W, H, 128, c_top, 45, hc)
    rows, labs, seps = [], [], []
    for m in mods:
        d = drv[drv.module == m].sort_values('summed_contribution', ascending=False).head(per_mod)
        rows += d.t.tolist();labs += d.target.tolist();seps.append(len(rows))
    yc = np.arange(len(rows))[::-1]
    c.barh(yc, rows, color=np.where(np.array(rows) > 0, RED, BLUE), height=.68)
    c.axvline(0, lw=.6, color='#888888')
    c.set_yticks(yc, labs, fontsize=5.2);c.set_ylim(-.6, len(rows) - .4);c.tick_params(axis='y', length=1.5, pad=1)
    c.tick_params(axis='x', labelsize=6);c.set_xlabel('Moderated t of driver target', fontsize=6.0)
    start = 0
    for m, end in zip(mods, seps):
        c.add_patch(Rectangle((-.34, 1 - end / len(rows)), .025, (end - start) / len(rows), transform=c.transAxes,
                              fc=MOD_COLOUR[m], ec='none', clip_on=False))
        c.text(-.38, 1 - (start + end) / 2 / len(rows), f'M{m}', transform=c.transAxes, ha='right', va='center', fontsize=5.4, color=DARK)
        if end < len(rows):
            c.axhline(len(rows) - end - .5, lw=.5, color='#CCCCCC')
        start = end
    letter(fig, W, H, 100, c_top - 9, 'c')
    fig.text(107 / W, 1 - (c_top - 7.4) / H, f'Top {per_mod} driver targets per module, by summed contribution', fontsize=5.8, color=DARK, va='top')

    # ---- d: structure nulls -------------------------------------------------------------------------------
    nul = nul.set_index('test')
    tests = [t for t in NULL_LABEL if t in nul.index]
    for i, t in enumerate(tests):
        r = nul.loc[t];ax = mm_axes(fig, W, H, 14 + (i % 2) * 46, d_top + 12 + (i // 2) * 19, 32, 4.5)
        lo, hi, obs = r.null_mean, r.null_q95, r.observed
        ax.plot([lo, hi], [0, 0], color=GREY, lw=3.2, solid_capstyle='butt')
        ax.scatter([lo], [0], s=14, color='#888888', zorder=3)
        ax.scatter([obs], [0], s=26, marker='D', color=DARK, zorder=4)
        span = [v for v in (lo, hi, obs) if np.isfinite(v)];pad = (max(span) - min(span)) * .25 or .1
        ax.set_xlim(min(span) - pad, max(span) + pad);ax.set_ylim(-1, 1);ax.set_yticks([])
        ax.spines['left'].set_visible(False);ax.tick_params(labelsize=5, length=1.5, pad=1)
        pv = r.empirical_p;txt = pfmt(pv, 'empirical P') if np.isfinite(pv) else 'not assessable'
        ax.set_title(NULL_LABEL[t], fontsize=5.6, loc='left', pad=12, weight='normal')
        ax.text(0, 1.3, f'{txt} ({int(r.null_iterations)} {NULL_REF[t]} nulls)', transform=ax.transAxes, ha='left', va='bottom', fontsize=5.0, color=DARK)
    fig.legend(handles=[Line2D([], [], marker='D', ls='', color=DARK, ms=3.5, label='Observed'),
                        Line2D([], [], marker='o', ls='-', color=GREY, mfc='#888888', mec='#888888', lw=3, ms=3.5, label='Null mean to 95th percentile')],
               loc='upper left', bbox_to_anchor=(14 / W, 1 - (d_top + 44) / H), ncol=2, fontsize=5.4, handletextpad=.5, handlelength=2.4,
               columnspacing=1.2, borderaxespad=0)
    letter(fig, W, H, 3, d_top, 'd')

    # ---- e: confounder removal per module ------------------------------------------------------------------
    progs = list(PROGRAMME);ye = np.arange(len(progs))[::-1]
    ex, ew = (34, 60) if small else (124, 52)
    w = (ew - 2 * (len(mods) - 1)) / len(mods)
    for j, m in enumerate(mods):
        ax = mm_axes(fig, W, H, ex + j * (w + 2), e_top, w, e_h)
        cc = conf[conf.module == m].set_index('programme').reindex(progs);u = uni[uni.module == m].iloc[0]
        ax.axvline(u.union_ulm_score, lw=.8, color=MOD_COLOUR[m]);ax.axvline(0, lw=.5, color='#BBBBBB')
        coll = cc.collapses.astype(str).eq('True').to_numpy()
        ax.scatter(cc.removed_score[~coll], ye[~coll], s=14, color=DARK, zorder=3)
        ax.scatter(cc.removed_score[coll], ye[coll], s=16, facecolors='white', edgecolors=RED, lw=.8, zorder=3)
        ax.set_yticks(ye, [PROGRAMME[p] for p in progs] if j == 0 else [], fontsize=5.3);ax.tick_params(axis='y', length=1.5, pad=1)
        ax.set_ylim(-.6, len(progs) - .4);ax.tick_params(axis='x', labelsize=5, length=1.5)
        hi = max(u.union_ulm_score, np.nanmax(cc.removed_score)) * 1.15;ax.set_xlim(min(-.5, np.nanmin(cc.removed_score) - .5), hi)
        ax.set_title(f'Module {m}\nperm. {pfmt(u.perm_p)}', fontsize=5.4, pad=2, weight='normal', color=DARK)
        if j == len(mods) // 2:
            ax.set_xlabel('Module target-union score after removal', fontsize=6.0)
    fig.legend(handles=[Line2D([], [], marker='o', ls='', color=DARK, ms=3.5, label='Same sign, P < 0.05'),
                        Line2D([], [], marker='o', ls='', mfc='white', mec=RED, ms=3.5, label='Collapses'),
                        Line2D([], [], ls='-', color='#777777', lw=.8, label='All targets')],
               loc='upper left', bbox_to_anchor=((ex - 22) / W, 1 - (e_top + e_h + 10) / H), ncol=3, fontsize=5.4, handletextpad=.3,
               columnspacing=.9, borderaxespad=0)
    letter(fig, W, H, 3 if small else 100, e_top - 11, 'e')

    tables = {'regulators': R.drop(columns=['k1', 'k2']), 'similarity_ordered': S.set_index('source').loc[order, order].reset_index(),
              'drivers': drv, 'nulls': nul.reset_index(), 'confounders': conf, 'target_union_tests': uni}
    export(fig, STEMS[line], tables, {'line': line, 'analysis_root': str(ROOTS[line]), 'figure_size_mm': [W, H],
           'rules': str(ROOTS[line] / 'provenance/DECISIONS.md'),
           'display_note': 'fry (estimated inter-gene correlation) is kept in source tables but not displayed, following the 27 Sep 2026 CAMERA reporting preference.',
           'inputs': {p for p in USED if f'{ROOTS[line]}' in p}})


def crossline_figure():
    Rs = {l: tab(l, 'TF_robustness') for l in LINES}
    keep = sorted(set().union(*[set(Rs[l].source[Rs[l].tier == 'robust']) for l in LINES]))
    score = pd.DataFrame({l: Rs[l].set_index('source').ulm_score.reindex(keep) for l in LINES})
    tier = pd.DataFrame({l: Rs[l].set_index('source').tier.reindex(keep) for l in LINES})
    good = tier.isin(['robust', 'supported'])
    key = pd.DataFrame({'n': good.sum(axis=1), 'first': [next((i for i, l in enumerate(LINES) if good.loc[s, l]), 9) for s in keep],
                        'z': [score.loc[s, LINES[next((i for i, l in enumerate(LINES) if good.loc[s, l]), 0)]] for s in keep]}, index=keep)
    order = key.sort_values(['n', 'first', 'z'], ascending=[False, True, True]).index.tolist()
    score, tier = score.loc[order], tier.loc[order]
    fig = plt.figure(figsize=(183 * MM, 170 * MM))

    # ---- a: ULM score for every regulator robust in at least one line ---------------------------------------
    a = fig.add_axes([.14, .06, .16, .87])
    lim = 8;norm = TwoSlopeNorm(0, -lim, lim)
    Z = np.ma.masked_invalid(score.to_numpy())
    cm = SIGNED.copy();cm.set_bad('#FFFFFF')
    im = a.imshow(np.clip(Z, -lim, lim), cmap=cm, norm=norm, aspect='auto', interpolation='nearest')
    for i in range(len(order)):
        for j, l in enumerate(LINES):
            t = tier.iloc[i, j]
            if pd.isna(t):
                a.text(j, i, '–', ha='center', va='center', fontsize=5, color='#999999')
            elif t == 'robust':
                a.scatter(j, i, s=9, color='#222222', lw=0)
            elif t == 'supported':
                a.scatter(j, i, s=9, facecolors='none', edgecolors='#222222', lw=.6)
    a.set_xticks(range(3), LINES, fontsize=6.4);a.xaxis.tick_top()
    for t_, l in zip(a.get_xticklabels(), LINES):
        t_.set_color(LINE_COLOUR[l]);t_.set_fontweight('semibold')
    a.set_yticks(range(len(order)), order, fontsize=5.0);a.tick_params(length=0)
    for s_ in a.spines.values():
        s_.set_visible(False)
    cb = fig.add_axes([.315, .08, .01, .14]);c_ = fig.colorbar(im, cax=cb, extend='both')
    c_.ax.tick_params(labelsize=5, length=1.5);c_.set_label('ULM score, Y-depleted\nvs sgCTRL (clipped ±8)', fontsize=5.4)
    a.legend(handles=[Line2D([], [], marker='o', ls='', color='#222222', ms=3.2, label='Robust'),
                      Line2D([], [], marker='o', ls='', mfc='none', mec='#222222', ms=3.2, label='Supported'),
                      Line2D([], [], marker='$–$', ls='', color='#999999', ms=4, label='Not scored')],
             loc='lower left', bbox_to_anchor=(1.12, .19), fontsize=5.4, handletextpad=.3, labelspacing=.35, borderaxespad=0)
    panel(fig, a, 'a', '', dx=-.12, dy=.02)

    # ---- b: tier counts per line ------------------------------------------------------------------------------
    b = fig.add_axes([.55, .76, .38, .16])
    tiers = ['robust', 'supported', 'consistent_only', 'overall_only']
    tcol = {'robust': '#333333', 'supported': '#8C8C8C', 'consistent_only': '#C8C8C8', 'overall_only': '#E8E8E8'}
    tlab = {'robust': 'Robust', 'supported': 'Supported', 'consistent_only': 'R1 + R2 only', 'overall_only': 'R1 only'}
    yb = np.arange(3)[::-1];left = np.zeros(3)
    counts = pd.DataFrame({l: Rs[l].tier.value_counts().reindex(tiers, fill_value=0) for l in LINES}).T
    for t in tiers:
        b.barh(yb, counts[t], left=left, color=tcol[t], height=.62, label=tlab[t], ec='white', lw=.4);left += counts[t].to_numpy()
    for yy, l, tot in zip(yb, LINES, left):
        b.text(tot + .8, yy, f'{int(tot)} of {len(Rs[l])} scored reach FDR < 0.05', va='center', fontsize=5.4, color=DARK)
    b.set_yticks(yb, LINES, fontsize=6.4)
    for t_, l in zip(b.get_yticklabels(), LINES):
        t_.set_color(LINE_COLOUR[l]);t_.set_fontweight('semibold')
    b.set_xlim(0, left.max() * 1.9);b.set_xlabel('Regulators (overall ULM FDR < 0.05)', fontsize=6.2)
    b.legend(loc='lower left', bbox_to_anchor=(0, 1.03), ncol=4, fontsize=5.4, handlelength=.9, columnspacing=1, borderaxespad=0)
    panel(fig, b, 'b', '', dx=-.08, dy=.035)

    # ---- c: structure null tests per line ---------------------------------------------------------------------
    c = fig.add_axes([.55, .42, .38, .22])
    tests = list(NULL_LABEL);rowsc = []
    for i, t in enumerate(tests):
        for j, l in enumerate(LINES):
            nl = tab(l, 'TF_module_nulls').set_index('test')
            p = nl.empirical_p.get(t, np.nan);rowsc.append({'line': l, 'test': t, 'empirical_p': p,
                                                          'null_iterations': nl.null_iterations.get(t, np.nan)})
            yy = (len(tests) - 1 - i) + (1 - j) * .22
            if np.isfinite(p):
                c.scatter(p, yy, s=18, color=LINE_COLOUR[l], zorder=3)
            else:
                c.text(1.25, yy, 'n.a.', fontsize=5, va='center', color=LINE_COLOUR[l])
    c.set_xscale('log');c.set_xlim(5e-4, 2);c.axvline(.05, lw=.6, color='#999999', ls=(0, (3, 2)))
    c.text(.05, len(tests) - .35, 'P = 0.05', fontsize=5, color='#777777', ha='left', va='bottom')
    c.set_yticks(np.arange(len(tests))[::-1], [NULL_LABEL[t].replace('\n', ' ') for t in tests], fontsize=5.6)
    c.set_ylim(-.6, len(tests) - .2);c.set_xlabel('Empirical P versus null', fontsize=6.2)
    c.legend(handles=[Line2D([], [], marker='o', ls='', color=LINE_COLOUR[l], ms=4, label=l) for l in LINES],
             loc='lower left', fontsize=5.4, handletextpad=.2, borderaxespad=.2)
    panel(fig, c, 'c', '', dx=-.19, dy=.012)

    # ---- d: modules, stability and dependence on confounding programmes -----------------------------------------
    d = fig.add_axes([.55, .06, .22, .26]);rowsd = []
    for l in LINES:
        m = tab(l, 'TF_modules')
        if m is None:
            continue
        for _, r in m.iterrows():
            cd = r.confounder_dependent if isinstance(r.confounder_dependent, str) and r.confounder_dependent else ''
            rowsd.append({'line': l, 'module': int(r.module), 'n': int(r.n), 'direction': r.direction, 'members': r.members,
                          'bootstrap_jaccard': r.bootstrap_jaccard, 'edge_drop_jaccard': r.edge_drop_jaccard,
                          'stable': r.stable, 'confounder_dependent': cd})
    D = pd.DataFrame(rowsd);yd = np.arange(len(D))[::-1]
    d.axvspan(.75, 1.02, color='#F2F2F2', lw=0);d.axvline(.6, lw=.5, color='#BBBBBB', ls=(0, (3, 2)))
    for yy, (_, r) in zip(yd, D.iterrows()):
        d.plot([r.bootstrap_jaccard, r.edge_drop_jaccard], [yy, yy], color='#CCCCCC', lw=.8, zorder=1)
        d.scatter(r.bootstrap_jaccard, yy, s=18, color=LINE_COLOUR[r.line], zorder=3)
        d.scatter(r.edge_drop_jaccard, yy, s=18, marker='s', facecolors='white', edgecolors=LINE_COLOUR[r.line], lw=.8, zorder=3)
        mem_ = r.members.split(';');mtxt = ', '.join(mem_[:4]) + (f' +{len(mem_) - 4}' if len(mem_) > 4 else '')
        conf_ = '\ncollapses without ' + ' or '.join(PROG_SHORT[p] for p in r.confounder_dependent.split(';')) + ' genes' if r.confounder_dependent else ''
        d.text(1.04, yy, f'{mtxt}; {r.direction}{conf_}', va='center', fontsize=4.9, color=DARK, clip_on=False, linespacing=1.15)
    d.set_yticks(yd, [f'{r.line} M{r.module} (n = {r.n})' for _, r in D.iterrows()], fontsize=5.6)
    for t_, l in zip(d.get_yticklabels(), D.line):
        t_.set_color(LINE_COLOUR[l])
    d.set_xlim(.3, 1.02);d.set_ylim(-.6, len(D) - .4);d.set_xlabel('Mean Jaccard with primary module', fontsize=6.2)
    d.text(.885, len(D) - .5, 'stable', fontsize=5, color='#777777', ha='center', va='bottom')
    d.text(1.04, len(D) - .5, 'Members; direction', fontsize=5.2, weight='semibold', va='bottom', clip_on=False)
    d.legend(handles=[Line2D([], [], marker='o', ls='', color='#777777', ms=4, label='Gene bootstrap'),
                      Line2D([], [], marker='s', ls='', mfc='white', mec='#777777', ms=4, label='Edge drop')],
             loc='upper left', fontsize=5.4, handletextpad=.2, borderaxespad=.2)
    panel(fig, d, 'd', '', dx=-.12, dy=.012)

    export(fig, 'ch03_main_crossline_tf_regulators',
           {'regulator_scores': score.reset_index().rename(columns={'index': 'source'}),
            'regulator_tiers': tier.reset_index().rename(columns={'index': 'source'}),
            'tier_counts': counts.reset_index().rename(columns={'index': 'line'}),
            'null_tests': pd.DataFrame(rowsc), 'modules': D},
           {'analysis_roots': {l: str(ROOTS[l]) for l in LINES}, 'inputs': set(USED)})


if __name__ == '__main__':
    for l in LINES:
        USED.clear();line_figure(l)
    USED.clear();crossline_figure()
