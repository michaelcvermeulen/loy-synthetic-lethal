"""Draft DU145 pathway heatmap matching A549 F03, using saved statistics only."""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from du145_rna_common import read, T, FIG, BLUE, RED, DARK, MM, export

STEM = 'ch03_main_du145_rna_pathways'
GROUPS = [
    ('Interferon and inflammatory signalling', [
        ('HALLMARK_INTERFERON_ALPHA_RESPONSE', '[H] Interferon-alpha response'),
        ('HALLMARK_INTERFERON_GAMMA_RESPONSE', '[H] Interferon-gamma response'),
        ('HALLMARK_INFLAMMATORY_RESPONSE', '[H] Inflammatory response'),
    ]),
    ('Epithelial-associated signatures', [
        ('KOHN_EMT_EPITHELIAL', '[CGP] Kohn: epithelial EMT component'),
        ('CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_UP', '[CGP] Charafe: basal vs mesenchymal, up'),
        ('ONDER_CDH1_TARGETS_2_DN', '[CGP] Onder: CDH1 targets 2, down'),
        ('HOLLERN_EMT_BREAST_TUMOR_DN', '[CGP] Hollern: breast-tumour EMT, down'),
    ]),
    ('Ribosomal and translation-associated sets', [
        ('GOMF_STRUCTURAL_CONSTITUENT_OF_RIBOSOME', '[MF] Structural constituent of ribosome'),
        ('REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION', '[R] Eukaryotic translation elongation'),
        ('REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE', '[R] SRP-dependent protein targeting'),
        ('REACTOME_RRNA_PROCESSING', '[R] rRNA processing'),
    ]),
    ('Histone and chromatin-associated sets', [
        ('GOCC_NUCLEOSOME', '[CC] Nucleosome'),
        ('GOMF_STRUCTURAL_CONSTITUENT_OF_CHROMATIN', '[MF] Structural constituent of chromatin'),
        ('REACTOME_HDACS_DEACETYLATE_HISTONES', '[R] HDACs deacetylate histones'),
        ('REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA', '[R] PRC2 methylates histones and DNA'),
    ]),
    ('Cell-cycle and RHO-associated sets', [
        ('REACTOME_M_PHASE', '[R] M phase'),
        ('REACTOME_RHO_GTPASE_EFFECTORS', '[R] RHO GTPase effectors'),
        ('HALLMARK_G2M_CHECKPOINT', '[H] G2M checkpoint'),
        ('HALLMARK_MITOTIC_SPINDLE', '[H] Mitotic spindle'),
    ]),
]

fg = read('fgsea_all_collections').set_index('pathway', verify_integrity=True)
cam = read('camera_all_collections').set_index('pathway', verify_integrity=True)
assert len(fg) == len(cam) == 23252
assert cam.correlation_mode.eq('fixed_0.01').all()
assert np.allclose(cam.Correlation, .01)
assert cam.residual_df.eq(2).all()
sig = fg.BH_across_collections < .05
assert int(sig.sum()) == 154

records, heat, heads = [], [], []
for theme, members in GROUPS:
    heads.append((len(heat), theme, np.mean([fg.loc[key, 'NES'] for key, _ in members])))
    heat.append([np.nan])
    for key, label in members:
        r, c = fg.loc[key], cam.loc[key]
        records.append(dict(
            row=len(heat), theme=theme, pathway=key, display_label=label,
            NES=r.NES, fgsea_FDR=r.BH_across_collections,
            CAMERA_FDR=c.BH_across_collections, CAMERA_direction=c.Direction,
            gene_set_size=r['size'], leadingEdge=r.leadingEdge,
            collection=r.collection, subcollection=r.subcollection,
            NES_CN_stable_only=r.NES_CN_stable_only,
            FDR_CN_stable_only=r.BH_CN_stable_only,
            survives_CN_screen=r.survives_CN_screen,
        ))
        heat.append([r.NES])
display = pd.DataFrame(records)
assert len(display) == 19 and display.pathway.is_unique
assert display.NES.abs().max() <= 2.6
assert ((display.NES > 0) == display.CAMERA_direction.eq('Up')).all()
assert set(fg[sig & fg.collection.eq('H')].index).issubset(display.pathway)

fig = plt.figure(figsize=(183 * MM, 166 * MM))
ax = fig.add_axes([.53, .15, .105, .65])
cmap = LinearSegmentedColormap.from_list('DU145_NES', [BLUE, '#FAFAFA', RED])
im = ax.imshow(heat, aspect='auto', cmap=cmap,
               norm=TwoSlopeNorm(vmin=-2.6, vcenter=0, vmax=2.6), interpolation='nearest')
ax.set_yticks(display.row, display.display_label, fontsize=7.0)
ax.set_xticks([0], ['Overall'], fontsize=7)
ax.tick_params(length=0)
ax.spines[:].set_visible(False)
for row, heading, direction in heads:
    ax.text(-.57, row, heading, ha='right', va='center', fontsize=7.2,
            weight='semibold', color=BLUE if direction < 0 else '#B5422B', clip_on=False)

stats = fig.add_axes([.66, .15, .31, .65])
stats.set(xlim=(-.5, 1.5), ylim=(len(heat)-.5, -.5))
stats.axis('off')
for r in display.itertuples():
    if r.fgsea_FDR < .05:
        ax.text(0, r.row, '*', ha='center', va='center', fontsize=10,
                color='white' if abs(r.NES) > 2.6*.55 else '#242424')
    for col, value in enumerate([r.fgsea_FDR, r.CAMERA_FDR]):
        stats.text(col, r.row, f'{value:.1e}', ha='center', va='center', fontsize=7,
                   weight='semibold' if value < .05 else 'normal', color=DARK)
for col, heading in enumerate(['fgsea', 'CAMERA\ndefault (0.01)']):
    stats.text(col, -1.0, heading, ha='center', va='bottom', fontsize=7.1)

# Match the panel/header structure of A549 F03; all counts are DU145-specific.
fig.text(.035, .95, 'a', fontsize=12, weight='semibold')
fig.text(.08, .95, 'Gene-set enrichment', fontsize=10, weight='semibold')
fig.text(.66, .95, 'b', fontsize=12, weight='semibold')
fig.text(.70, .95, 'Global FDR', fontsize=9, weight='semibold')
fig.text(.08, .91, f'All collections | {len(fg):,} sets tested', fontsize=7.5)
fig.text(.08, .881,
         f'{sig.sum():,} significant: {(sig & (fg.NES < 0)).sum():,} lower, '
         f'{(sig & (fg.NES > 0)).sum():,} higher | selected sets shown',
         fontsize=7.2, color=DARK)
cax = fig.add_axes([.34, .083, .35, .017])
cb = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=[-2,-1,0,1,2])
cb.ax.tick_params(labelsize=7)
cb.ax.set_title('Normalised enrichment score (NES)', fontsize=7.5, pad=6)
fig.text(.04, .044, '* fgsea FDR < 0.05. Bold values: FDR < 0.05 in the indicated analysis.', fontsize=6.6)
fig.text(.04, .024, 'H, Hallmark; R, Reactome; CGP, chemical and genetic perturbations.', fontsize=6.6)
fig.text(.04, .008, 'MF, GO Molecular Function; CC, GO Cellular Component.', fontsize=6.6)

export(fig, STEM, {'pathway_display': display}, {
    'inputs': [T/'fgsea_all_collections.csv', T/'camera_all_collections.csv'],
    'status': 'Draft matching the design of A549 CH3-S3.4.6.2-F03; no assigned figure ID.',
    'notes': ('Saved effect-ranked fgsea NES and global BH values; saved paired-passage CAMERA '
              'with fixed correlation 0.01 and residual df 2. Both testing families contain '
              '23,252 MSigDB 2026.1.Hs sets. No fit, ranking or correction was repeated. '
              'All three globally significant Hallmarks are included, alongside illustrative '
              'epithelial, ribosomal, histone/chromatin and cell-cycle sets. Hallmark G2M '
              'checkpoint and mitotic spindle are included for context despite nonsignificant '
              'fgsea results. Full source names and leading edges are in the display table. '
              'The shared NES scale matches A549 F03 (-2.6 to +2.6); no displayed score is clipped. '
              'Theme labels describe selected annotations, not direct pathway activity. '
              'Current selected DU145 figures are unchanged.'),
})

src = FIG/'source_data/du145_rna_pathways'
summary = {'sets_tested': len(fg), 'fgsea_significant': int(sig.sum()),
           'negative': int((sig & (fg.NES < 0)).sum()),
           'positive': int((sig & (fg.NES > 0)).sum()),
           'CAMERA_significant': int((cam.BH_across_collections < .05).sum()),
           'displayed_sets': len(display),
           'displayed_fgsea_significant': int((display.fgsea_FDR < .05).sum()),
           'displayed_CAMERA_significant': int((display.CAMERA_FDR < .05).sum()),
           'same_enrichment_direction_between_methods': True}
(src/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(json.dumps(summary))
