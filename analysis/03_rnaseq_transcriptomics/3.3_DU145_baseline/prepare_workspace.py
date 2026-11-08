"""Prepare the DU145 (filed as LNCAP) thesis reanalysis workspace.

Copies frozen resources from HT1080 SEPT_11_2026_REANALYZE (hash-verified) and the
unchanged core modules from A549 SEPT_18_2026_REANALYZE, then writes the sample
manifest and configuration. Original count files and earlier analyses are read only.
"""
from pathlib import Path
import csv, hashlib, json, shutil

root = Path(__file__).resolve().parent
repo = root.parents[2]                      # .../projects/RNAseq
ht = repo / 'analysis/HT1080/SEPT_11_2026_REANALYZE'
a549 = repo / 'analysis/A549/SEPT_18_2026_REANALYZE'
for folder in ['R', 'config', 'resources', 'provenance/reference_code', 'outputs/tables', 'outputs/models',
               'outputs/figures', 'outputs/report', 'outputs/captions', 'outputs/source_data']:
    (root / folder).mkdir(parents=True, exist_ok=True)

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

copied = []
for name in ['collectri_network.csv', 'hg38_cytoBand.txt', 'msigdb_all_collections.rds', 'msigdb_all_set_metadata.csv',
             'msigdb_all_collection_catalogue.csv', 'progeny_model_top100.rds']:
    source = ht / 'resources' / name
    dest = root / 'resources' / name
    if dest.exists():
        assert digest(dest) == digest(source), f'Existing resource differs: {dest}'
    else:
        shutil.copy2(source, dest)
    copied.append({'path': f'resources/{name}', 'source': str(source), 'sha256': digest(dest)})
with (root / 'resources/manifest.csv').open('w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['path', 'source', 'sha256']); w.writeheader(); w.writerows(copied)

shutil.copy2(a549 / 'config/hg38_PAR.bed', root / 'config/hg38_PAR.bed')
# Core modules whose logic is unchanged from A549/HT1080.
for name in ['00_common.R', '02_normalize.R', '03_adjust.R']:
    shutil.copy2(a549 / 'R' / name, root / 'R' / name)
for name in ['00_common.R', '01_ingest.R', '02_normalize.R', '03_adjust.R', '04_expression.R', '05_pathways.R',
             '06_figures.R', '07_report.R', '08_par_branch.R']:
    shutil.copy2(a549 / 'R' / name, root / 'provenance/reference_code' / f'A549_{name}')
shutil.copy2(ht / 'R/02_batch_aware_de.R', root / 'provenance/reference_code/HT1080_02_batch_aware_de.R')

counts_dir = repo / 'counts/STAR_v2_PARfix/MAY_2026_LNCAP_LOY'
rows = []
for treatment, label in [('CTRL', 'sgCTRL'), ('LOY', 'LOY')]:
    for p in ['P12', 'P13', 'P14']:
        key = f'LNCaP-{label}-{p}'
        rows.append(dict(sample_key=key, cell_line='DU145', historical_label='LNCaP', treatment=treatment,
                         batch='MAY_2026', passage=p, passage_block=p,
                         population='CP + Q53 + P47 2x' if treatment == 'LOY' else 'sgCTRL',
                         population_source='author_2026-09-23' if treatment == 'LOY' else 'source_label',
                         source_count_file=(counts_dir / f'{key}.txt').as_posix()))
with (root / 'config/sample_manifest.csv').open('w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

cfg = dict(seed=42,
           raw_counts='analysis/LNCAP/count_matrix/LNCAP_counts_raw.csv',
           gtf='annotation/gencode.v49.primary_assembly.annotation.gtf.gz',
           sample_manifest='config/sample_manifest.csv',
           prior_count=2, filter_min_count=10, filter_min_total=15,
           de_fdr=0.05, de_abs_log2fc=1, de_pairs_required=3,
           pathway_min_size=15, pathway_max_size=500, fgsea_eps=1e-10, fgsea_nPermSimple=10000,
           camera_inter_gene_cor=0.01, tf_min_targets=5,
           cn_abs_delta=0.2,
           wgs_bins='../../WGS/LNCaP_QDNAseq_cn_df.rds',
           wgs_ctrl='LNCaP-sgCTRL', wgs_loy='LNCaP-LOY-2X', wgs_loy_alt='LNCaP-DMSO-LOY-1X',
           t2t_genes='../../WGS/T2T_genes.gz',
           t2t_cytobands='../../WGS/chm13v2.0_cytobands_allchrs.bed',
           depmap_root='H:/PRISM/data/CCLE/26Q1',
           depmap_expression='raw/depmap_public_26q1/expression/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv',
           depmap_model='raw/depmap_public_26q1/metadata/Model.csv',
           depmap_top_variable=8000,
           normalization='TMM_expressed_autosomes_nonPAR_X',
           primary_inference='paired_voom_limma_passage_plus_treatment',
           pairing='passage_matched_P12_P13_P14',
           independent_derivation='not_claimed_one_population_per_arm',
           count_unit='reads', msigdb_version='2026.1.Hs')
(root / 'config/analysis.json').write_text(json.dumps(cfg, indent=2) + '\n', encoding='utf-8')
print(root)
print(f'Prepared {len(rows)} sample records and {len(copied)} frozen resources')
