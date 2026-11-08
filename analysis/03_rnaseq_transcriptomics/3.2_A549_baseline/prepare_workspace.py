from pathlib import Path
import csv, hashlib, json, shutil

root = Path(__file__).resolve().parent
repo = root.parents[2]
reference = repo / 'analysis/HT1080/SEPT_11_2026_REANALYZE'
for folder in ['R', 'config', 'resources', 'provenance/reference_code', 'outputs/tables', 'outputs/models', 'outputs/figures', 'outputs/report', 'outputs/captions', 'outputs/source_data']:
    (root / folder).mkdir(parents=True, exist_ok=True)

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

copied = []
for name in ['collectri_network.csv', 'hg38_cytoBand.txt', 'msigdb_all_collections.rds', 'msigdb_all_set_metadata.csv', 'msigdb_all_collection_catalogue.csv', 'progeny_model_top100.rds']:
    source = reference/'resources'/name
    dest = root/'resources'/name
    if dest.exists():
        assert digest(dest) == digest(source), f'Existing resource differs: {dest}'
    else:
        shutil.copy2(source, dest)
    copied.append({'path': str(dest.relative_to(root)).replace('\\','/'), 'source':str(source), 'sha256':digest(dest)})
with (root/'resources/manifest.csv').open('w', newline='', encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=['path','source','sha256']); w.writeheader(); w.writerows(copied)
shutil.copy2(reference/'config/hg38_PAR.bed', root/'config/hg38_PAR.bed')
for name in ['00_common.R','01_ingest.R','02_effects.R','03_enrichment_helpers.R','03_pathways.R']:
    shutil.copy2(reference/'R'/name, root/'provenance/reference_code'/name)
shutil.copy2(reference/'R/00_common.R', root/'R/00_common.R')
# Ingest is adapted only for sample count and explicitly unknown pairing/derivation.
ingest=(reference/'R/01_ingest.R').read_text(encoding='utf-8-sig')
ingest=ingest.replace("message('All 15 source tables match; importing exact GENCODE v49 gene records')", "message('All 11 source tables match; importing exact GENCODE v49 gene records')")
start=ingest.index('  assert(nrow(core) == 12')
end=ingest.index("  tab(a, 'gene_annotation_gencode49')", start)
ingest=ingest[:start]+"  assert(nrow(core) == 8 && all(table(core$treatment) == 4), 'Expected 4 LOY and 4 CTRL libraries')\n  pairs <- NULL # Author states pairing/derivation history is uncertain.\n  tab(meta, 'sample_manifest')\n"+ingest[end:]
(root/'R/01_ingest.R').write_text(ingest,encoding='utf-8')
effects=(reference/'R/02_effects.R').read_text(encoding='utf-8-sig')
(root/'R/02_normalize.R').write_text(effects.split('expression_analysis <- function')[0],encoding='utf-8')
helpers=(reference/'R/03_enrichment_helpers.R').read_text(encoding='utf-8-sig')
adjust=helpers[helpers.index('enrichment_adjust <- function'):helpers.index('camera_analysis <- function')]
(root/'R/03_adjust.R').write_text(adjust,encoding='utf-8')
metadata=list(csv.DictReader((repo/'analysis/A549/count_matrix/A549_sample_metadata.csv').open(encoding='utf-8-sig')))
rows=[dict(sample_key=z['colname'],treatment=z['treatment'],batch='AUG_2025',sample_index=z['replicate'],pairing_status='uncertain_author_2026-09-18',derivation_status='uncertain_author_2026-09-18',source_count_file=z['file_path']) for z in metadata]
with (root/'config/sample_manifest.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
cfg=dict(seed=42,raw_counts='analysis/A549/count_matrix/A549_counts_raw.csv',gtf='annotation/gencode.v49.primary_assembly.annotation.gtf.gz',sample_manifest='config/sample_manifest.csv',prior_count=2,filter_min_count=10,filter_min_total=15,effect_threshold=0.5,pathway_min_size=15,pathway_max_size=500,fgsea_eps=1e-10,fgsea_nPermSimple=10000,camera_inter_gene_cor=0.01,tf_min_targets=5,normalization='TMM_expressed_autosomes_nonPAR_X',primary_inference='descriptive_unpaired_population_difference',pairing='uncertain_do_not_pair',independent_derivation='uncertain_do_not_claim_replication',count_unit='reads',msigdb_version='2026.1.Hs')
(root/'config/analysis.json').write_text(json.dumps(cfg,indent=2)+'\n',encoding='utf-8')
print(root)
print(f'Prepared {len(rows)} sample records and {len(copied)} frozen resources')
