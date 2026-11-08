"""Describe parental HT1080 using mean Chronos across 20 reference SAGA genes."""
from pathlib import Path
import hashlib
import json
import pandas as pd

FIG = Path(__file__).resolve().parents[1]
OUT = FIG / 'source_data/depmap_model_positions'
SOURCES = {
    'chronos': Path('A:/DEPMAP/DepMap/CRISPRGeneEffect.csv'),
    'reference': FIG / 'source_data/ht1080_candidate_discovery/saga_reference_scores.csv',
    'cohort': FIG / 'source_data/xy_paralog_ranked_partners_y_status/model_cohort_and_exclusions.csv',
}

def main():
    reference = pd.read_csv(SOURCES['reference'])
    genes = reference.loc[reference.reference_20, 'gene'].tolist()
    assert len(genes) == len(set(genes)) == 20
    header = pd.read_csv(SOURCES['chronos'], nrows=0).columns
    mapping = {c.split(' (')[0]: c for c in header if c.split(' (')[0] in genes}
    assert set(mapping) == set(genes)
    values = pd.read_csv(SOURCES['chronos'], usecols=[header[0], *mapping.values()], index_col=0)
    values = values.rename(columns={v: k for k, v in mapping.items()})[genes]
    assert values.index.is_unique and values.notna().all().all() and len(values) == 1186
    cohort = pd.read_csv(SOURCES['cohort']).set_index('ModelID')
    xy = cohort.index[(cohort.group == 'XY') & cohort.included].intersection(values.index)
    ht = 'ACH-000054'
    assert len(xy) == 324 and ht in xy and cohort.loc[ht, 'CellLineName'] == 'HT-1080'
    scores = values.mean(axis=1)
    model_table = values.copy()
    model_table.insert(0, 'SAGA_mean_Chronos', scores)
    model_table.insert(1, 'Y_retaining_male', values.index.isin(xy))
    model_table.index.name = 'ModelID'
    model_table.to_csv(OUT / 'saga_module_model_scores.csv')
    rows = []
    for name, ids in [('all_models', scores.index), ('Y_retaining_male', xy)]:
        comparison = scores.loc[ids]
        rows.append(dict(cohort=name, n_models=len(comparison), n_genes=20,
                         HT1080_mean_Chronos=scores[ht],
                         HT1080_least_dependent_rank=1 + int((comparison > scores[ht]).sum()),
                         models_more_dependent=int((comparison < scores[ht]).sum())))
    summary = pd.DataFrame(rows)
    assert summary.HT1080_least_dependent_rank.tolist() == [16, 5]
    summary.to_csv(OUT / 'saga_module_summary.csv', index=False)
    manifest = dict(
        definition='Unweighted arithmetic mean of DepMap Public 25Q3 Chronos gene effects across all 20 Herbst et al. 2021 reference SAGA subunits. No selection by screen significance, no gene standardisation, no missing values.',
        ranking='Descending score: rank 1 is least dependent. Descriptive aggregate, not a fitted latent dependency or a test of complex activity.',
        figure='Current parental DepMap-position figure shows individual genes; this aggregate is recorded separately.',
        genes=genes,
        sources={name: dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                 for name, path in SOURCES.items()},
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    (OUT / 'saga_module_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print(summary.to_string(index=False))

if __name__ == '__main__':
    main()
