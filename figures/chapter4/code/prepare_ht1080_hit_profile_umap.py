"""Embed genome-wide correlation fingerprints defined by 56 HT1080 candidates.

Frozen pairwise Pearson correlations are reused, not refitted. The diagonal
self-comparisons are omitted from distances (no value imputation). Exact kNN
distances are Euclidean, rescaled by 56 / number of shared non-self features.
"""
from pathlib import Path
import argparse, hashlib, importlib.metadata, json, os, sys, time

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_hit_profile_umap'
SOURCE = OUT / 'source_data/ht1080_depmap_codependency'
os.environ.setdefault('NUMBA_CACHE_DIR', str(DATA / 'numba_cache'))
os.environ.setdefault('NUMBA_NUM_THREADS', '4')
import numpy as np
import pandas as pd


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare():
    DATA.mkdir(exist_ok=True)
    annot = pd.read_csv(SOURCE / 'seed_annotations.csv').set_index('gene')
    hits = annot.loc[annot.screen_hit].copy()
    features = hits.loc[hits.depmap_available].copy()
    assert len(hits) == 57 and len(features) == 56
    assert hits.screen_direction.value_counts().to_dict() == {'Depleted': 48, 'Enriched': 9}
    assert features.depmap_n.eq(1186).all()
    assert hits.loc[~hits.depmap_available].index.tolist() == ['ZCCHC16']
    z = np.load(SOURCE / 'genomewide_correlations.npz')
    idx = [z['seeds'].tolist().index(g) for g in features.index]
    raw = z['r'][idx].T.copy()
    n = z['n'][idx].T.copy()
    eligible = np.isfinite(raw).all(axis=1)
    assert eligible.sum() == 17931
    genes = z['genes'][eligible]
    raw, n = raw[eligible], n[eligible]
    profiles = raw.copy()
    for j, g in enumerate(features.index):
        i = genes.tolist().index(g)
        assert abs(raw[i, j] - 1) < 1e-9
        profiles[i, j] = np.nan
    assert np.isnan(profiles).sum() == 56
    assert np.isfinite(profiles).sum(1).min() == 55
    np.savez_compressed(DATA / 'input_profiles.npz', genes=genes,
                        features=features.index.to_numpy(dtype=str), raw_r=raw,
                        distance_r=profiles, n_models=n)
    features.to_csv(DATA / 'feature_candidates.csv', index_label='gene')
    pd.DataFrame(profiles, index=genes, columns=features.index).to_csv(
        DATA / 'gene_correlation_profiles.csv.gz', index_label='gene')
    pd.DataFrame(n, index=genes, columns=features.index).to_csv(
        DATA / 'correlation_model_counts.csv.gz', index_label='gene')
    excluded = pd.DataFrame({'gene': z['genes'][~eligible],
                            'available_correlations': np.isfinite(z['r'][idx].T[~eligible]).sum(1),
                            'min_shared_models': z['n'][idx].T[~eligible].min(1),
                            'max_shared_models': z['n'][idx].T[~eligible].max(1)})
    excluded['reason'] = np.where(excluded.max_shared_models < 200,
                                  'Fewer than 200 shared models',
                                  'Undefined correlation (see original source)')
    excluded.to_csv(DATA / 'excluded_genes.csv', index=False)
    tab = pd.DataFrame(index=pd.Index(genes, name='gene'))
    tab['screen_direction'] = hits.screen_direction.reindex(genes).fillna('Not nominated')
    tab['candidate_normZ'] = hits.normZ.reindex(genes)
    tab['candidate_fdr'] = hits.screen_fdr.reindex(genes)
    tab['saga_module'] = annot.saga_module.reindex(genes).fillna('')
    tab['saga_reference_20'] = annot.saga_reference_20.reindex(genes).fillna(False).astype(bool)
    xy = pd.read_csv(OUT / 'source_data/xy_paralog_codependency/reference_gametolog_pairs.csv')
    tab['gametolog'] = xy.drop_duplicates('common_name').set_index('common_name').Gametolog.reindex(genes).fillna('')
    tab['n_models'] = n.min(1)
    tab['profile_rms'] = np.sqrt(np.nanmean(profiles ** 2, axis=1))
    tab['mean_r'] = np.nanmean(profiles, axis=1)
    tab['max_positive_r'] = np.nanmax(profiles, axis=1)
    tab['strongest_positive_hit'] = features.index.to_numpy()[np.nanargmax(profiles, axis=1)]
    tab.to_csv(DATA / 'gene_annotations.csv')
    manifest = dict(
        source_release='DepMap Public 25Q3', source_models=1186,
        candidate_selection='Original directional drugZ FDR <0.25: 48 depleted + 9 enriched',
        unavailable_candidate='ZCCHC16', feature_count=56,
        feature_directions=features.screen_direction.value_counts().to_dict(),
        source_genes=len(z['genes']), embedded_genes=len(genes),
        excluded_genes=len(excluded), exclusion_reasons=excluded.reason.value_counts().to_dict(),
        correlations='Frozen unadjusted pairwise-complete Pearson; minimum 200 models',
        n_models_range=[int(n.min()), int(n.max())],
        cohort='All available models; no sex, lineage or Y-status restriction',
        feature_scaling='None; raw signed r values, no PCA and no top-partner cutoff',
        self_handling='Each candidate self-comparison is missing for distance calculation, never set to zero or imputed',
        metric='sqrt(56 / shared_feature_count * sum((r_i - r_j)^2)), excluding self entries of either gene',
        shared_feature_range=[54, 56],
        interpretation='A genome-wide map conditioned on HT1080-selected candidates, not independent validation of those candidates or a general DepMap atlas',
        sources={str(p): digest(p) for p in [SOURCE / 'genomewide_correlations.npz',
                                            SOURCE / 'seed_annotations.csv',
                                            SOURCE / 'analysis_manifest.json']})
    (DATA / 'analysis_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(manifest, indent=2), flush=True)


def exact_knn():
    z = np.load(DATA / 'input_profiles.npz')
    profiles = z['distance_r']
    # Float64 arithmetic; no full genome x genome distance matrix is retained.
    mask = np.isfinite(profiles).astype(float)
    x = np.nan_to_num(profiles)
    x2 = x * x
    m, p = x.shape
    k = 51  # includes self; supplies all planned 15/30/50-neighbour fits
    indices = np.empty((m, k), dtype=np.int32)
    distances = np.empty((m, k), dtype=np.float32)
    t = time.monotonic()
    for lo in range(0, m, 384):
        hi = min(m, lo + 384)
        common = mask[lo:hi] @ mask.T
        d2 = (x2[lo:hi] @ mask.T + mask[lo:hi] @ x2.T - 2 * x[lo:hi] @ x.T)
        d2 *= p / common
        d2 = np.maximum(d2, 0)
        d2[np.arange(hi-lo), np.arange(lo, hi)] = 0
        pick = np.argpartition(d2, k-1, axis=1)[:, :k]
        ds = np.take_along_axis(d2, pick, axis=1)
        order = np.argsort(ds, axis=1, kind='stable')
        indices[lo:hi] = np.take_along_axis(pick, order, axis=1)
        distances[lo:hi] = np.sqrt(np.take_along_axis(ds, order, axis=1))
        if lo % (384*8) == 0 or hi == m:
            print(f'Exact profile neighbours: {hi:,}/{m:,} ({time.monotonic()-t:.1f}s)', flush=True)
    assert (indices[:, 0] == np.arange(m)).all()
    assert (distances[:, 0] == 0).all()
    assert (np.diff(distances, axis=1) >= 0).all()
    errors = []
    # Check distances independently, including candidate-candidate pairs.
    ids = [z['genes'].tolist().index(g) for g in ['DDX3X', 'KDM5C', 'TADA2B', 'TAF6L', 'USP7', 'TP53']]
    ids += list(np.random.default_rng(42).choice(m, 40, replace=False))
    for i in ids:
        for j, stored in zip(indices[i, 1:11], distances[i, 1:11]):
            ok = np.isfinite(profiles[i]) & np.isfinite(profiles[j])
            direct = np.sqrt(p / ok.sum() * ((profiles[i, ok] - profiles[j, ok]) ** 2).sum())
            errors.append(abs(float(stored)-direct))
    assert max(errors) < 1e-6
    np.savez_compressed(DATA / 'exact_knn.npz', indices=indices, distances=distances)
    rows = []
    selected = sorted(set(z['features']) | set(['KDM5C','KDM6A','KDM5D','UTY','EIF1AX','RPS4X','ZFX','KIF1A']))
    for g in selected:
        if g not in z['genes']: continue
        i = z['genes'].tolist().index(g)
        for rank, (j, d) in enumerate(zip(indices[i, 1:21], distances[i, 1:21]), 1):
            rows.append(dict(gene=g, rank=rank, neighbour=z['genes'][j], profile_distance=float(d)))
    pd.DataFrame(rows).to_csv(DATA / 'selected_profile_neighbours.csv', index=False)
    (DATA / 'distance_validation.json').write_text(json.dumps(dict(
        independent_pairs=len(errors), max_absolute_error=max(errors),
        exact_neighbours=True, k_including_self=k, time_seconds=time.monotonic()-t), indent=2))


def embed():
    from umap import UMAP
    from scipy.sparse import coo_matrix, save_npz
    from scipy.sparse.csgraph import connected_components
    z = np.load(DATA / 'input_profiles.npz')
    knn = np.load(DATA / 'exact_knn.npz')
    indices, distances = knn['indices'], knn['distances']
    # Complete raw input only satisfies UMAP's validation. Geometry is exclusively
    # the supplied exact kNN graph; self r=1 values never enter those distances.
    raw = z['raw_r'].astype(np.float32)
    edge = coo_matrix((np.ones(len(raw)*29),
                      (np.repeat(np.arange(len(raw)),29),indices[:,1:30].ravel())),
                     shape=(len(raw),len(raw))).tocsr()
    n_components, _ = connected_components(edge, directed=False)
    assert n_components == 1, 'Review spectral component placement if graph disconnects'
    fits = [('primary',30,.3,42), ('neighbours15',15,.3,42),
            ('neighbours50',50,.3,42), ('seed43',30,.3,43)]
    parameters = []
    for name, neighbours, min_dist, seed in fits:
        parameters.append(dict(name=name, n_neighbors=neighbours, min_dist=min_dist, seed=seed,
                               n_epochs=500, init='spectral', spread=1,
                               graph_connected_components=n_components))
        path = DATA / f'embedding_{name}.csv'
        if path.exists():
            print(f'Using saved {name} embedding', flush=True)
            continue
        print(f'Fitting {name}: n_neighbors={neighbours}, min_dist={min_dist}, seed={seed}', flush=True)
        reducer = UMAP(n_neighbors=neighbours, min_dist=min_dist, spread=1,
                       n_components=2, metric='euclidean', init='spectral',
                       n_epochs=500, random_state=seed, transform_seed=seed,
                       n_jobs=1, precomputed_knn=(indices.copy(), distances.copy(), None),
                       force_approximation_algorithm=True, verbose=False)
        result = reducer.fit_transform(raw)
        assert result.shape == (len(raw),2) and np.isfinite(result).all()
        # Check the graph really used the frozen self-excluded neighbour distances.
        assert np.array_equal(reducer._knn_indices, indices[:,:neighbours])
        assert np.allclose(reducer._knn_dists, distances[:,:neighbours], rtol=0, atol=0)
        pd.DataFrame(result,index=z['genes'],columns=['UMAP1','UMAP2']).to_csv(path,index_label='gene')
        if name == 'primary': save_npz(DATA / 'primary_fuzzy_graph.npz',reducer.graph_)
        print(f'Saved {name}', flush=True)
    report=dict(parameters=parameters, geometry='Exact masked Euclidean profile distances; see analysis_manifest',
                versions={p:importlib.metadata.version(p) for p in ['numpy','scipy','scikit-learn','umap-learn','numba','pynndescent']},
                graph_validation='Supplied neighbour indices and distances exactly match the fitted UMAP graph inputs',
                no_formal_clusters='No discrete cluster labels or significance inferred from 2D geometry')
    (DATA / 'embedding_manifest.json').write_text(json.dumps(report,indent=2),encoding='utf-8')


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['prepare','knn','embed','all'],default='all')
    args=parser.parse_args()
    if args.stage in ['prepare','all']: prepare()
    if args.stage in ['knn','all']: exact_knn()
    if args.stage in ['embed','all']: embed()
