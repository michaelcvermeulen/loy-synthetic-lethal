"""Reproduce the screen notebook's scaled log2(CPM+1) sample PCA."""
from pathlib import Path
import hashlib
import json
import sys

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
import numpy as np
import pandas as pd

DATA = THESIS / 'figures/chapter4/source_data/ht1080_figure1_qc'
SOURCE = DATA / 'guide_log2_cpm.csv'
NOTEBOOK = THESIS.parent / 'crispr_screen/ht1080_loy_screen_2/code/ht1080_screen2_analysis.Rmd'


def main():
    log = pd.read_csv(SOURCE)
    samples = list(log.columns[2:])
    x = log[samples].to_numpy().T
    assert np.isfinite(x).all()
    sd = x.std(axis=0, ddof=1)
    keep = sd > 0
    z = (x[:, keep] - x[:, keep].mean(axis=0)) / sd[keep]
    u, singular, vt = np.linalg.svd(z, full_matrices=False)
    # Deterministic signs: the largest-magnitude sample score is positive.
    scores = u * singular
    for j in range(scores.shape[1]):
        if scores[np.argmax(np.abs(scores[:, j])), j] < 0:
            scores[:, j] *= -1
            vt[j] *= -1
    pct = 100 * singular**2 / np.sum(singular**2)
    out = pd.DataFrame(scores[:, :5], index=samples, columns=[f'PC{i}' for i in range(1, 6)])
    out.index.name = 'sample'
    out['genotype'] = ['LOY' if 'LOY' in s else 'CTRL' for s in samples]
    out['stage'] = ['D0' if s.startswith('D0') else 'Endpoint' for s in samples]
    out.to_csv(DATA / 'pca_sample_scores.csv')
    pd.DataFrame({'component': range(1, 7), 'singular_value': singular, 'variance_percent': pct}).to_csv(DATA / 'pca_variance.csv', index=False)
    features = log.loc[keep, ['sgRNA', 'Gene']].copy()
    features['mean_log2_cpm_plus_1'] = x[:, keep].mean(axis=0)
    features['sd_log2_cpm_plus_1'] = sd[keep]
    for j in range(5):
        features[f'PC{j+1}_loading'] = vt[j]
    features.to_csv(DATA / 'pca_guide_loadings.csv.gz', index=False)
    # Independent eigendecomposition of the sample Gram matrix checks the SVD.
    eigen = np.linalg.eigvalsh(z @ z.T)[::-1]
    assert np.allclose(singular[:5]**2, eigen[:5], rtol=1e-10)
    assert np.allclose(scores.sum(axis=0), 0, atol=1e-8)
    manifest = {
        'method': 'PCA of samples by guide log2(CPM+1); each guide centered and scaled to sample SD=1 (ddof=1), matching prcomp(scale.=TRUE) in the existing screen notebook.',
        'samples': samples, 'guides_input': len(log), 'guides_retained': int(keep.sum()),
        'zero_variance_guides_excluded': int((~keep).sum()),
        'variance_percent': pct.tolist(), 'batch_adjustment': 'none',
        'sign_convention': 'Largest absolute sample score positive per component; axis signs are arbitrary.',
        'independent_gram_eigenvalue_check': True,
        'sources': {p.name: {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in [SOURCE, NOTEBOOK]},
    }
    (DATA / 'pca_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in manifest.items() if k != 'sources'}, indent=2))


if __name__ == '__main__':
    main()
