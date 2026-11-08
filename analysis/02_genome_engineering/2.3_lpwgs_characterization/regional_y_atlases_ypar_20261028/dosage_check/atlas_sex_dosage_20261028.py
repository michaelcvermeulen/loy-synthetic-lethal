"""Fixed-bin sex-chromosome dosage check for atlas caption interpretation.

Read depth is not a cytogenetic karyotype. Predictions are conditional on the
non-PAR dosage representing PAR dosage and complete loss of Y-derived PAR.
Run in the pilot WSL environment; existing alignments and fits are read only.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import csv, json, statistics, subprocess, hashlib

ROOT = Path(__file__).resolve().parent
THESIS = ROOT.parents[2] / 'SEPT_11_2026_THESIS'
OUT = ROOT / 'atlas_sex_dosage_20261028'
OUT.mkdir(exist_ok=True)
SAMTOOLS = ROOT / 'runtime/env/bin/samtools'
BED = ROOT / 'a549_par_dosage_20261028/fixed_strict_nonPAR_XY_bins.bed'
PAIRS = {'HT1080': ('HT1080-sgCTRL-D0', 'HT1080-LOY-D0'),
         'DU145': ('LNCaP-sgCTRL', 'LNCaP-LOY-2X')}

def read(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def save(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

def process(cell):
    samples = PAIRS[cell]
    fit_path = THESIS / f'research/model_development/ht1080_du145_ypar_20261027/{cell}/masked/bins.csv'
    fit = {(r['sample'], r['bin_id']): float(r['copynumber']) for r in read(fit_path)
           if r['copynumber'] and r['use'].upper() == 'TRUE'}
    norm = {r['sample']: float(r['median_autosomal_depth']) for r in
            read(ROOT / f'atlas_pair_coverage/{cell}/masked/autosomal_normalization.csv')}
    data, summary = [], []
    for sample in samples:
        dest = OUT / f'{sample}.nonPAR_XY.bedcov.tsv'
        with dest.open('w') as f:
            subprocess.run([str(SAMTOOLS), 'bedcov', '-j', '-Q', '37', '-G', '0x800', str(BED),
                            str(ROOT / f'alignment/{sample}.YPARmasked.bam')], stdout=f, check=True)
        for line in dest.read_text().splitlines():
            chrom, start, end, name, bases = line.split('\t')
            region, bin_id = name.split(':', 1)
            depth = int(bases) / (int(end) - int(start))
            data.append(dict(cell=cell, sample=sample, region=region, bin_id=bin_id,
                             relative_primary_depth=depth/norm[sample], corrected_CN=fit[sample, bin_id]))
        for region in ['X_nonPAR', 'MSY_Yp', 'euchromatic_Yq']:
            rows = [r for r in data if r['sample'] == sample and r['region'] == region]
            summary.append(dict(cell=cell, sample=sample, region=region, n_fixed_bins=len(rows),
                                median_relative_primary_depth=statistics.median(r['relative_primary_depth'] for r in rows),
                                median_corrected_CN_same_bins=statistics.median(r['corrected_CN'] for r in rows)))
    return data, summary, str(fit_path), hashlib.sha256(fit_path.read_bytes()).hexdigest()

with ThreadPoolExecutor(max_workers=2) as pool:
    results = list(pool.map(process, PAIRS))
data = [r for d, s, p, h in results for r in d]
summary = [r for d, s, p, h in results for r in s]
for r in read(ROOT / 'a549_par_dosage_20261028/nonPAR_XY_dosage.csv'):
    summary.append(dict(cell='A549', sample=r['sample'], region=r['region'], n_fixed_bins=int(r['n_fixed_bins']),
                        median_relative_primary_depth=float(r['median_relative_primary_depth']),
                        median_corrected_CN_same_bins=float(r['median_corrected_CN_same_bins'])))
save(OUT / 'nonPAR_XY_dosage.csv', summary)
save(OUT / 'HT1080_DU145_per_bin_dosage.csv', data)
predictions = []
for cell, (control, loy) in dict(PAIRS, A549=('A549-CTRL2', 'A549-LOY')).items():
    lookup = {(r['sample'], r['region']): r for r in summary if r['cell'] == cell}
    for y_region in ['MSY_Yp', 'euchromatic_Yq']:
        for method in ['median_relative_primary_depth', 'median_corrected_CN_same_bins']:
            xc = lookup[control, 'X_nonPAR'][method]
            yc = lookup[control, y_region][method]
            xl = lookup[loy, 'X_nonPAR'][method]
            predictions.append(dict(cell=cell, Y_region=y_region, method=method,
                                    control_X_to_Y=xc/yc, LOY_to_control_X=xl/xc,
                                    expected_PAR_percent_fixed_X=100*xc/(xc+yc),
                                    expected_PAR_percent_observed_X=100*xl/(xc+yc)))
save(OUT / 'conditional_PAR_predictions.csv', predictions)
(OUT / 'provenance.json').write_text(json.dumps(dict(
    fixed_bins='Original strict reference mask, wholly within non-PAR X, MSY Yp or euchromatic Yq; no sample-dependent selection',
    BED=str(BED), BED_sha256=hashlib.sha256(BED.read_bytes()).hexdigest(),
    read_filter='samtools bedcov -j -Q37 -G0x800: primary mapped, no QC-failed/marked duplicate records',
    normalizer='The same 21,446 strict autosomal bins as in the PAR plots',
    fits={p:h for d,s,p,h in results},
    limitations='Bulk relative dosage does not establish absolute chromosome counts or structural karyotype. X-adjusted predictions assume comparable autosomal scaling and complete loss of the Y PAR contribution; they are not estimates of actual Y-PAR loss.'
), indent=2))
print(json.dumps(predictions, indent=2))
