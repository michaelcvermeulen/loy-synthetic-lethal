"""Recover target-specific DU145 drug comparisons from the March 2026 export."""
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'figures/chapter3/source_data/du145_drug_optimisation'
SOURCE = Path('A:/Dropbox/LOY_paralogs/data/ddPCR/AMELX Craig_MARCH_25_2026_2026-03-25-22-38/data_table.csv')
MAPPING = {
    'LNCaP_LOY_DMSO': ('DMSO vehicle', 0, False),
    'LNCaP_LOY_KU': ('KU-55933', 0, True),
    'LNCaP_LOY_Pano_30nM': ('Panobinostat 30 nM', 30, False),
    'LNCaP_LOY_Pano_70nM': ('Panobinostat 75 nM', 75, False),
    'LNCaP_LOY_KU_Pano_30nM': ('KU-55933 + panobinostat 30 nM', 30, True),
    'LNCaP_LOY_KU_Pano_70nM': ('KU-55933 + panobinostat 75 nM', 75, True),
    'LNCaP_sgCTRL': ('sgCTRL', None, None),
    'LNCaP_LOY_pass1': ('LOY pass1 (source label)', None, None),
}

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    frozen = OUT / '2026-03-25_data_table.csv'
    if not frozen.exists():
        shutil.copy2(SOURCE, frozen)
    d = pd.read_csv(frozen)
    concentration = next(c for c in d if c.startswith('Conc(copies/'))
    rows = []
    for _, y in d[d.Target.isin(['AMELY', 'KDM5D', 'DDX3Y'])].iterrows():
        source_name = y['Sample description 1']
        condition, pano, ku = MAPPING[source_name]
        ref = 'AMELX' if y.Target == 'AMELY' else 'RNaseP'
        x = d[d.Well.eq(y.Well) & d.Target.str.upper().eq(ref.upper())]
        assert len(x) == 1
        x = x.iloc[0]
        yc, xc = float(y[concentration]), float(x[concentration])
        assert xc > 0
        ratio = yc/xc
        assert np.isclose(ratio, float(y.Ratio))
        rows.append(dict(source_sample=source_name, cell_line='DU145', condition=condition,
            panobinostat_nM=pano, ku_present=ku, well=y.Well,
            target=y.Target, reference=ref, target_concentration=yc,
            reference_concentration=xc, concentration_ratio=ratio,
            instrument_ratio=float(y.Ratio), target_status=y.Status,
            reference_status=x.Status, assay_date='2026-03-25'))
    result = pd.DataFrame(rows)
    assert len(result) == 24 and not result.duplicated(['condition','target']).any()
    vehicle = result[result.condition.eq('DMSO vehicle')].set_index('target').concentration_ratio
    primary = result[~result.condition.isin(['sgCTRL','LOY pass1 (source label)'])].copy()
    primary['vehicle_ratio'] = primary.target.map(vehicle)
    primary['percent_reduction_vs_vehicle'] = 100*(1-primary.concentration_ratio/primary.vehicle_ratio)
    # Preserve negative reductions; never pool distinct genes or historical controls.
    result.to_csv(OUT/'all_target_measurements.csv', index=False)
    primary.to_csv(OUT/'drug_vs_vehicle.csv', index=False)
    provenance = dict(source=str(SOURCE), frozen_sha256=hashlib.sha256(frozen.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        target_rows=24, treatment_target_rows=18, observations_per_condition_target=1,
        author_corrections=['DU145, historical LNCaP names retained',
            'Higher panobinostat dose 75 nM; raw 70 nM labels unchanged',
            'KU-55933 before nucleofection; panobinostat afterwards'],
        unresolved=['KU-55933 dose', 'Exposure and recovery durations',
            'Treatment-population derivation and sequencing links'],
        analysis='Concentration ratios; matched vehicle by target; no absolute-copy scaling, historical pooling, clipping or inferential tests')
    (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    print(primary.pivot(index='condition',columns='target',values='percent_reduction_vs_vehicle').round(1).to_string())

if __name__ == '__main__':
    main()
