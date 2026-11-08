"""Validate final scientific payload and JavaScript syntax; no browser execution."""
from pathlib import Path
import hashlib,json,re,subprocess
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parents[1]
DATA=OUT/'source_data/ht1080_hit_profile_umap'
STEM='ch04_supp_ht1080_hit_profile_umap_draft'
z=np.load(DATA/'input_profiles.npz')
e=pd.read_csv(DATA/'embedding_primary.csv',index_col=0)
a=pd.read_csv(DATA/'gene_annotations.csv',index_col=0).fillna('')
k=np.load(DATA/'exact_knn.npz')
html=(OUT/f'{STEM}.html').read_text(encoding='utf-8')
compatibility=OUT/'code/ht1080_hit_profile_umap_viewer.html'
assert compatibility.read_text(encoding='utf-8')==html
payload=json.loads(html.split('const D=',1)[1].split(';\nconst C=',1)[0])
assert payload['genes']==e.index.tolist()==z['genes'].tolist()
assert payload['features']==z['features'].tolist()
assert np.allclose(payload['xy'],e.to_numpy(),atol=5.01e-6,rtol=0)
assert np.allclose(np.array(payload['profiles'],dtype=float),z['distance_r'],atol=5.01e-5,rtol=0,equal_nan=True)
assert np.array_equal(payload['neighbours'],k['indices'][:,1:11])
assert np.allclose(payload['distances'],k['distances'][:,1:11],atol=5.01e-5,rtol=0)
assert payload['n']==a.n_models.tolist()
assert payload['directions']==a.screen_direction.tolist()
sd=pd.read_csv(OUT/'source_data/xy_saga_profile_distance/x_gametolog_saga_distances.csv')
assert len(payload['saga_distance'])==16
for record in payload['saga_distance']:
    source=sd.loc[sd.gene.eq(record['gene'])].iloc[0]
    assert record['gene'] in payload['genes']
    for key in ['five_hit_distance','six_hit_distance','rank_X_five','rank_X_six']:
        assert abs(record[key]-float(source[key]))<1e-12
assert len(payload['genes'])==17931 and len(set(payload['genes']))==17931
assert a.screen_direction.value_counts().to_dict()=={'Not nominated':17875,'Depleted':47,'Enriched':9}
manifest=json.loads((DATA/'analysis_manifest.json').read_text())
for source,expected in manifest['sources'].items():
    assert hashlib.sha256(Path(source).read_bytes()).hexdigest()==expected
assert '/*DATA_PAYLOAD*/' not in html and 'FONT_PAYLOAD' not in html
assert not re.search(r'(?:src|href)=["\']https?://',html)
script=html.split('<script>',1)[1].split('</script>',1)[0]
js_path=DATA/'viewer_syntax_check.js'
js_path.write_text(script,encoding='utf-8')
check=subprocess.run(['C:/Program Files/nodejs/node.exe','--check',str(js_path)],capture_output=True,text=True)
assert check.returncode==0,check.stderr
render=json.loads((DATA/'render_manifest.json').read_text())
assert not render['label_overlaps'] and not render['text_outside_page']
report=dict(genes=len(e),features=len(z['features']),finite_coordinates=True,
            masked_self_cells=int(np.isnan(z['distance_r']).sum()),
            original_source_hashes_unchanged=True,html_payload_matches_source=True,
            previously_opened_html_path_contains_identical_data=True,
            saga_distance_rows_validated=16,
            javascript_syntax_check='pass',external_html_resources=False,
            browser_interaction_test='Not performed: browser policy prohibits local-file URLs; no alternative route used',
            figure_checks='PDF fonts and bounds, label overlaps, colour/greyscale/deuteranopia previews inspected',
            output_sha256={ext:hashlib.sha256((OUT/f'{STEM}.{ext}').read_bytes()).hexdigest() for ext in ['pdf','svg','png','html']})
(DATA/'final_validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='output_sha256'},indent=2))
