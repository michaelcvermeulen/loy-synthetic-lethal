"""Descriptive 20-aa window scan of the saved canonical X/Y alignment."""
from pathlib import Path
import csv, json, sys

THESIS = Path(__file__).resolve().parents[3]
PROJECT = THESIS.parent
BASE = THESIS / 'figures/chapter4/source_data/ddx3x_ddx3y_inhibitor_site_conservation'
OUT = BASE / 'divergent_regions'
OUT.mkdir(exist_ok=True)
sys.path[:0] = [str(PROJECT / '.tmp/ddx3_bio'), str(PROJECT / '.figure_python_deps')]
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

with (BASE / 'data/alignment_residue_map.tsv').open(encoding='utf-8') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
for r in rows:
    r['x_pos'] = int(r['x_pos']) if r['x_pos'] else None
    r['y_pos'] = int(r['y_pos']) if r['y_pos'] else None
xmap = {r['x_pos']: r for r in rows if r['x_pos']}

def window(start, length=20):
    rr = [xmap[p] for p in range(start, start + length)]
    n = len(rr)
    identical = sum(r['status'] == 'identical' for r in rr)
    positive = sum(r['status'] in ['identical', 'conservative'] for r in rr)
    end = start + length - 1
    return dict(x_start=start, x_end=end,
                y_start=min(r['y_pos'] for r in rr if r['y_pos']),
                y_end=max(r['y_pos'] for r in rr if r['y_pos']),
                columns=n, identical=identical, positive=positive,
                differences=n-identical, gap_columns=sum(r['status']=='gap' for r in rr),
                identity_pct=100*identical/n, similarity_pct=100*positive/n,
                within_uniprot_helicase_domain=(211<=start and end<=403) or (414<=start and end<=575),
                x_sequence=''.join(r['x_aa'] for r in rr), y_sequence=''.join(r['y_aa'] for r in rr))

windows = [window(s) for s in range(1, 662-20+2)]
global_best = max(w['differences'] for w in windows)
domain_best = max(w['differences'] for w in windows if w['within_uniprot_helicase_domain'])
global_ties = [w for w in windows if w['differences']==global_best]
domain_ties = [w for w in windows if w['within_uniprot_helicase_domain'] and w['differences']==domain_best]
selected = [dict(region='N-terminal IDR', selection='highest mismatch count across all 20-aa windows', **global_ties[0]),
            dict(region='RecA2 domain', selection='highest mismatch count within either UniProt helicase domain; smallest start among ties', **domain_ties[0])]
assert [(w['x_start'],w['x_end'],w['identical']) for w in selected] == [(9,28,10),(421,440,15)]
assert all(w['gap_columns']==0 for w in selected)

def write(name, values):
    with (OUT / name).open('w', encoding='utf-8', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(values[0]), delimiter='\t')
        w.writeheader(); w.writerows(values)

write('all_20aa_windows.tsv', windows)
write('selected_windows.tsv', selected)
write('selected_residue_map.tsv', [dict(region=w['region'], **xmap[p]) for w in selected for p in range(w['x_start'],w['x_end']+1)])

# Check whether the variable folded segment is present in an experimental structure.
canonical = json.loads((BASE/'sources/O00571.json').read_text())['sequence']['value']
chain=PDBParser(QUIET=True).get_structure('2I4I',BASE/'sources/2I4I.pdb')[0]['A']
residues=[r for r in chain if r.id[0]==' ' and 'CA' in r]
sequence=''.join(seq1(r.resname) for r in residues)
aligner=PairwiseAligner(mode='global',substitution_matrix=substitution_matrices.load('BLOSUM62'),open_gap_score=-10,extend_gap_score=-.5)
aligner.end_gap_score=0
alignment=aligner.align(canonical,sequence)[0]
xpos=ppos=0; mapping=[]
for a,b in zip(str(alignment[0]),str(alignment[1])):
    xpos+=a!='-'; ppos+=b!='-'
    if a!='-' and b!='-' and 421<=xpos<=440:
        assert a==b
        mapping.append(dict(x_pos=xpos,x_aa=a,pdb='2I4I',chain='A',pdb_residue=residues[ppos-1].id[1]))
assert len(mapping)==20
write('folded_window_pdb_mapping.tsv',mapping)
audit=dict(window_length=20,step=1,total_windows=len(windows),
           denominator='20 alignment columns anchored to consecutive DDX3X residues; no X gaps in saved alignment',
           difference_definition='all nonidentical columns, including gaps',
           domain_intervals=[[211,403],[414,575]],
           global_best_ties=[[w['x_start'],w['x_end']] for w in global_ties],
           domain_best_ties=[[w['x_start'],w['x_end']] for w in domain_ties],
           tie_rule='smallest DDX3X start coordinate',selected=selected,
           folded_window_resolved_in_2I4I=True,
           limitation='Descriptive sequence-based prioritisation. No cavity, accessibility, ligand binding or PROTAC ternary-complex suitability was established.',
           disorder_source='Owens et al. 2024, doi:10.1038/s41467-024-53636-0; N-terminal IDR 1–131')
(OUT/'selection_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps(audit,indent=2))
