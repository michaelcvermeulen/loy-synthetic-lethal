"""Convert pinned JASPAR PFMs to HOMER known motifs with explicit score thresholds."""
from pathlib import Path
import numpy as np
import MOODS.tools
ROOT=Path(__file__).resolve().parents[1]
cfg=ROOT/'config/TF_extension'
lines=(cfg/'JASPAR2024_CORE_vertebrates.jaspar').read_text().splitlines()
output=[]; accounting=[]
for i in range(0,len(lines),5):
    header=lines[i][1:].split(); mid=header[0];name=header[1]
    pfm=np.array([[float(x) for x in row[row.index('[')+1:row.index(']')].split()] for row in lines[i+1:i+5]])
    # Positive probabilities avoid undefined log odds. Use the same probabilities for thresholding and HOMER.
    probs=(pfm+.25)/(pfm.sum(axis=0,keepdims=True)+1)
    probs=np.maximum(probs,.001);probs=probs/probs.sum(axis=0,keepdims=True)
    pwm=MOODS.tools.log_odds(probs.tolist(),[.25]*4,0)
    threshold=MOODS.tools.threshold_from_p(pwm,[.25]*4,5e-5)
    consensus=''.join('ACGT'[k] for k in probs.argmax(axis=0))
    output.append(f'>{consensus}\t{mid}_{name}\t{threshold:.10f}')
    output.extend('\t'.join(f'{x:.10f}' for x in col) for col in probs.T)
    accounting.append(f'{mid}\t{name}\t{probs.shape[1]}\t{threshold:.10f}')
(cfg/'JASPAR2024_CORE_vertebrates.homer').write_text('\n'.join(output)+'\n')
(cfg/'HOMER_motif_thresholds.tsv').write_text('motif_id\tmotif_name\twidth\tln_odds_threshold\n'+'\n'.join(accounting)+'\n')
print('Converted',len(accounting),'motifs at uniform-background site P=5e-5')
