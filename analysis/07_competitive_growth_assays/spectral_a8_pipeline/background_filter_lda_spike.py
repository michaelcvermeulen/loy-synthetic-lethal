from pathlib import Path
import numpy as np
from spectral.a8_pipeline import competition_report as c

repo = Path(__file__).parents[2]
d4 = repo / "spectral/20260901_125230_Mike_Run1/Mike_Run1"
g = c.load_cef_gates(d4 / "20260901_125230_Mike_Run1.cef")
files = c._file_map(d4 / "1")
roles = {"A03-CTRL_GFP": 1, "A05-LOY_GFP": 1, "A07-CTRL_cherry": 0, "A09-LOY_cherry": 0}
train_x=[]; train_y=[]
for n,label in roles.items():
    x=c._primary_features(files[n],g)[0]; train_x.append(x); train_y.append(np.full(len(x),label))
lda=c.fit_lda(np.vstack(train_x),np.concatenate(train_y))
u=np.vstack([c._primary_features(files[n],g)[0] for n in ("A01-unstained_CTRL","A02-unstained_LOY")])
th=c.background_gate_thresholds(u,95)
print("95th-percentile thresholds",th.tolist())
print("sample | retained | LDA GFP-like among retained | LDA mCherry-like among retained")
for n in sorted(files):
    if not n.startswith(("A","B","C","D","E")): continue
    x=c._primary_features(files[n],g)[0]; keep=np.any(x>th,axis=1); pred=lda.predict(x[keep]); total=len(x)
    print(n, f"{100*keep.mean():.1f}%", f"{100*pred.mean():.1f}%", f"{100*(1-pred.mean()):.1f}%", sep=" | ")
