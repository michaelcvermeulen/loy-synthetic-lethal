"""Scientific and delivery checks for the batch-aware gene-level extension."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
ROOT=Path(__file__).resolve().parents[1];T=ROOT/"outputs/tables"
inputs=pd.read_csv(ROOT/"provenance/batchDE_input_hashes.csv")
for row in inputs.itertuples(): assert hashlib.sha256((ROOT/row.path).read_bytes()).hexdigest()==row.sha256,row.path
ids=pd.read_csv(T/"HT1080_LOY_vs_sgCTRL_primary_effects.csv",usecols=["gene_id"]).gene_id.tolist()
assert len(ids)==25526 and len(set(ids))==25526
frames={}
for model,contrasts,df in [("passage_period",["combined","AUG_2025","MAY_2026","MAY_minus_AUG"],4),("batch_condition",["combined"],9)]:
    D=pd.read_csv(T/f"batchDE_design_{model}.csv").drop(columns="sample_key").to_numpy()
    assert D.shape[0]==12 and np.linalg.matrix_rank(D)==D.shape[1] and 12-D.shape[1]==df
    for contrast in contrasts:
        z=pd.read_csv(T/f"batchDE_{model}__{contrast}.csv");frames[(model,contrast)]=z
        assert z.gene_id.tolist()==ids
        assert np.isfinite(z[["logFC","PValue","FDR","CI.L","CI.R"]].to_numpy()).all()
        assert ((z.PValue>=0)&(z.PValue<=1)&(z.FDR>=0)&(z.FDR<=1)).all()
        assert ((z['CI.L']<=z.logFC)&(z.logFC<=z['CI.R'])).all()
        order=np.argsort(z.PValue.to_numpy());p=z.PValue.to_numpy()[order]
        bh=np.minimum(1,np.minimum.accumulate((p*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
        assert np.allclose(z.FDR.to_numpy()[order],bh,atol=1e-12,rtol=1e-10)
aug=frames[("passage_period","AUG_2025")].logFC.to_numpy();may=frames[("passage_period","MAY_2026")].logFC.to_numpy()
assert np.allclose(frames[("passage_period","combined")].logFC,.5*(aug+may),atol=1e-10)
assert np.allclose(frames[("passage_period","MAY_minus_AUG")].logFC,may-aug,atol=1e-10)
for report in ["analysis.html","batch_aware_de.html"]:
    html=ROOT/"outputs/report"/report;assert html.exists()
    soup=BeautifulSoup(html.read_text(encoding="utf-8"),"html.parser");text=soup.get_text(" ",strip=True)
    assert "Batch-aware gene-level differential expression" in text
    assert "4674" in text and "1003" in text
    assert "Relation to the network modules" in text
    assert any(str(img.get("src","")).startswith("data:image") for img in soup.find_all("img"))
manifest=pd.read_csv(ROOT/"outputs/figures/batchDE_figure_manifest.csv")
for row in manifest.itertuples(): assert hashlib.sha256((ROOT/"outputs/figures"/row.file).read_bytes()).hexdigest()==row.sha256
print("PASS: frozen original inputs, exact gene universe, full-rank designs, all finite results, independently checked BH, equal-period/interaction contrasts, rendered reports and figure hashes.")
