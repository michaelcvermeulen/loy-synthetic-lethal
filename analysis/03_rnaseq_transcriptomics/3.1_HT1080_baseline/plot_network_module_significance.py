"""Plot the existing top-annotation P and BH q values for every RNA network module.

This is a presentation-only extension. No network, annotation or RNA test is rerun.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent
TABLES = ROOT / "outputs/tables/network_modules"
FIGURES = ROOT / "outputs/figures/network_modules"
STEM = "HT1080_LOY_all_module_annotation_significance"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def probability_text(value):
    return f"{value:.2e}" if value < .001 else f"{value:.3g}"


def panel(fig, data, rect, color, letter, title):
    ax = fig.add_axes(rect)
    y = np.arange(len(data))
    ax.set_ylim(len(data)-.5, -.5)
    ax.set_xlim(-.6, 27)
    ax.set_xticks([0, 5, 10, 15, 20, 25])
    ax.tick_params(axis="x", labelsize=6, length=2.5, width=.5, color="#58636B")
    ax.set_yticks([])
    for side in ["left", "right", "top"]:
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(.5)
    ax.spines["bottom"].set_color("#58636B")
    ax.grid(axis="x", color="#E8ECEF", linewidth=.45)
    ax.set_axisbelow(True)
    ax.axvline(-np.log10(.05), color="#747F85", linewidth=.6,
               linestyle=(0, (3, 2)), zorder=1)
    trans = blended_transform_factory(fig.transFigure, ax.transData)
    for i, row in enumerate(data.itertuples()):
        if i % 2 == 0:
            ax.axhspan(i-.5, i+.5, color="#F5F7F8", zorder=0)
        ax.plot([row.neglog10_q, row.neglog10_P], [i, i], color=color,
                alpha=.6, linewidth=.7, zorder=2)
        ax.scatter(row.neglog10_P, i, s=16, facecolors="white", edgecolors=color,
                   linewidths=.8, marker="o", zorder=4)
        ax.scatter(row.neglog10_q, i, s=13, facecolors=color, edgecolors=color,
                   linewidths=.4, marker="D", zorder=5)
        label = f"{row.display_id}  {row.module_label}"
        ax.text(.03, i, label, transform=trans, fontsize=6,
                va="center", ha="left", color="#202B32", clip_on=False)
        ax.text(.410, i, str(row.module_genes), transform=trans, fontsize=6,
                va="center", ha="right", color="#566169", clip_on=False)
        ax.text(.866, i, probability_text(row.P_value), transform=trans,
                fontsize=6, va="center", ha="right", clip_on=False)
        ax.text(.978, i, probability_text(row.BH_all_modules_and_terms), transform=trans,
                fontsize=6, va="center", ha="right", clip_on=False)
    top = rect[1]+rect[3]
    fig.text(.011, top+.028, letter, fontsize=8, fontweight="bold", color="#202B32")
    fig.text(.03, top+.028, title, fontsize=7, fontweight="bold", color=color)
    for x, label, align in [(.03, "Module", "left"), (.410, "Genes", "right"),
                            (.866, "P value", "right"), (.978, "BH q", "right")]:
        fig.text(x, top+.009, label, fontsize=6, fontweight="bold", ha=align, color="#39464E")
    return ax


def main():
    inputs = [TABLES / "module_top_annotations.csv", TABLES / "modules.csv",
              ROOT / "config/network_module_labels.json",
              ROOT / "provenance/network_modules/analysis_summary.json"]
    before = {p.relative_to(ROOT).as_posix(): sha(p) for p in inputs}
    source = pd.read_csv(inputs[0], float_precision="round_trip")
    modules = pd.read_csv(inputs[1], float_precision="round_trip")
    labels = json.loads(inputs[2].read_text(encoding="utf-8"))
    summary = json.loads(inputs[3].read_text())
    data = source.merge(modules[["module_id", "rank", "nodes", "mcode_score", "physical_edges", "edges"]],
                        on="module_id", validate="one_to_one")
    assert len(data) == 42 and set(data.module_id) == set(modules.module_id)
    assert data.groupby("direction").size().to_dict() == {"DOWN":24, "UP":18}
    assert data.background.eq("expressed_STRING_mapped").all()
    assert data.background_genes.eq(14564).all() and data.module_genes.eq(data.nodes).all()
    assert data.overlap_genes.ge(2).all()
    assert data.P_value.between(0,1,inclusive="right").all()
    assert data.BH_all_modules_and_terms.between(0,1,inclusive="right").all()
    assert (data.BH_all_modules_and_terms >= data.P_value).all()
    assert summary["annotation_tests_per_background"] == 260526
    assert summary["annotation_terms"] == 6203
    data["module_label"] = data.module_id.map(labels)
    assert data.module_label.notna().all()
    data["display_id"] = data.module_id.str.replace("DOWN_", "D", regex=False).str.replace("UP_", "U", regex=False)
    data["neglog10_P"] = -np.log10(data.P_value)
    data["neglog10_q"] = -np.log10(data.BH_all_modules_and_terms)
    data["annotation_q_below_005"] = data.BH_all_modules_and_terms < .05
    data = data.sort_values(["direction", "rank"]).reset_index(drop=True)
    assert data.neglog10_P.max() < 27
    data_file = TABLES / "annotation_significance_plot_data.csv"
    data.to_csv(data_file, index=False, float_format="%.17g")

    plt.rcParams.update({"font.family":"Arial", "font.size":6,
                         "pdf.fonttype":42, "ps.fonttype":42, "svg.fonttype":"none",
                         "svg.hashsalt":"HT1080_module_significance_42",
                         "mathtext.fontset":"custom", "mathtext.rm":"Arial",
                         "mathtext.it":"Arial:italic", "mathtext.bf":"Arial:bold"})
    fig = plt.figure(figsize=(183/25.4,235/25.4), facecolor="white")
    fig.text(.03,.978,"Functional annotation significance of all network modules",
             fontsize=8.5,fontweight="bold",va="top")
    fig.text(.03,.957,"Existing top annotation for each module; RNA direction defines the two panels",
             fontsize=6.5,color="#53616A",va="top")
    handles = [Line2D([0],[0],marker="o",markersize=4,markerfacecolor="white",markeredgecolor="#293A44",linestyle="none"),
               Line2D([0],[0],marker="D",markersize=3.5,markerfacecolor="#293A44",markeredgecolor="#293A44",linestyle="none"),
               Line2D([0],[0],color="#747F85",linewidth=.6,linestyle=(0,(3,2)))]
    fig.legend(handles,["Raw P value","BH-adjusted q","P or q = 0.05"],loc="upper left",
               bbox_to_anchor=(.024,.939),ncol=3,frameon=False,fontsize=6.3,
               columnspacing=2,handletextpad=.6)
    panel(fig,data[data.direction.eq("DOWN")],[.445,.459,.312,.422],"#0072B2","a","Lower RNA in LOY · 24 modules")
    panel(fig,data[data.direction.eq("UP")],[.445,.080,.312,.3165],"#D55E00","b","Higher RNA in LOY · 18 modules")
    fig.text(.601,.048,r"$-\log_{10}$(P value or BH q)",fontsize=7,ha="center")
    fig.text(.03,.026,"Annotation BH adjustment: 42 modules × 6,203 terms (260,526 tests); expressed STRING-mapped background.",
             fontsize=5.7,color="#53616A")
    fig.text(.03,.012,"These are functional-annotation P values, not tests of module-level RNA change or protein interaction; MCODE has no P value.",
             fontsize=5.7,color="#53616A")
    fig.canvas.draw()
    for ext in ["pdf","svg","png"]:
        fig.savefig(FIGURES / f"{STEM}.{ext}",dpi=600,facecolor="white")
    plt.close(fig)

    caption = ("**Functional-annotation significance for all lower- and higher-RNA network modules in HT1080 LOY.** "
        "a, All 24 lower-RNA modules. b, All 18 higher-RNA modules. Module identifiers and their order match the MCODE atlas; "
        "row labels are descriptive module summaries. Each row uses the previously selected top GO biological-process or Reactome annotation "
        "with at least two overlapping genes. Open circles show raw one-sided hypergeometric P values and filled diamonds show the corresponding "
        "Benjamini–Hochberg-adjusted q values, plotted as −log10 values on a common scale. Connecting lines show the change after adjustment; "
        "the vertical dashed line marks 0.05. Numeric columns report both probabilities. BH correction spans all 260,526 module–term tests "
        "(42 modules × 6,203 terms), using 14,564 expressed STRING-mapped genes as background. The original values were reused without new testing "
        "or a smaller correction family. All modules are retained, including D20 with q = 0.657. Exact tested terms and full-precision probabilities "
        "are supplied in the source CSV. These exploratory annotation statistics describe functional coherence; they do not test module-level RNA "
        "change, protein interaction or causal LOY effects in independently generated populations.")
    caption_file = ROOT / "outputs/captions/network_module_significance.md"
    caption_file.write_text("# All-module annotation significance\n\n"+caption+"\n\n"
        "[Figure PDF](../figures/network_modules/"+STEM+".pdf) · "
        "[Figure PNG](../figures/network_modules/"+STEM+".png) · "
        "[Figure SVG](../figures/network_modules/"+STEM+".svg) · "
        "[Source data](../tables/network_modules/annotation_significance_plot_data.csv)\n",encoding="utf-8")

    # Confirm that presentation preserves the saved probabilities and all rows.
    reread = pd.read_csv(data_file, float_precision="round_trip").set_index("module_id")
    original = source.set_index("module_id")
    for field in ["P_value", "BH_all_modules_and_terms", "candidate_background_BH"]:
        assert np.allclose(reread.loc[original.index,field], original[field],rtol=1e-13,atol=1e-30),field
    pdf = PdfReader(FIGURES / f"{STEM}.pdf")
    assert len(pdf.pages)==1
    page = pdf.pages[0]
    assert abs(float(page.mediabox.width)*25.4/72-183)<.01
    assert abs(float(page.mediabox.height)*25.4/72-235)<.01
    chunks=[]
    page.extract_text(visitor_text=lambda text,cm,tm,font,size: chunks.append(text))
    text="\n".join(chunks)
    for row in data.itertuples():
        assert re.search(r"\b"+row.display_id+r"\b",text),row.display_id
        assert probability_text(row.P_value) in text
        assert probability_text(row.BH_all_modules_and_terms) in text
    fonts=set()
    for ref in page["/Resources"]["/Font"].values():
        font=ref.get_object()
        for item in font.get("/DescendantFonts",[font]):
            f=item.get_object(); descriptor=f["/FontDescriptor"].get_object()
            assert any(k in descriptor for k in ["/FontFile","/FontFile2","/FontFile3"])
            assert "Arial" in str(f.get("/BaseFont"))
            fonts.add(str(f.get("/BaseFont")))
    with Image.open(FIGURES / f"{STEM}.png") as image:
        assert abs(image.width-183/25.4*600)<1.1 and abs(image.height-235/25.4*600)<1.1
        assert all(abs(d-600)<1 for d in image.info["dpi"])
    svg=(FIGURES / f"{STEM}.svg").read_text(encoding="utf-8")
    assert "<text" in svg and "Arial" in svg
    assert before=={p.relative_to(ROOT).as_posix():sha(p) for p in inputs}
    outputs=[FIGURES/f"{STEM}.{ext}" for ext in ["pdf","svg","png"]]+[data_file,caption_file]
    result={"completed_utc":datetime.now(timezone.utc).isoformat(),"modules":42,"lower_RNA_modules":24,"higher_RNA_modules":18,
            "raw_P_values":42,"BH_q_values":42,"annotation_q_below_005":int(data.annotation_q_below_005.sum()),
            "unchanged_input_sha256":before,"script_sha256":sha(__file__),"annotation_tests_per_family":260526,
            "new_tests_performed":False,"embedded_fonts":sorted(fonts),"caption_words":len(caption.split()),
            "exports":{p.relative_to(ROOT).as_posix():sha(p) for p in outputs},"visual_review":"pending"}
    (ROOT/"provenance/network_modules/annotation_significance_completion.json").write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ["exports","unchanged_input_sha256"]},indent=2))


if __name__=="__main__":
    main()
