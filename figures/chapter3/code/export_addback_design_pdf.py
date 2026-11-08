"""Export the saved add-back SVG to PDF with bundled, embedded Open Sans fonts.

Windows Cairo substitutes Arial on this host. Run after the R renderer:
    python figures/chapter3/code/export_addback_design_pdf.py
Only the local SVG and bundled fonts are rendered; no data analysis is rerun.
"""

import base64
import argparse
import os
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET

from pypdf import PdfReader


def main():
    root = Path(__file__).resolve().parents[3]
    figure_dir = root / "figures/chapter3"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stem', nargs='?', default='ch03_main_ht1080_addback_design_expression',
                        help='Existing chapter3 SVG stem, without extension')
    stem = parser.parse_args().stem
    if not re.fullmatch(r'ch03_[A-Za-z0-9_]+', stem):
        raise ValueError('Expected a chapter3 figure filename stem')
    font_dir = root / "skills/thesis-figure-style/assets/fonts"
    temp_dir = root / ".tmp"
    temp_dir.mkdir(exist_ok=True)
    css = []
    for face, weight, style in [
        ("Regular", 400, "normal"),
        ("SemiBold", 700, "normal"),
        ("Italic", 400, "italic"),
        ("SemiBoldItalic", 700, "italic"),
    ]:
        encoded = base64.b64encode((font_dir / f"OpenSans-{face}.ttf").read_bytes()).decode("ascii")
        css.append(
            "@font-face{font-family:'Open Sans';"
            f"font-weight:{weight};font-style:{style};"
            f"src:url(data:font/ttf;base64,{encoded}) format('truetype');}}"
        )
    svg = (figure_dir / f"{stem}.svg").read_text(encoding="utf-8")
    svg = re.sub(r"<\?xml.*?\?>", "", svg, flags=re.S)
    # The SVG viewBox uses points; follow the current figure dimensions.
    _, _, width, height = map(float, ET.fromstring(svg).attrib['viewBox'].split())
    page_css = (
        f"@page{{size:{width}pt {height}pt;margin:0}}"
        f"html,body{{width:{width}pt;height:{height}pt;margin:0;padding:0;overflow:hidden;"
        "print-color-adjust:exact;-webkit-print-color-adjust:exact}"
        f"svg{{display:block;width:{width}pt;height:{height}pt}}"
    )
    html = (
        '<!DOCTYPE html><html><head><meta charset="UTF-8"><style>'
        + "".join(css)
        + page_css
        + "</style></head><body>"
        + svg
        + "</body></html>"
    )
    html_file = temp_dir / f"{stem}_pdf_render.html"
    html_file.write_text(html, encoding="utf-8")
    browser_candidates = [
        Path(os.environ.get("PROGRAMFILES", "C:/Program Files")) / "Google/Chrome/Application/chrome.exe",
        Path(os.environ.get("PROGRAMFILES(X86)", "C:/Program Files (x86)")) / "Microsoft/Edge/Application/msedge.exe",
    ]
    browser = next((p for p in browser_candidates if p.exists()), None)
    if browser is None:
        raise RuntimeError("Chrome or Edge is required for this Windows PDF export.")
    pdf_file = figure_dir / f"{stem}.pdf"
    result = subprocess.run(
        [
            str(browser), "--headless", "--no-sandbox", "--disable-gpu",
            "--disable-software-rasterizer", "--no-first-run", "--no-default-browser-check",
            f"--user-data-dir={temp_dir / 'addback_pdf_chrome_profile'}",
            f"--print-to-pdf={pdf_file}", "--no-pdf-header-footer",
            "--virtual-time-budget=3000", html_file.as_uri(),
        ],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        creationflags=subprocess.CREATE_NO_WINDOW, timeout=45,
    )
    if result.returncode:
        raise RuntimeError(result.stderr[-2000:])
    pdf = PdfReader(pdf_file)
    assert len(pdf.pages) == 1
    fonts = pdf.pages[0]["/Resources"]["/Font"].get_object()
    for reference in fonts.values():
        font = reference.get_object()
        assert "OpenSans" in str(font["/BaseFont"]), font["/BaseFont"]
        leaf = font["/DescendantFonts"][0].get_object() if "/DescendantFonts" in font else font
        descriptor = leaf["/FontDescriptor"].get_object()
        assert any(
            key in descriptor and len(descriptor[key].get_object().get_data()) > 0
            for key in ("/FontFile", "/FontFile2", "/FontFile3")
        )
    print(f"Exported {pdf_file.name}: one page, embedded Open Sans.")


if __name__ == "__main__":
    main()
