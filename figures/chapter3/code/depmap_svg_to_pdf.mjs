// Render a local SVG with embedded bundled Open Sans into an editable-text PDF.
import { readFile } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import path from 'node:path';
import { createRequire } from 'node:module';

const [svgPath, pdfPath, fontDir] = process.argv.slice(2);
if (!svgPath || !pdfPath || !fontDir) throw new Error('Expected SVG, PDF and bundled font directory');
const modules = process.env.CODEX_NODE_MODULES || path.join(process.env.USERPROFILE,
  '.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules');
const require = createRequire(path.join(modules, '_depmap_pdf_renderer.cjs'));
const { chromium } = require('playwright');
let svg = await readFile(svgPath, 'utf8');
const width = svg.match(/\bwidth=['"]([\d.]+)pt['"]/)[1];
const height = svg.match(/\bheight=['"]([\d.]+)pt['"]/)[1];
let fonts = '';
for (const [file, weight, style] of [
  ['OpenSans-Regular.ttf', '400', 'normal'],
  ['OpenSans-SemiBold.ttf', '700', 'normal'],
  ['OpenSans-Italic.ttf', '400', 'italic'],
  ['OpenSans-SemiBoldItalic.ttf', '700', 'italic'],
]) {
  const data = (await readFile(path.join(fontDir, file))).toString('base64');
  fonts += `@font-face{font-family:'Open Sans';font-weight:${weight};font-style:${style};src:url(data:font/ttf;base64,${data}) format('truetype');}`;
}
svg = svg.replace(/<\?xml[^>]*\?>/, '');
const html = `<!doctype html><html><head><style>${fonts}
@page{size:${width}pt ${height}pt;margin:0}html,body{margin:0;padding:0;width:${width}pt;height:${height}pt}
svg{display:block}</style></head><body>${svg}</body></html>`;
const edge = path.join(process.env['ProgramFiles(x86)'], 'Microsoft/Edge/Application/msedge.exe');
const browser = await chromium.launch({headless:true, ...(existsSync(edge) ? {executablePath:edge} : {})});
try {
  const page = await browser.newPage();
  await page.route('**/*', route => route.abort());
  await page.setContent(html, {waitUntil:'load'});
  await page.evaluate(async () => { await document.fonts.ready; });
  const fontLoaded = await page.evaluate(() => document.fonts.check('9px "Open Sans"'));
  if (!fontLoaded) throw new Error('Bundled Open Sans failed to load');
  await page.pdf({path:pdfPath, preferCSSPageSize:true, printBackground:true, margin:{top:0,right:0,bottom:0,left:0}});
  console.log(`PDF exported with embedded Open Sans: ${pdfPath}`);
} finally {
  await browser.close();
}
