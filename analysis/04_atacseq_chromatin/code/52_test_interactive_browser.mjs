import {spawn} from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath,pathToFileURL} from 'node:url';
import assert from 'node:assert/strict';
const root=path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const portable=process.argv.includes('--portable');
const drivers=process.argv.includes('--drivers');
const out=path.join(root,drivers?'integration/pathway_GSEA/gene_drivers':portable?'integration/portable_atac_rna':'integration/interactive_atac_rna');
fs.mkdirSync(out,{recursive:true});
const packageRoot=process.argv.includes('--package')?path.resolve(process.argv[process.argv.indexOf('--package')+1]):path.join(root,'exports/HT1080_ATAC_RNA_portable');
const sleep=ms=>new Promise(r=>setTimeout(r,ms));
const chrome=spawn('C:/Program Files/Google/Chrome/Application/chrome.exe',[
  '--headless=new','--no-first-run','--no-default-browser-check','--disable-background-networking','--disable-component-update',
  '--remote-debugging-address=127.0.0.1','--remote-debugging-port=9231',
  '--user-data-dir='+path.join(out,'browser_test_profile'),'--use-angle=swiftshader','--enable-unsafe-swiftshader','about:blank'
],{windowsHide:true,stdio:['ignore','ignore','pipe']});
const chromeLog=fs.createWriteStream(path.join(out,'browser_test_chrome.log'));
chrome.stderr.pipe(chromeLog);
chrome.on('error',e=>{console.error(e);process.exit(1);});
let target;
for(let i=0;i<60;i++){try{target=await (await fetch('http://127.0.0.1:9231/json/new?about:blank',{method:'PUT'})).json();break;}catch{await sleep(500);}}
if(!target)throw Error('Headless Chrome did not start');
const ws=new WebSocket(target.webSocketDebuggerUrl);await new Promise((resolve,reject)=>{ws.onopen=resolve;ws.onerror=reject;});
let seq=0;const pending=new Map(),errors=[],networkRequests=[];
ws.onmessage=e=>{const m=JSON.parse(e.data);if(m.id){const p=pending.get(m.id);if(p){pending.delete(m.id);clearTimeout(p.timer);m.error?p.reject(Error(JSON.stringify(m.error))):p.resolve(m.result);}}else if(m.method==='Runtime.exceptionThrown')errors.push(m.params.exceptionDetails);else if(m.method==='Network.requestWillBeSent'&&/^https?:/.test(m.params.request.url))networkRequests.push(m.params.request.url);};
ws.onclose=e=>{for(const p of pending.values()){clearTimeout(p.timer);p.reject(Error('CDP closed: '+e.code+' '+e.reason));}pending.clear();};
function cdp(method,params={}){return new Promise((resolve,reject)=>{if(ws.readyState!==1){reject(Error('CDP not open'));return;}const id=++seq;const timer=setTimeout(()=>{pending.delete(id);reject(Error('CDP timeout: '+method));},20000);pending.set(id,{resolve,reject,timer});ws.send(JSON.stringify({id,method,params}));});}
async function evaluate(expression){const r=await cdp('Runtime.evaluate',{expression,returnByValue:true,awaitPromise:true});if(r.exceptionDetails)throw Error(JSON.stringify(r.exceptionDetails));return r.result.value;}
async function until(expression,label){for(let i=0;i<120;i++){if(await evaluate(expression))return;await sleep(250);}throw Error('Timeout: '+label+'; '+await evaluate('document.body.innerText.slice(-1000)'));}
try{
  await cdp('Runtime.enable');await cdp('Page.enable');
  if(portable||drivers){
    await cdp('Network.enable');
    await cdp('Network.emulateNetworkConditions',{offline:true,latency:0,downloadThroughput:0,uploadThroughput:0});
    await cdp('Page.addScriptToEvaluateOnNewDocument',{source:`window.fetch=()=>{throw Error('Portable report must not use fetch or a server.');};`});
  }
  await cdp('Emulation.setDeviceMetricsOverride',{width:1500,height:1100,deviceScaleFactor:1,mobile:false});
  if(drivers){
    const {testDrivers}=await import('./66_test_gene_drivers.mjs');
    await testDrivers({cdp,evaluate,until,out,root,errors,networkRequests});
  }else{
  await cdp('Page.navigate',{url:portable?pathToFileURL(path.join(packageRoot,'START_HERE.html')).href:'http://127.0.0.1:8765/report/interactive_ATAC_RNA.html'});
  await until('window.viewerReady && window.viewerState?.locus?.gene === "25850"','initial viewer');
  const initial=await evaluate('({plotly:Plotly.version,title:document.getElementById("gene-heading").textContent,summary:document.getElementById("scatter-summary").textContent,traces:document.getElementById("coverage").data.length,expression:document.getElementById("expression").data.reduce((n,t)=>n+t.x.length,0)})');
  if(!initial.summary.includes('82 overlap')||initial.expression!==12)throw Error(JSON.stringify(initial));
  // Use an actual mouse click on the PDPN scatter point, not a direct selection function call.
  const point=await evaluate(`(()=>{const gd=document.getElementById('scatter');let p;for(const t of gd.data.slice(0,4)){const i=t.customdata.findIndex(d=>d[0]==='10630');if(i>=0)p={x:t.x[i],y:t.y[i]};}const b=gd.getBoundingClientRect(),l=gd._fullLayout;return{x:b.x+l.xaxis._offset+l.xaxis.l2p(p.x),y:b.y+l.yaxis._offset+l.yaxis.l2p(p.y)};})()`);
  await cdp('Input.dispatchMouseEvent',{type:'mouseMoved',x:point.x,y:point.y});
  await cdp('Input.dispatchMouseEvent',{type:'mousePressed',x:point.x,y:point.y,button:'left',clickCount:1});
  await cdp('Input.dispatchMouseEvent',{type:'mouseReleased',x:point.x,y:point.y,button:'left',clickCount:1});
  await until('window.viewerState.gene === "10630" && window.viewerState.locus?.gene === "10630"','scatter point click');
  await evaluate(`document.getElementById('dataset').value='AUG_2025';document.getElementById('dataset').dispatchEvent(new Event('change'));`);
  await until('window.viewerState.dataset === "AUG_2025" && document.getElementById("scatter-summary").textContent.includes("29 overlap") && window.viewerState.locus?.gene === "10630"','August switch');
  const six=await evaluate('document.getElementById("expression").data.reduce((n,t)=>n+t.x.length,0)');if(six!==6)throw Error('August expression sample count');
  await evaluate(`document.getElementById('dataset').value='MAY_2026';document.getElementById('dataset').dispatchEvent(new Event('change'));`);
  await until('document.getElementById("scatter-summary").textContent.includes("83 overlap") && window.viewerState.locus?.gene === "10630"','May switch');
  await evaluate(`document.getElementById('feature').value='gene_body';document.getElementById('feature').dispatchEvent(new Event('change'));`);
  await until('window.viewerState.locus?.feature === "gene_body" && document.getElementById("scatter-summary").textContent.includes("115 overlap")','feature switch');
  // Select a gene outside the prepared atlas to verify genuinely on-demand tracks.
  const outside=await evaluate(`(()=>{const row=ATAC_RNA_DATA.panels['MAY_2026:gene_body'].rows.find(r=>!ATAC_RNA_DATA.genes[r[0]].atlas);document.getElementById('gene-search').value=row[0];document.getElementById('find-gene').click();return row[0];})()`);
  await until(`window.viewerState.locus?.gene === '${outside}'`,'non-atlas gene coverage');
  await evaluate(`document.getElementById('locus-view').value='overview';document.getElementById('locus-view').dispatchEvent(new Event('change'));`);
  await until('window.viewerState.view === "overview" && window.viewerState.locus !== null','gene overview');
  const view=await evaluate('({gene:viewerState.gene,atlas:ATAC_RNA_DATA.genes[viewerState.gene].atlas,coverageSamples:viewerState.locus.coverage.length,linkedPeaks:viewerState.locus.linked_peaks.length,models:viewerState.locus.genes.length})');
  if(view.atlas||view.coverageSamples!==6)throw Error(JSON.stringify(view));
  // Return to a readable promoter example and capture both the scatter/expression and locus sections.
  await evaluate(`document.getElementById('dataset').value='Combined';document.getElementById('dataset').dispatchEvent(new Event('change'));`);
  await until('viewerState.locus !== null','combined return');
  await evaluate(`document.getElementById('feature').value='promoter';document.getElementById('feature').dispatchEvent(new Event('change'));`);
  await until('viewerState.locus?.feature === "promoter"','promoter return');
  await evaluate(`document.getElementById('gene-search').value='ZNF345';document.getElementById('find-gene').click();`);
  await until('viewerState.locus?.gene === "25850"','final example');await sleep(750);
  await cdp('Input.dispatchMouseEvent',{type:'mouseMoved',x:20,y:20});
  await evaluate('window.scrollTo(0,0)');
  let screenshot=await cdp('Page.captureScreenshot',{format:'png'});fs.writeFileSync(path.join(out,'interactive_scatter_preview.png'),Buffer.from(screenshot.data,'base64'));
  await evaluate('document.getElementById("coverage").scrollIntoView()');await sleep(250);
  screenshot=await cdp('Page.captureScreenshot',{format:'png'});fs.writeFileSync(path.join(out,'interactive_tracks_preview.png'),Buffer.from(screenshot.data,'base64'));
  if(errors.length)throw Error(JSON.stringify(errors));
  const anchor=await evaluate('document.getElementById("coverage").layout.xaxis.anchor');if(anchor!=='y9')throw Error('Coverage genomic axis is not below gene models');
  let portableValidation=null;
  if(portable){
    const fixtures=JSON.parse(fs.readFileSync(path.join(out,'reference_loci.json'),'utf8'));
    let maxError=0;
    for(const f of fixtures){
      const actual=await evaluate(`window.ATAC_PORTABLE.locus(${JSON.stringify(f.gene)},${JSON.stringify(f.feature)},${JSON.stringify(f.view)})`);
      assert.equal(actual.coverage.length,6);
      for(let s=0;s<6;s++){
        assert.equal(actual.coverage[s].length,f.expected.coverage[s].length);
        for(let j=0;j<actual.coverage[s].length;j++){
          const error=Math.abs(actual.coverage[s][j]-f.expected.coverage[s][j]);maxError=Math.max(maxError,error);
          assert.ok(error<1.1e-7,`Coverage mismatch ${f.gene} ${f.view} sample ${s} bin ${j}: ${error}`);
        }
      }
      const {coverage:unusedA,...metaA}=actual,{coverage:unusedB,...metaB}=f.expected;
      assert.deepStrictEqual(metaA,metaB,`Locus metadata mismatch ${f.gene} ${f.view}`);
    }
    // Explicit peak selection as well as the default/overview aliases.
    await evaluate(`loadLocus(viewerState.locus.linked_peaks.at(-1).id)`);
    assert.ok(await evaluate('viewerState.locus !== null && viewerState.view === viewerState.locus.selected_peak'));
    assert.equal(await evaluate('location.protocol'),'file:');
    assert.deepStrictEqual(networkRequests,[],'Portable viewer attempted HTTP requests');
    if(errors.length)throw Error(JSON.stringify(errors));
    portableValidation={mode:'file:// with browser network offline and fetch disabled',reference_loci:fixtures.length,chromosomes:new Set(fixtures.map(f=>f.expected.chrom)).size,max_coverage_absolute_error:maxError,HTTP_requests:networkRequests.length,explicit_peak_selection:'passed',package_root:packageRoot};
  }
  const result={status:'passed',initial,actual_mouse_click_gene:'PDPN',August_expression_samples:six,non_atlas_test:view,checks:['Actual scatter click updates gene, RNA and ATAC panels','Combined/August/May overlap counts match static reports','Feature switching works','All plotted genes can load tracks, including non-atlas genes','Overview and peak views work','Genomic axis is beneath the gene models','No browser JavaScript exceptions']};
  if(portable)result.portable_validation=portableValidation;
  fs.writeFileSync(path.join(out,'browser_verification.json'),JSON.stringify(result,null,2));console.log(JSON.stringify(result,null,2));
  }
}finally{await cdp('Browser.close').catch(()=>{});ws.close();}
