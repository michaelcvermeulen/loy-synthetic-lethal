import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import {pathToFileURL} from 'node:url';
export async function testDrivers({cdp,evaluate,until,out,root,errors,networkRequests}){
  await cdp('Page.navigate',{url:pathToFileURL(path.join(root,'report/pathway_gene_drivers.html')).href});
  await until('window.driverReady===true','driver viewer ready');
  const check=()=>evaluate(`(()=>{const g=PATHWAY_DRIVERS.genes[driverState.gene],p=g.peaks.find(p=>p.id===driverState.peak),e=document.getElementById('driver-expression').data,a=document.getElementById('driver-atac-counts').data,f=document.getElementById('driver-effects').data[0];return {gene:driverState.gene,rnaSamples:e.reduce((n,t)=>n+t.x.length,0),atacSamples:a.reduce((n,t)=>n+t.x.length,0),rnaIntervals:f.error_x.array.length,rawRNA:e.flatMap(t=>t.customdata.map(d=>d[1])),rawATAC:a.flatMap(t=>t.customdata.map(d=>d[1])),expectedATAC:p.raw,shownEffects:f.x,expectedEffects:g.rna.map(r=>r.logFC)};})()`);
  let state=await check();assert.equal(state.rnaSamples,12);assert.equal(state.atacSamples,6);assert.equal(state.rnaIntervals,3);assert.deepEqual(state.shownEffects,state.expectedEffects);assert.deepEqual(state.rawATAC,state.expectedATAC);
  for(const key of ['MYC_V1','MYC_V2','IFN_alpha','IFN_gamma']){
    await evaluate(`document.getElementById('driver-pathway').value='${key}';document.getElementById('driver-pathway').dispatchEvent(new Event('change'));`);
    await until(`driverState.pathway==='${key}' && document.getElementById('driver-status').textContent.startsWith('Loaded '+PATHWAY_DRIVERS.genes[driverState.gene].symbol)`,'pathway '+key);
    const count=await evaluate('document.getElementById("driver-gene").options.length');assert.ok(count>20);
  }
  // Actual mouse click on a different gene row in the current pathway table.
  const target=await evaluate(`(()=>{const row=[...document.querySelectorAll('#driver-table tr[data-gene]')].find(r=>r.dataset.gene!==driverState.gene);row.scrollIntoView({block:'center'});const b=row.getBoundingClientRect();return {gene:row.dataset.gene,x:b.x+30,y:b.y+b.height/2};})()`);
  await cdp('Input.dispatchMouseEvent',{type:'mousePressed',x:target.x,y:target.y,button:'left',clickCount:1});
  await cdp('Input.dispatchMouseEvent',{type:'mouseReleased',x:target.x,y:target.y,button:'left',clickCount:1});
  await until(`driverState.gene==='${target.gene}' && document.getElementById('driver-status').textContent.startsWith('Loaded '+PATHWAY_DRIVERS.genes['${target.gene}'].symbol)`,'gene table click');
  const multi=await evaluate(`Object.keys(PATHWAY_DRIVERS.genes).find(g=>PATHWAY_DRIVERS.genes[g].peaks.length>1)`);
  await evaluate(`document.getElementById('driver-search').value='${multi}';document.getElementById('driver-find').click();`);
  await until(`driverState.gene==='${multi}' && document.getElementById('driver-status').textContent.startsWith('Loaded '+PATHWAY_DRIVERS.genes['${multi}'].symbol)`,'gene search');
  await evaluate(`const s=document.getElementById('driver-peak');s.value=s.options[s.options.length-1].value;s.dispatchEvent(new Event('change'));`);
  state=await check();assert.deepEqual(state.rawATAC,state.expectedATAC);assert.deepEqual(state.shownEffects,state.expectedEffects);
  await evaluate(`document.getElementById('driver-search').value='BST2';document.getElementById('driver-find').click();`);
  await until(`document.getElementById('driver-status').textContent.startsWith('Loaded BST2')`,'BST2 example');
  await evaluate(`document.getElementById('driver-heading').scrollIntoView({block:'start'})`);
  const shot=await cdp('Page.captureScreenshot',{format:'png'});fs.writeFileSync(path.join(out,'driver_browser_preview.png'),Buffer.from(shot.data,'base64'));
  assert.deepEqual(errors,[]);assert.deepEqual(networkRequests,[]);
  const result={status:'passed',mode:'file:// with networking disabled',genes:await evaluate('Object.keys(PATHWAY_DRIVERS.genes).length'),pathways:5,RNA_samples:12,ATAC_samples:6,checks:['All five pathway selectors','Actual gene-table mouse click','Gene search across pathways','Peak selection updates sample counts','Displayed RNA effects match source payload','No JavaScript errors or HTTP requests']};
  fs.writeFileSync(path.join(out,'browser_verification.json'),JSON.stringify(result,null,2));console.log(JSON.stringify(result,null,2));
}
