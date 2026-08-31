import"./styles-B8WF9G-F.js";import{p as I}from"./parseXYZ-ByfvMHmk.js";import{g as G,h as U,z as Ce,o as Re,f as Le,n as Ie,C as Ae,F as je,M as Ne}from"./ri-DT3jfd0v.js";import{c as Oe}from"./gradient-Bmg7fDyY.js";import{c as Te}from"./hessian-DpHI2yIv.js";import"./properties-Odo1T7fO.js";const He="modulepreload",ze=function(e){return"/GANSU-Lite/"+e},ae={},Pe=function(t,n,s){let o=Promise.resolve();if(n&&n.length>0){let u=function(l){return Promise.all(l.map(c=>Promise.resolve(c).then(b=>({status:"fulfilled",value:b}),b=>({status:"rejected",reason:b}))))};document.getElementsByTagName("link");const r=document.querySelector("meta[property=csp-nonce]"),i=r?.nonce||r?.getAttribute("nonce");o=u(n.map(l=>{if(l=ze(l),l in ae)return;ae[l]=!0;const c=l.endsWith(".css"),b=c?'[rel="stylesheet"]':"";if(document.querySelector(`link[href="${l}"]${b}`))return;const d=document.createElement("link");if(d.rel=c?"stylesheet":He,c||(d.as="script"),d.crossOrigin="",d.href=l,i&&d.setAttribute("nonce",i),document.head.appendChild(d),c)return new Promise((p,h)=>{d.addEventListener("load",p),d.addEventListener("error",()=>h(new Error(`Unable to preload CSS for ${l}`)))})}))}function a(r){const i=new Event("vite:preloadError",{cancelable:!0});if(i.payload=r,window.dispatchEvent(i),!i.defaultPrevented)throw r}return o.then(r=>{for(const i of r||[])i.status==="rejected"&&a(i.reason);return t().catch(a)})};function Q(){const e=navigator;return{ua:e.userAgent,platform:e.platform??"unknown",deviceMemoryGB:e.deviceMemory??null,hardwareConcurrency:e.hardwareConcurrency??0,screen:`${screen.width}×${screen.height}`,language:e.language||"unknown",timestamp:new Date().toISOString()}}function pe(){try{const e=performance.getEntriesByType("navigation")[0];return e?Math.round(e.loadEventEnd-e.fetchStart):null}catch{return null}}const N=.9572,O=104.52*Math.PI/180,T=3,We=[["O",0,0,0],["H",N*Math.sin(O/2),0,N*Math.cos(O/2)],["H",-N*Math.sin(O/2),0,N*Math.cos(O/2)]];function Ge(e){const t=Math.ceil(Math.cbrt(e))+1,n=[];for(let s=-t;s<=t;s++)for(let o=-t;o<=t;o++)for(let a=-t;a<=t;a++)n.push([s,o,a]);return n.sort((s,o)=>{const a=s[0]*s[0]+s[1]*s[1]+s[2]*s[2],r=o[0]*o[0]+o[1]*o[1]+o[2]*o[2];return a!==r?a-r:s[2]!==o[2]?s[2]-o[2]:s[1]!==o[1]?s[1]-o[1]:s[0]-o[0]}),n.slice(0,e)}function Ue(e){const t=Math.PI*(3-Math.sqrt(5));return[e*t,e*1.2345678,e*.7654321]}function Je(e,[t,n,s]){const o=(r,i)=>[r[0]*Math.cos(i)-r[1]*Math.sin(i),r[0]*Math.sin(i)+r[1]*Math.cos(i),r[2]];return o(((r,i)=>[r[0]*Math.cos(i)+r[2]*Math.sin(i),r[1],-r[0]*Math.sin(i)+r[2]*Math.cos(i)])(o([e[0],e[1],e[2]],s),n),t)}function De(e){const t=Ge(e),n=[];for(let s=0;s<e;s++){const[o,a,r]=t[s],i=Ue(s),u=o*T,l=a*T,c=r*T;for(const[b,d,p,h]of We){const[y,v,g]=Je([d,p,h],i);n.push(`${b.padEnd(2)} ${(u+y).toFixed(6).padStart(12)} ${(l+v).toFixed(6).padStart(12)} ${(c+g).toFixed(6).padStart(12)}`)}}return`${e*3}
(H2O)${e} — cubic lattice a=${T} A, golden-angle orientations
${n.join(`
`)}`}const re=1.526,B=1.09,qe=109.47*Math.PI/180;function Ye(e){const t=qe/2,n=re*Math.sin(t),s=re*Math.cos(t),o=[];for(let c=0;c<e;c++)o.push([c*n,0,c%2===0?0:s]);const a=[];for(const c of o)a.push(_("C",c));for(let c=0;c<e;c++){const[b,d,p]=o[c],h=c%2===0?1:-1,y=p+h*B*Math.cos(t),v=B*Math.sin(t);a.push(_("H",[b,d+v,y])),a.push(_("H",[b,d-v,y]))}const r=o[0],i=o[e-1],u=[-Math.sin(t),0,e>1&&o[1][2]>r[2]?-Math.cos(t):Math.cos(t)],l=[Math.sin(t),0,e>1&&o[e-2][2]>i[2]?-Math.cos(t):Math.cos(t)];return a.push(_("H",[r[0]+B*u[0],r[1],r[2]+B*u[2]])),a.push(_("H",[i[0]+B*l[0],i[1],i[2]+B*l[2]])),`${e+2*e+2}
C${e}H${2*e+2} all-anti
${a.join(`
`)}`}function _(e,[t,n,s]){return`${e.padEnd(2)} ${t.toFixed(6).padStart(12)} ${n.toFixed(6).padStart(12)} ${s.toFixed(6).padStart(12)}`}const E=[{id:"water",label:"Water cluster (H₂O)ₙ",name:e=>`(H2O)${e}`,natoms:e=>3*e,xyz:De,ladder:[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32]},{id:"alkane",label:"n-Alkane CₙH₂ₙ₊₂",name:e=>`C${e}H${2*e+2}`,natoms:e=>3*e+2,xyz:Ye,ladder:[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20]}],ee="gansu-stress-v1",te="gansu-stress-probe-v1",be=180,X="6-31g(d,p)",Ve=3e3;function Xe(){try{const e=localStorage.getItem(ee);return e?JSON.parse(e):null}catch{return null}}function K(e){try{localStorage.setItem(ee,JSON.stringify(e))}catch{}}function Ke(){try{localStorage.removeItem(ee)}catch{}}function k(e){return e==null?"—":e<1024**2?`${(e/1024).toFixed(0)} KB`:e<1024**3?`${(e/1024**2).toFixed(0)} MB`:`${(e/1024**3).toFixed(2)} GB`}function fe(e){return e==null?"—":(e/1e3).toFixed(2)}const se={ok:"ok",timeout:"timeout (exceeded budget)",error:"error / OOM","worker-died":"worker killed (OOM)","tab-crash":"tab crash"};function ie(e,t,n,s,o,a){return new Promise(r=>{const i=new Worker(new URL("/GANSU-Lite/assets/stressWorker-XVKoLPu9.js",import.meta.url),{type:"module"});a(i);let u=!1;const l=performance.now();let c=l,b="no progress reported",d=null;const p=w=>{u||(u=!0,clearTimeout(v),i.terminate(),a(null),r(w))},h=w=>{const f=w*(w+1)/2;return f*(f+1)/2*8},y=()=>({nbasis:d,totalMs:null,scfMs:null,iterations:null,converged:null,energy:null,eriBytes:d==null?null:h(d),heapUsedBytes:null,heapLimitBytes:null}),v=setTimeout(()=>{const w=(performance.now()-c)/1e3;p({...y(),status:"timeout",totalMs:performance.now()-l,detail:`exceeded ${(n/1e3).toFixed(0)} s budget; silent for ${w.toFixed(0)} s after "${b}"`+(w>30?" — the worker had gone quiet, so it may have been killed rather than merely slow":" — the worker was still reporting progress, so this is genuine slowness, not a crash")})},n);i.onmessage=w=>{const f=w.data;if(c=performance.now(),f.type==="progress"){const S=/(\d+) basis functions/.exec(f.message);S&&(d=Number(S[1])),b=f.message,o(f.message,f.elapsedMs);return}if(f.type==="done"){p({status:"ok",nbasis:f.nbasis,totalMs:f.totalMs,scfMs:f.scfMs,iterations:f.iterations,converged:f.converged,energy:f.energy,eriBytes:f.eriBytes,heapUsedBytes:f.heapUsedBytes,heapLimitBytes:f.heapLimitBytes,...f.converged?{}:{detail:"SCF did not converge"}});return}f.type==="error"&&p({...y(),status:"error",nbasis:f.nbasis??d,totalMs:f.elapsedMs,detail:`${f.name}: ${f.message} (phase: ${f.phase})`})},i.onerror=w=>{p({...y(),status:"worker-died",totalMs:performance.now()-l,detail:w.message||"worker terminated by the browser"})};const g={type:"stress-run",xyzText:e,basisGBS:t,charge:0,eriBackend:"auto",baseUrl:"/GANSU-Lite/",...s?{dftConfig:{functional:"B3LYP",gridLevel:"medium"}}:{}};i.postMessage(g)})}function Ze(){const e=new Set;for(const t of E[0].ladder)e.add(t*25);return[...e].filter(t=>t<=300).sort((t,n)=>t-n)}function Qe(e){return new Promise(t=>{const n=new Worker(new URL("/GANSU-Lite/assets/stressWorker-XVKoLPu9.js",import.meta.url),{type:"module"}),s=[],o=()=>{n.terminate(),t(s)};n.onmessage=r=>{const i=r.data;if(i.type==="probe-result"){const u={nbasis:i.nbasis,bytes:i.bytes,peakBytes:i.peakBytes,ok:i.ok,ms:i.ms,verdict:i.verdict,...i.failedStage?{failedStage:i.failedStage}:{},...i.detail?{detail:i.detail}:{}};s.push(u),e(u);return}i.type==="probe-done"&&o()},n.onerror=()=>o();const a={type:"memory-probe",nbasisLadder:Ze(),baseUrl:"/GANSU-Lite/",holdMs:Ve};n.postMessage(a)})}function he(e,t){const n=t.verdict==="ok",s=t.verdict==="ok"?"held: wasm array + JS copy":t.verdict==="unavailable"?"WASM unavailable":t.failedStage==="js-copy"?"FAILED on the JS copy — device memory cap":"FAILED in wasm — address space",o=document.createElement("tr");o.innerHTML=`<td class="num">${t.nbasis||"—"}</td><td class="num">${t.bytes?k(t.bytes):"—"}</td><td class="num">${t.peakBytes?k(t.peakBytes):"—"}</td><td class="num">${t.ms.toFixed(0)}</td><td${n?"":' class="cb-err"'}>${s}${t.detail?`<br><span class="cb-err">${ge(t.detail)}</span>`:""}</td>`,e.appendChild(o)}let x=[];function et(e){try{localStorage.setItem(te,JSON.stringify(e))}catch{}}function tt(){try{const e=localStorage.getItem(te);return e?JSON.parse(e):[]}catch{return[]}}async function st(){const e=document.getElementById("cb-stress-probe"),t=document.getElementById("cb-probe-body"),n=document.getElementById("cb-probe-section");e&&(e.disabled=!0),n&&(n.style.display=""),t&&(t.innerHTML=""),x=[],M("probing memory ceiling…"),await Qe(r=>{x.push(r),et(x),t&&he(t,r)});const s=x.filter(r=>r.verdict==="ok"),o=x.find(r=>r.verdict!=="ok"),a=document.getElementById("cb-probe-summary");if(a){const r=s.length?s[s.length-1]:null,i=x.find(u=>u.verdict==="unavailable");if(i)a.textContent=`No WASM on this device, so there is no wasm32 ceiling to measure — it runs the JS backend, whose limit is device RAM and shows up in the ladder above. ${i.detail??""}`;else if(!r)a.textContent="Could not allocate even the smallest probe array — something is wrong.";else{const u=o?o.failedStage==="js-copy"?`${o.nbasis} failed while taking the JS copy — this device's own memory cap, not the wasm32 address space, is what stops it (needed ${k(o.peakBytes)})`:`${o.nbasis} failed inside wasm — the module's 4 GiB address space ran out (needed ${k(o.bytes)} in one block)`:"the probe ladder topped out at 300 basis functions without hitting a ceiling";a.textContent=`Memory ceiling: ${o?"":"≥ "}${r.nbasis} basis functions — the SCF's peak of ${k(r.peakBytes)} (a ${k(r.bytes)} ERI array in wasm plus its JS copy) is held successfully; ${u}.`}}M("memory probe done"),m.output&&ne(),e&&(e.disabled=!1)}function ye(e,t){let n=null;for(const s of e)if(!(s.status!=="ok"||s.totalMs==null)){if(s.totalMs/1e3>t)break;n=s}return n}function nt(e,t){if(e.length<2||t.nbasis==null)return null;const n=e[e.length-2],s=e[e.length-1];if(n.nbasis==null||s.nbasis==null||n.totalMs==null||s.totalMs==null||n.nbasis>=s.nbasis||s.nbasis>=t.nbasis||n.totalMs<=0)return null;const o=Math.log(s.totalMs/n.totalMs)/Math.log(s.nbasis/n.nbasis);return!Number.isFinite(o)||o<=0?null:{seconds:s.totalMs*Math.pow(t.nbasis/s.nbasis,o)/1e3,exponent:o}}function ot(e,t){const n=e.filter(r=>r.status==="ok"),s=n.length?n[n.length-1]:null,o=e.find(r=>r.status!=="ok")??null;let a=null;if(o&&o.status==="timeout"){const r=nt(n,o);r&&(a={estimatedSeconds:r.seconds,budgetSeconds:t,exponent:r.exponent,fragile:r.seconds<1.25*t})}return{borderline:a,largestCompleted:s&&s.nbasis!=null&&s.totalMs!=null?{name:s.name,natoms:s.natoms,nbasis:s.nbasis,seconds:+(s.totalMs/1e3).toFixed(1)}:null,firstFailure:o?{name:o.name,natoms:o.natoms,nbasis:o.nbasis,mode:se[o.status]}:null}}let m={series:E[0],points:[],running:!1,stopRequested:!1,activeWorker:null,output:null};function at(){return`
    <section class="cb-panel">
      <h2>Stress test — how large a system does this device survive?</h2>
      <p class="cb-note">
        Runs a monotonically growing series at RHF/6-31G(d,p) until the device fails, and records
        the failure mode. Each point runs in a fresh Web Worker with a wall-clock budget, so an
        out-of-memory abort or a run that is simply too slow is captured instead of taking the page
        down. Progress is checkpointed to this browser's local storage after every point — if the
        tab does die, reload and the lost point is reported as a tab crash.
      </p>
      <div id="cb-stress-recovered" class="cb-note" style="display:none;color:var(--color-error,#e05050)"></div>
      <div class="cb-stress-controls">
        <label>Series
          <select id="cb-stress-series">
            ${E.map(e=>`<option value="${e.id}">${e.label}</option>`).join("")}
          </select>
        </label>
        <label>Theory
          <select id="cb-stress-theory">
            <option value="rhf">RHF/6-31G(d,p)</option>
            <option value="b3lyp">B3LYP/6-31G(d,p)</option>
          </select>
        </label>
        <label>Budget per point (s)
          <input id="cb-stress-budget" type="number" min="10" max="3600" step="10" value="${be}">
        </label>
        <label>Ladder (n)
          <input id="cb-stress-ladder" type="text" value="${E[0].ladder.join(", ")}">
        </label>
      </div>
      <div style="margin:10px 0 8px">
        <button id="cb-stress-run">Run stress test</button>
        <button id="cb-stress-stop" disabled>Stop</button>
        <button id="cb-stress-reset">Clear saved progress</button>
        <span id="cb-stress-status" class="cb-wl-status" style="margin-left:10px"></span>
      </div>
      <div style="overflow-x:auto">
        <table class="cb-results-table" style="font-size:0.8rem">
          <thead><tr>
            <th>System</th><th class="num">Atoms</th><th class="num">Basis fns</th>
            <th class="num">ERI array</th><th class="num">Iters</th>
            <th class="num">Time (s)</th><th>Status</th>
          </tr></thead>
          <tbody id="cb-stress-body"></tbody>
        </table>
      </div>
      <h3 style="font-size:0.95rem;margin:16px 0 4px">Hard memory ceiling</h3>
      <p class="cb-note">
        The SCF's peak footprint is <em>two</em> copies of the unique-ERI array (N⁴/8 doubles each):
        one inside the WASM module's linear memory, plus the JS <code>Float64Array</code> it is
        sliced into — and wasm memory never shrinks, so both stay resident. The probe reproduces
        exactly that pair, at each size, and holds it resident for 3 s — seconds instead of the tens
        of minutes an SCF would take, but long enough that a pressure-based tab killer gets the same
        chance at it. It reports which half gave out, because the two answer different questions:
        failing <em>inside wasm</em> means the 4 GiB address space ran out (what stops a desktop),
        while failing <em>on the JS copy</em> means the device's own memory cap did (what stops a
        phone). Synthetic substitutes do not work — a lone JS array or a bare
        <code>WebAssembly.Memory</code> both claim success where the real calculation dies.
      </p>
      <button id="cb-stress-probe">Probe memory ceiling</button>
      <div id="cb-probe-section" style="display:none;margin-top:10px">
        <table class="cb-results-table" style="font-size:0.8rem;max-width:520px">
          <thead><tr><th class="num">Basis fns</th><th class="num">ERI array</th><th class="num">Peak held</th><th class="num">ms</th><th>Result</th></tr></thead>
          <tbody id="cb-probe-body"></tbody>
        </table>
        <p id="cb-probe-summary" class="cb-note"></p>
      </div>

      <div id="cb-stress-summary" style="display:none">
        <h3 style="font-size:0.95rem;margin:14px 0 6px">Row for the paper table</h3>
        <label style="font-size:0.78rem;color:var(--color-text-dim);display:block;margin-bottom:8px">
          Practical limit measured at a patience threshold of
          <input id="cb-stress-threshold" type="number" min="1" max="3600" step="10" value="180"
                 style="width:80px;background:var(--color-input);color:var(--color-text);border:1px solid var(--color-border);border-radius:5px;padding:3px 6px">
          s — applied to the recorded timings, so every device can be compared at the
          same threshold no matter what budget its ladder ran with.
        </label>
        <pre id="cb-stress-row" class="cb-stress-row"></pre>
        <button id="cb-stress-copy-row">Copy row</button>
        <button id="cb-stress-copy-json">Copy JSON</button>
        <button id="cb-stress-download-json">Download JSON</button>
        <button id="cb-stress-download-csv">Download CSV</button>
      </div>
    </section>
  `}function rt(e){const t=a=>document.getElementById(a),n=t("cb-stress-series"),s=t("cb-stress-ladder");n?.addEventListener("change",()=>{const a=E.find(r=>r.id===n.value);a&&s&&(m.series=a,s.value=a.ladder.join(", "))}),t("cb-stress-run")?.addEventListener("click",()=>{dt(e)}),t("cb-stress-stop")?.addEventListener("click",()=>{m.stopRequested=!0,m.activeWorker?.terminate(),M("stopping…")}),t("cb-stress-reset")?.addEventListener("click",()=>{Ke();try{localStorage.removeItem(te)}catch{}x=[],m.points=[],m.output=null,W();const a=t("cb-stress-recovered");a&&(a.style.display="none");const r=t("cb-stress-summary");r&&(r.style.display="none"),M("saved progress cleared")}),t("cb-stress-probe")?.addEventListener("click",()=>{st()}),t("cb-stress-threshold")?.addEventListener("input",ne),t("cb-stress-copy-row")?.addEventListener("click",()=>ce(ve())),t("cb-stress-copy-json")?.addEventListener("click",()=>ce(JSON.stringify(m.output,null,2))),t("cb-stress-download-json")?.addEventListener("click",()=>le("stress-test.json",JSON.stringify(m.output,null,2),"application/json")),t("cb-stress-download-csv")?.addEventListener("click",()=>le("stress-test.csv",mt(m.output),"text/csv")),it();const o=tt();if(o.length){x=o;const a=document.getElementById("cb-probe-section"),r=document.getElementById("cb-probe-body");if(a&&r){a.style.display="";for(const i of o)he(r,i)}}}function it(){const e=Xe();if(!e)return;m.points=e.points??[];const t=E.find(n=>n.id===e.seriesId);if(t&&(m.series=t),e.inFlight){const n=e.inFlight;if(!m.points.some(o=>o.n===n.n)){const o=n.nbasis??null,a=o==null?null:o*(o+1)/2;m.points.push({n:n.n,name:n.name,natoms:n.natoms,nbasis:o,status:"tab-crash",totalMs:null,scfMs:null,iterations:null,converged:null,energy:null,eriBytes:a==null?null:a*(a+1)/2*8,heapUsedBytes:null,heapLimitBytes:null,detail:"tab died while this point was running"})}const s=document.getElementById("cb-stress-recovered");s&&(s.style.display="",s.textContent=`Recovered from a previous session: the tab died while running ${n.name}. Results up to that point are restored below and ${n.name} is recorded as a tab crash.`),K({...e,points:m.points,inFlight:null})}m.points.length&&W()}function M(e){const t=document.getElementById("cb-stress-status");t&&(t.textContent=e)}async function ct(){const e=await fetch(`/GANSU-Lite/basis/${X}.gbs`);if(!e.ok)throw new Error(`Cannot load basis ${X}: HTTP ${e.status}`);return e.text()}function lt(e,t){const n=e.split(/[,\s]+/).map(s=>parseInt(s,10)).filter(s=>Number.isFinite(s)&&s>0);return n.length?n:t}async function dt(e){if(m.running)return;const t=document.getElementById("cb-stress-run"),n=document.getElementById("cb-stress-stop"),s=Math.max(10,Number(document.getElementById("cb-stress-budget")?.value)||be),o=document.getElementById("cb-stress-theory")?.value==="b3lyp",a=lt(document.getElementById("cb-stress-ladder")?.value??"",m.series.ladder);m.running=!0,m.stopRequested=!1,m.points=[],m.output=null,t&&(t.disabled=!0),n&&(n.disabled=!1),W();let r;try{M("loading basis set…"),r=await ct()}catch(c){M(`FAILED: ${c instanceof Error?c.message:String(c)}`),m.running=!1,t&&(t.disabled=!1),n&&(n.disabled=!0);return}const i=m.series,u=o?"B3LYP":"RHF",l=s*1e3;try{M("warm-up…"),await ie(i.xyz(a[0]),r,l,o,()=>{},c=>{m.activeWorker=c})}catch{}for(const c of a){if(m.stopRequested){M("stopped");break}const b=i.name(c),d=i.natoms(c),p={n:c,name:b,natoms:d,startedAt:Date.now(),nbasis:null},h=()=>K({seriesId:i.id,method:u,budgetSeconds:s,points:m.points,inFlight:p});h(),M(`${b} — starting…`);const y=await ie(i.xyz(c),r,l,o,(g,w)=>{if(p.nbasis==null){const f=/(\d+) basis functions/.exec(g);f&&(p.nbasis=Number(f[1]),h())}M(`${b} — ${g} (${(w/1e3).toFixed(1)} s)`)},g=>{m.activeWorker=g}),v={n:c,name:b,natoms:d,...y};if(m.points.push(v),W(),K({seriesId:i.id,method:u,budgetSeconds:s,points:m.points,inFlight:null}),y.status!=="ok"){M(`${b} failed: ${se[y.status]} — stopping ladder`);break}M(`${b} ok — ${fe(y.totalMs)} s`),await new Promise(g=>setTimeout(g,200))}m.output={device:Q(),backend:e(),seriesId:i.id,seriesLabel:i.label,method:u,basis:X,budgetSeconds:s,points:m.points,summary:ot(m.points,s),memoryProbe:x,toolUrl:location.href},ne(),m.running=!1,m.activeWorker=null,t&&(t.disabled=!1),n&&(n.disabled=!0),m.stopRequested||M("done")}function W(){const e=document.getElementById("cb-stress-body");if(e){e.innerHTML="";for(const t of m.points){const n=document.createElement("tr"),s=t.status!=="ok";n.innerHTML=`
      <td>${t.name}</td>
      <td class="num">${t.natoms}</td>
      <td class="num">${t.nbasis??"—"}</td>
      <td class="num">${k(t.eriBytes)}</td>
      <td class="num">${t.iterations??"—"}</td>
      <td class="num">${fe(t.totalMs)}</td>
      <td${s?' class="cb-err"':""}>${se[t.status]}${t.detail?`<br><span class="cb-err">${ge(t.detail)}</span>`:""}</td>
    `,e.appendChild(n)}}}function ge(e){return e.replace(/[&<>"]/g,t=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"})[t])}function we(){const e=document.getElementById("cb-stress-threshold"),t=Number(e?.value);return Number.isFinite(t)&&t>0?t:180}function ve(){const e=m.output;if(!e)return"";const t=we(),n=ye(e.points,t),s=e.summary.largestCompleted,o=e.summary.firstFailure,a=x.filter(l=>l.verdict==="ok"),r=a.length?a[a.length-1].nbasis:null,i=x.find(l=>l.verdict==="failed");return`| ${[ut(e.device),n?n.name:"—",n?String(n.natoms):"—",n?.nbasis!=null?String(n.nbasis):"—",n&&n.totalMs!=null?(n.totalMs/1e3).toFixed(1):"—",s&&(n?.nbasis==null||s.nbasis>n.nbasis)?`${s.name} (${s.seconds}s)`:"—",o?`${o.name}: ${o.mode}`:"(no failure in ladder)",r==null?"(probe not run)":`${i?"":"≥"}${r}`].join(" | ")} |`}function ut(e){const t=e.ua;return/iPhone/.test(t)?"iPhone":/iPad/.test(t)?"iPad":/Mac OS X|Macintosh/.test(t)?"Apple-silicon laptop":/Android/.test(t)?"Android":/Windows/.test(t)?"Windows":e.platform}function ne(){const e=document.getElementById("cb-stress-summary"),t=document.getElementById("cb-stress-row");if(!e||!t||!m.output)return;e.style.display="";const n=m.output.summary,s=we(),o=`| Device | Largest within ${s}s | Atoms | Basis fns | Time (s) | Ran further (over threshold) | Ladder ended | WASM ceiling |`,a="|---|---|--:|--:|--:|---|---|--:|",r=[];if(ye(m.output.points,s)||r.push(`# Nothing completed within ${s} s — lower the threshold or read the timings directly.`),n.firstFailure||r.push("# The ladder finished without a failure — extend it to find where this device stops."),n.borderline){const l=n.borderline;r.push(`# ${n.firstFailure?.name} is extrapolated to need ~${l.estimatedSeconds.toFixed(0)} s (fitted N^${l.exponent.toFixed(1)}) against the ${l.budgetSeconds} s run budget.`),l.fragile&&r.push("# FRAGILE: that is within 25% of the budget, so run-to-run noise alone could flip it.",`# The "${s}s" column is unaffected — it is read off the recorded timings — but the`,`# ladder stopped early, so re-run at ${Math.ceil(l.estimatedSeconds*1.5/60)*60} s if you need rungs beyond this one.`)}const u=m.points.filter(l=>l.status==="ok"&&l.converged===!1);u.length&&r.push(`# SCF did not converge for: ${u.map(l=>l.name).join(", ")}`),x.length||r.push('# Memory ceiling not measured — press "Probe memory ceiling".'),t.textContent=[o,a,ve(),...r].join(`
`)}function ce(e){navigator.clipboard.writeText(e).then(()=>alert("Copied to clipboard")).catch(()=>alert("Copy failed — open browser console"))}function le(e,t,n){const s=new Blob([t],{type:n}),o=URL.createObjectURL(s),a=document.createElement("a");a.href=o,a.download=e,a.click(),setTimeout(()=>URL.revokeObjectURL(o),1e3)}function mt(e){if(!e)return"";const t=["n","system","atoms","nbasis","eri_bytes","iterations","converged","setup_plus_scf_ms","scf_ms","energy_hartree","status","detail"],n=e.points.map(o=>[o.n,o.name,o.natoms,o.nbasis??"",o.eriBytes??"",o.iterations??"",o.converged==null?"":String(o.converged),o.totalMs?.toFixed(2)??"",o.scfMs?.toFixed(2)??"",o.energy?.toFixed(8)??"",o.status,o.detail??""].map(a=>`"${String(a).replace(/"/g,'""')}"`).join(","));return[`# device.ua,"${e.device.ua.replace(/"/g,'""')}"`,`# device.platform,"${e.device.platform}"`,`# device.deviceMemoryGB,"${e.device.deviceMemoryGB??""}"`,`# device.hardwareConcurrency,"${e.device.hardwareConcurrency}"`,`# backend,"${e.backend}"`,`# series,"${e.seriesLabel}"`,`# method,"${e.method}/${e.basis}"`,`# budget_seconds,"${e.budgetSeconds}"`,`# timestamp,"${e.device.timestamp}"`].join(`
`)+`
`+t.join(",")+`
`+n.join(`
`)+`
`}const D=4,pt=1,q=1.8897259886,Me=`3
H2O standard
O   0.000000   0.000000   0.117790
H   0.000000   0.755450  -0.471160
H   0.000000  -0.755450  -0.471160`,bt=`3
H2O displaced
O   0.000000   0.000000   0.000000
H   0.000000   0.850000   0.620000
H   0.000000  -0.850000   0.620000`,ft=`3
CO2 linear
C   0.000000   0.000000   0.000000
O   0.000000   0.000000   1.160000
O   0.000000   0.000000  -1.160000`,Y=new Map;async function A(e){if(Y.has(e))return Y.get(e);const t=`/GANSU-Lite/basis/${e}.gbs`,n=await fetch(t);if(!n.ok)throw new Error(`Failed to load basis ${e}: ${n.status}`);const s=Le.fromGBS(await n.text());return Y.set(e,s),s}async function ht(e){const t=await A("6-31g(d,p)"),n=I(Me),s=new G(n,t,0);await U(s,"RHF").solve({eriBackend:e})}async function yt(e){const t=await A("sto-3g"),n=I(bt);let s=new Float64Array(n.length*3);for(let r=0;r<n.length;r++)s[3*r]=n[r].coordinate.x,s[3*r+1]=n[r].coordinate.y,s[3*r+2]=n[r].coordinate.z;const o=.4,a=25;for(let r=0;r<a;r++){const i=n.map((d,p)=>({...d,coordinate:{x:s[3*p],y:s[3*p+1],z:s[3*p+2]}})),u=new G(i,t,0),l=U(u,"RHF");await l.solve({eriBackend:e});const c=Oe(u.primitiveShells,u.atoms,u.cgtoNormalizationFactors,u.numBasis,u.numAlphaSpins,l.density,l.coefficients,l.orbitalEnergies).total;let b=0;for(let d=0;d<c.length;d++)Math.abs(c[d])>b&&(b=Math.abs(c[d]));if(b<5e-4)break;for(let d=0;d<c.length;d++)s[d]-=o*c[d]}}async function gt(e){const t=await A("sto-3g"),n=I(ft),s=new Float64Array(n.length*3);for(let o=0;o<n.length;o++)s[3*o]=n[o].coordinate.x*q,s[3*o+1]=n[o].coordinate.y*q,s[3*o+2]=n[o].coordinate.z*q;await Te(n.map(o=>o.atomicNumber),s,t,0,5e-4,void 0,void 0)}async function wt(e){const t=await A("6-31g(d,p)"),n=I(Me),s=new G(n,t,0);await U(s,"RHF",{functional:"B3LYP",gridLevel:"medium"}).solve({eriBackend:e})}const xe=[{id:"h2o_hf_631gdp",label:"H2O HF/6-31G(d,p)",run:ht},{id:"h2o_geomopt_sto3g",label:"H2O geom opt HF/STO-3G",run:yt},{id:"co2_hess_sto3g",label:"CO2 Hessian HF/STO-3G",run:gt},{id:"h2o_b3lyp_631gdp",label:"H2O B3LYP/6-31G(d,p)",run:wt}];function Z(e){const t=[...e].sort((s,o)=>s-o),n=t.length;return n===0?NaN:n%2===1?t[(n-1)/2]:(t[n/2-1]+t[n/2])/2}async function Se(e){const t=await fetch(e);if(!t.ok)throw new Error(`Cannot fetch ${e}`);const n=await t.arrayBuffer(),s=await WebAssembly.compile(n);let o;const a={"./wasm_eri_bg.js":{__wbindgen_init_externref_table:()=>{const d=o.exports.__wbindgen_externrefs,p=d.grow(4);d.set(0,void 0),d.set(p+0,void 0),d.set(p+1,null),d.set(p+2,!0),d.set(p+3,!1)}}};o=await WebAssembly.instantiate(s,a);const r=o.exports;r.__wbindgen_start&&r.__wbindgen_start();const i=r.memory,u=r.__wbindgen_malloc,l=r.__wbindgen_free,c=d=>{const p=u(d.length*8,8)>>>0;return new Float64Array(i.buffer).set(d,p/8),[p,d.length]},b=d=>{const[p,h]=d,y=new Float64Array(i.buffer).slice(p/8,p/8+h);return l(p,h*8,8),y};return{computeERIs:(d,p,h,y)=>{const[v,g]=c(d),[w,f]=c(p);return b(r.compute_eris_wasm(v,g,w,f,h,y))},computeFockRhf:(d,p,h,y)=>{const[v,g]=c(d),[w,f]=c(p),[S,j]=c(h);return b(r.compute_fock_rhf(v,g,w,f,S,j,y))}}}let C,R;async function vt(){if(C!==void 0)return C;try{C=await Se("/GANSU-Lite/wasm/wasm_eri_bg.wasm")}catch{C=null}return C}async function Mt(){if(R!==void 0)return R;try{R=await Se("/GANSU-Lite/wasm/wasm_eri_simd_bg.wasm")}catch{R=null}return R}async function $e(){return await Ce("/GANSU-Lite/"),Re()?Ie()==="wasm-simd"?{backend:"wasm",label:"WASM-SIMD"}:{backend:"wasm",label:"WASM"}:{backend:"js",label:"JS"}}function de(e,t,n,s){let o=e,a=t;if(o<a){const c=o;o=a,a=c}let r=n,i=s;if(r<i){const c=r;r=i,i=c}let u=o*(o+1)/2+a,l=r*(r+1)/2+i;if(u<l){const c=u;u=l,l=c}return u*(u+1)/2+l}function xt(e,t,n,s){const o=new Float64Array(s*s);for(let a=0;a<s;a++)for(let r=a;r<s;r++){let i=0;for(let l=0;l<s;l++)for(let c=0;c<s;c++){const b=t[l*s+c];Math.abs(b)<1e-10||(i+=b*(e[de(a,r,l,c)]-.5*e[de(a,l,r,c)]))}const u=n[a*s+r]+i;o[a*s+r]=u,o[r*s+a]=u}return o}async function St(e,t,n,s){const o=e.numBasis,a=new je(t,o),r=U(e,"RHF");r.computeNuclearRepulsionEnergy(),r.computeCoreHamiltonianMatrix(),r.computeTransformMatrix(),r.eri=a,r.guessInitialFockMatrix();const i=r.coreHamiltonianMatrix.data;let u=0,l=0,c=0;const b=performance.now();for(let d=0;d<200;d++){r.computeCoefficientMatrix(),r.computeDensityMatrix();const p=r.densityMatrix.data,h=await s(t,p,i,o);if(r.fockMatrix=new Ne(o,o,h),c=r.computeEnergy(),d>0&&Math.abs(c-u)<1e-8){l=d+1;break}r.updateFockMatrix(),u=c}return{eriMs:n,scfMs:performance.now()-b,iters:l,energy:c}}const V=[{molecule:"H2O",xyzFile:"H2O.xyz",basis:"6-31g(d,p)"},{molecule:"CO2",xyzFile:"CO2.xyz",basis:"6-31g(d,p)"},{molecule:"Benzene",xyzFile:"Benzene.xyz",basis:"6-31g(d,p)"}],ue=3;async function $t(e){const t=await A(e.basis),n=`/GANSU-Lite/xyz/${e.xyzFile}`,s=await fetch(n);if(!s.ok)throw new Error(`Cannot load ${n}`);const o=await s.text(),a=I(o);return{mol:new G(a,t,0),xyz:o}}async function kt(e,t,n,s){const{mol:o}=await $t(e),a=o.numBasis,r=o.primitiveShells,i=o.cgtoNormalizationFactors,u=(await Pe(async()=>{const{packShells:h}=await import("./ri-DT3jfd0v.js").then(y=>y.W);return{packShells:h}},[])).packShells(r),l=new Float64Array(i);async function c(h,y,v){const g=[],w=[];let f=0,S;try{for(let $=0;$<ue;$++){s(`  ${h} run ${$+1}/${ue}…`);const Be=performance.now(),Ee=await y(),Fe=performance.now()-Be,J=await St(o,Ee,Fe,v);g.push(J.eriMs),w.push(J.scfMs),$===0&&(f=J.iters),await new Promise(_e=>setTimeout(_e,30))}}catch($){S=$ instanceof Error?$.message:String($)}const j=g.length?Z(g):NaN,oe=w.length?Z(w):NaN;return{eriMedianMs:j,scfMedianMs:oe,totalMedianMs:j+oe,iters:f,eriRuns:g,scfRuns:w,...S?{error:S}:{}}}const b=await c("JS",()=>Ae(r,i,a,1e-10,void 0,"js"),(h,y,v,g)=>xt(h,y,v,g)),d=t?await c("WASM baseline",()=>t.computeERIs(u,l,a,1e-10),(h,y,v,g)=>t.computeFockRhf(h,y,v,g)):{eriMedianMs:NaN,scfMedianMs:NaN,totalMedianMs:NaN,iters:0,eriRuns:[],scfRuns:[],error:"WASM baseline not available"},p=n?await c("WASM+SIMD",()=>n.computeERIs(u,l,a,1e-10),(h,y,v,g)=>n.computeFockRhf(h,y,v,g)):{eriMedianMs:NaN,scfMedianMs:NaN,totalMedianMs:NaN,iters:0,eriRuns:[],scfRuns:[],error:"SIMD not available"};return{molecule:e.molecule,basis:e.basis,nbasis:a,js:b,wasmBase:d,wasmSimd:p}}let z=null;async function Bt(){const e=document.getElementById("cb-run-compare"),t=document.getElementById("cb-stop-compare");e&&(e.disabled=!0),t&&(t.disabled=!1);const n=r=>{const i=document.getElementById("cb-compare-status");i&&(i.textContent=r)};n("Loading WASM binaries…");const s=await vt(),o=await Mt(),a=[];for(let r=0;r<V.length;r++){const i=V[r];n(`[${r+1}/${V.length}] ${i.molecule}/${i.basis}…`);try{a.push(await kt(i,s,o,n)),Et(a)}catch(u){const l=u instanceof Error?u.message:String(u);n(`  ERROR: ${l}`)}}z=a,n("Done."),e&&(e.disabled=!1),t&&(t.disabled=!0)}function Et(e){const t=document.getElementById("cb-compare-section");if(!t)return;t.style.display="";const n=document.getElementById("cb-compare-body");if(!n)return;n.innerHTML="";const s=a=>isFinite(a)?a.toFixed(0):"—",o=(a,r)=>isFinite(a)&&isFinite(r)&&r>0?`${(a/r).toFixed(1)}×`:"—";for(const a of e){const r=document.createElement("tr");r.innerHTML=`
      <td>${a.molecule}</td><td>${a.basis}</td><td class="num">${a.nbasis}</td><td class="num">${a.js.iters||"—"}</td>
      <td class="num">${s(a.js.eriMedianMs)}</td>
      <td class="num">${s(a.js.scfMedianMs)}</td>
      <td class="num"><b>${s(a.js.totalMedianMs)}</b></td>
      <td class="num">${s(a.wasmBase.eriMedianMs)}</td>
      <td class="num">${s(a.wasmBase.scfMedianMs)}</td>
      <td class="num"><b>${s(a.wasmBase.totalMedianMs)}</b></td>
      <td class="num">${o(a.js.totalMedianMs,a.wasmBase.totalMedianMs)}</td>
      <td class="num">${s(a.wasmSimd.eriMedianMs)}</td>
      <td class="num">${s(a.wasmSimd.scfMedianMs)}</td>
      <td class="num"><b>${s(a.wasmSimd.totalMedianMs)}</b></td>
      <td class="num">${o(a.js.totalMedianMs,a.wasmSimd.totalMedianMs)}</td>
    `,n.appendChild(r)}}const Ft=document.getElementById("app");let F=null,me="detecting…";function _t(){const e=Q(),t=pe();Ft.innerHTML=`
    <div class="cb-page">
      <header class="cb-header">
        <h1>GANSU Lite — Cross-device Benchmark</h1>
        <p class="cb-subtitle">Reproducible end-to-end timing for the paper. Median of 3 timed runs (+ 1 warmup discarded).</p>
      </header>

      <section class="cb-panel">
        <h2>Device & Backend</h2>
        <table class="cb-info">
          <tr><th>User-Agent</th><td>${e.ua}</td></tr>
          <tr><th>Platform</th><td>${e.platform}</td></tr>
          <tr><th>Device memory</th><td>${e.deviceMemoryGB?`${e.deviceMemoryGB} GB`:"unknown"}</td></tr>
          <tr><th>Cores (navigator)</th><td>${e.hardwareConcurrency}</td></tr>
          <tr><th>Screen</th><td>${e.screen}</td></tr>
          <tr><th>Initial load</th><td>${t!=null?`${t} ms`:"unavailable"}</td></tr>
          <tr><th>Backend</th><td id="cb-backend">detecting…</td></tr>
        </table>
      </section>

      <section class="cb-panel">
        <h2>Workloads</h2>
        <p class="cb-note">All fixed: SCF tol = 1e-7, DIIS, default initial guess. Single Web Worker thread (no parallelism — single-core perf + memory bandwidth).</p>
        <ul class="cb-workloads">
          ${xe.map(n=>`<li id="cb-row-${n.id}"><span class="cb-wl-label">${n.label}</span><span class="cb-wl-status">pending</span></li>`).join("")}
        </ul>
        <button id="cb-run">Run all workloads</button>
        <button id="cb-stop" disabled>Stop</button>
      </section>

      <section class="cb-panel">
        <h2>Backend comparison (reference device)</h2>
        <p class="cb-note">
          SCF wall-clock time for H₂O, CO₂, and benzene at 6-31G(d,p), broken down by kernel
          (ERI computation + SCF iteration loop), comparing TypeScript / baseline WebAssembly /
          WebAssembly +simd128 on this device. Medians over 3 runs each.
        </p>
        <div style="margin:6px 0 8px">
          <button id="cb-run-compare">Run backend comparison</button>
          <button id="cb-stop-compare" disabled>Stop</button>
          <span id="cb-compare-status" style="margin-left:10px;font-size:0.85rem;color:var(--color-text-dim);font-family:'Cascadia Code',monospace;"></span>
        </div>
        <div id="cb-compare-section" style="display:none;overflow-x:auto;">
          <table class="cb-results-table" style="font-size:0.78rem;">
            <thead><tr>
              <th rowspan="2">Molecule</th><th rowspan="2">Basis</th>
              <th rowspan="2" class="num">N</th><th rowspan="2" class="num">Iters</th>
              <th colspan="3" style="text-align:center;border-bottom:1px solid var(--color-border)">TypeScript (JS)</th>
              <th colspan="4" style="text-align:center;border-bottom:1px solid var(--color-border)">Baseline WASM</th>
              <th colspan="4" style="text-align:center;border-bottom:1px solid var(--color-border)">WASM +simd128</th>
            </tr><tr>
              <th class="num">ERI</th><th class="num">SCF</th><th class="num">Total</th>
              <th class="num">ERI</th><th class="num">SCF</th><th class="num">Total</th><th class="num">↑vs JS</th>
              <th class="num">ERI</th><th class="num">SCF</th><th class="num">Total</th><th class="num">↑vs JS</th>
            </tr></thead>
            <tbody id="cb-compare-body"></tbody>
          </table>
          <p class="cb-note">All times in ms. ↑vs JS = TypeScript total / WASM total.</p>
          <div style="margin-top:8px">
            <button id="cb-copy-compare-json">Copy JSON</button>
            <button id="cb-download-compare-json">Download JSON</button>
            <button id="cb-download-compare-csv">Download CSV</button>
          </div>
        </div>
      </section>

      ${at()}

      <section class="cb-panel" id="cb-results" style="display:none">
        <h2>Results</h2>
        <table class="cb-results-table">
          <thead><tr>
            <th>Workload</th>
            <th class="num">Warmup (ms)</th>
            <th class="num">Run 1 (ms)</th>
            <th class="num">Run 2 (ms)</th>
            <th class="num">Run 3 (ms)</th>
            <th class="num">Median (ms)</th>
          </tr></thead>
          <tbody id="cb-results-body"></tbody>
        </table>
        <div class="cb-export">
          <button id="cb-copy-json">Copy JSON</button>
          <button id="cb-download-json">Download JSON</button>
          <button id="cb-download-csv">Download CSV</button>
        </div>
        <details class="cb-json">
          <summary>Raw JSON</summary>
          <pre id="cb-json-pre"></pre>
        </details>
      </section>
    </div>
  `,jt(),document.getElementById("cb-run").addEventListener("click",Ct),document.getElementById("cb-stop").addEventListener("click",()=>{P=!0}),document.getElementById("cb-run-compare")?.addEventListener("click",Bt),document.getElementById("cb-copy-compare-json")?.addEventListener("click",()=>ke(JSON.stringify(z,null,2))),document.getElementById("cb-download-compare-json")?.addEventListener("click",()=>H("backend-comparison.json",JSON.stringify(z,null,2),"application/json")),document.getElementById("cb-download-compare-csv")?.addEventListener("click",()=>H("backend-comparison.csv",It(z),"text/csv")),document.getElementById("cb-copy-json")?.addEventListener("click",Lt),document.getElementById("cb-download-json")?.addEventListener("click",()=>H("benchmark.json",JSON.stringify(F,null,2),"application/json")),document.getElementById("cb-download-csv")?.addEventListener("click",()=>H("benchmark.csv",At(F),"text/csv")),rt(()=>me),$e().then(({label:n})=>{me=n;const s=document.getElementById("cb-backend");s&&(s.textContent=n)})}let P=!1;function L(e,t){const n=document.getElementById(`cb-row-${e}`);if(n){const s=n.querySelector(".cb-wl-status");s&&(s.textContent=t)}}async function Ct(){P=!1;const e=document.getElementById("cb-run"),t=document.getElementById("cb-stop");e.disabled=!0,t.disabled=!1;const{backend:n,label:s}=await $e(),o=Q(),a=pe(),r=[];for(const i of xe){if(P){L(i.id,"stopped");continue}const u=[];let l;try{for(let d=0;d<D&&!P;d++){L(i.id,`run ${d+1}/${D}…`);const p=performance.now();await i.run(n);const h=performance.now()-p;u.push(h),L(i.id,`run ${d+1}/${D}: ${h.toFixed(0)} ms`),await new Promise(y=>setTimeout(y,50))}}catch(d){l=d instanceof Error?d.message:String(d),L(i.id,`ERROR: ${l}`)}const c=u.slice(pt),b=c.length>0?Z(c):NaN;r.push({id:i.id,label:i.label,warmupMs:u[0]??NaN,runs:c,medianMs:b,...l?{error:l}:{}}),!l&&c.length>0&&L(i.id,`median ${b.toFixed(0)} ms`)}F={device:o,backend:s,initialLoadMs:a,workloads:r,toolUrl:location.href},Rt(F),e.disabled=!1,t.disabled=!0}function Rt(e){const t=document.getElementById("cb-results");t.style.display="";const n=document.getElementById("cb-results-body");n.innerHTML="";for(const s of e.workloads){const o=document.createElement("tr");o.innerHTML=`
      <td>${s.label}${s.error?` <span class="cb-err">(${s.error})</span>`:""}</td>
      <td class="num">${isFinite(s.warmupMs)?s.warmupMs.toFixed(0):"—"}</td>
      <td class="num">${s.runs[0]!=null?s.runs[0].toFixed(0):"—"}</td>
      <td class="num">${s.runs[1]!=null?s.runs[1].toFixed(0):"—"}</td>
      <td class="num">${s.runs[2]!=null?s.runs[2].toFixed(0):"—"}</td>
      <td class="num"><b>${isFinite(s.medianMs)?s.medianMs.toFixed(0):"—"}</b></td>
    `,n.appendChild(o)}document.getElementById("cb-json-pre").textContent=JSON.stringify(e,null,2)}function Lt(){F&&ke(JSON.stringify(F,null,2))}function ke(e){navigator.clipboard.writeText(e).then(()=>alert("Copied to clipboard")).catch(()=>alert("Copy failed — open browser console"))}function It(e){if(!e||e.length===0)return"";const t=["molecule","basis","nbasis","iters","js_eri_ms","js_scf_ms","js_total_ms","wasm_eri_ms","wasm_scf_ms","wasm_total_ms","speedup_wasm_vs_js","simd_eri_ms","simd_scf_ms","simd_total_ms","speedup_simd_vs_js"],n=a=>isFinite(a)?a.toFixed(2):"",s=(a,r)=>isFinite(a)&&isFinite(r)&&r>0?(a/r).toFixed(2):"",o=e.map(a=>[a.molecule,a.basis,String(a.nbasis),String(a.js.iters),n(a.js.eriMedianMs),n(a.js.scfMedianMs),n(a.js.totalMedianMs),n(a.wasmBase.eriMedianMs),n(a.wasmBase.scfMedianMs),n(a.wasmBase.totalMedianMs),s(a.js.totalMedianMs,a.wasmBase.totalMedianMs),n(a.wasmSimd.eriMedianMs),n(a.wasmSimd.scfMedianMs),n(a.wasmSimd.totalMedianMs),s(a.js.totalMedianMs,a.wasmSimd.totalMedianMs)].map(r=>`"${r}"`).join(","));return t.join(",")+`
`+o.join(`
`)+`
`}function H(e,t,n){const s=new Blob([t],{type:n}),o=URL.createObjectURL(s),a=document.createElement("a");a.href=o,a.download=e,a.click(),setTimeout(()=>URL.revokeObjectURL(o),1e3)}function At(e){if(!e)return"";const t=["workload_id","workload_label","backend","warmup_ms","run1_ms","run2_ms","run3_ms","median_ms","error"],n=e.workloads.map(o=>[o.id,o.label,e.backend,isFinite(o.warmupMs)?o.warmupMs.toFixed(2):"",o.runs[0]?.toFixed(2)??"",o.runs[1]?.toFixed(2)??"",o.runs[2]?.toFixed(2)??"",isFinite(o.medianMs)?o.medianMs.toFixed(2):"",o.error??""].map(a=>`"${String(a).replace(/"/g,'""')}"`).join(","));return[`# device.ua,"${e.device.ua.replace(/"/g,'""')}"`,`# device.platform,"${e.device.platform}"`,`# device.deviceMemoryGB,"${e.device.deviceMemoryGB??""}"`,`# device.hardwareConcurrency,"${e.device.hardwareConcurrency}"`,`# device.screen,"${e.device.screen}"`,`# backend,"${e.backend}"`,`# initialLoadMs,"${e.initialLoadMs??""}"`,`# timestamp,"${e.device.timestamp}"`].join(`
`)+`
`+t.join(",")+`
`+n.join(`
`)+`
`}function jt(){if(document.getElementById("cb-styles"))return;const e=document.createElement("style");e.id="cb-styles",e.textContent=`
    .cb-page { max-width: 900px; margin: 24px auto; padding: 0 16px; color: var(--color-text); }
    .cb-header h1 { margin: 0 0 4px; font-size: 1.6rem; }
    .cb-subtitle { color: var(--color-text-dim); font-size: 0.92rem; margin-bottom: 24px; }
    .cb-panel { background: var(--color-surface); border: 1px solid var(--color-border); border-radius: 8px; padding: 16px 20px; margin-bottom: 16px; }
    .cb-panel h2 { font-size: 1.1rem; margin: 0 0 12px; color: var(--color-accent); }
    .cb-note { font-size: 0.82rem; color: var(--color-text-dim); margin: 0 0 12px; }
    table.cb-info { width: 100%; font-size: 0.85rem; }
    table.cb-info th { text-align: left; padding: 4px 12px 4px 0; color: var(--color-text-dim); font-weight: 500; vertical-align: top; width: 160px; }
    table.cb-info td { padding: 4px 0; font-family: 'Cascadia Code', monospace; font-size: 0.82rem; word-break: break-all; }
    .cb-workloads { list-style: none; padding: 0; margin: 0 0 16px; }
    .cb-workloads li { display: flex; justify-content: space-between; padding: 6px 0; border-bottom: 1px solid var(--color-border); font-size: 0.9rem; }
    .cb-wl-label { font-weight: 500; }
    .cb-wl-status { color: var(--color-text-dim); font-family: 'Cascadia Code', monospace; font-size: 0.82rem; }
    .cb-panel button { background: var(--color-accent); color: #fff; border: none; padding: 8px 18px; border-radius: 6px; cursor: pointer; font-size: 0.9rem; margin-right: 8px; }
    .cb-panel button:disabled { opacity: 0.4; cursor: default; }
    .cb-panel button:hover:not([disabled]) { opacity: 0.85; }
    .cb-results-table { width: 100%; font-size: 0.85rem; border-collapse: collapse; margin-bottom: 12px; }
    .cb-results-table th, .cb-results-table td { padding: 6px 10px; border-bottom: 1px solid var(--color-border); }
    .cb-results-table th { background: var(--color-input); text-align: left; }
    .cb-results-table td.num, .cb-results-table th.num { text-align: right; font-family: 'Cascadia Code', monospace; }
    .cb-export { margin-bottom: 12px; }
    .cb-json pre { background: var(--color-input); padding: 12px; border-radius: 6px; font-size: 0.78rem; overflow-x: auto; max-height: 400px; }
    .cb-err { color: var(--color-error, #e05050); font-size: 0.75rem; }
    .cb-stress-controls { display: flex; flex-wrap: wrap; gap: 12px 18px; align-items: flex-end; margin-bottom: 4px; }
    .cb-stress-controls label { display: flex; flex-direction: column; gap: 4px; font-size: 0.78rem; color: var(--color-text-dim); }
    .cb-stress-controls select, .cb-stress-controls input { background: var(--color-input); color: var(--color-text); border: 1px solid var(--color-border); border-radius: 5px; padding: 5px 8px; font-size: 0.85rem; }
    .cb-stress-controls input[type="text"] { min-width: 260px; font-family: 'Cascadia Code', monospace; font-size: 0.78rem; }
    .cb-stress-row { background: var(--color-input); padding: 10px 12px; border-radius: 6px; font-size: 0.75rem; overflow-x: auto; white-space: pre; margin: 0 0 10px; }
  `,document.head.appendChild(e)}_t();
