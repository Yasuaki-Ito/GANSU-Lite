import"./styles-B8WF9G-F.js";import{p as g}from"./parseXYZ-BE3rbSP_.js";import{g as k,h as $,z as _,o as R,f as L,n as I}from"./ri-Cuj4t-H2.js";import{c as j}from"./gradient-gHuNZ-_b.js";import{c as N}from"./hessian-Co-zz1lw.js";import"./properties-DST8hwXs.js";const y=4,A=1,v=1.8897259886,O=`3
H2O standard
O   0.000000   0.000000   0.117790
H   0.000000   0.755450  -0.471160
H   0.000000  -0.755450  -0.471160`,z=`3
H2O displaced
O   0.000000   0.000000   0.000000
H   0.000000   0.850000   0.620000
H   0.000000  -0.850000   0.620000`,G=`3
CO2 linear
C   0.000000   0.000000   0.000000
O   0.000000   0.000000   1.160000
O   0.000000   0.000000  -1.160000`,x=new Map;async function w(t){if(x.has(t))return x.get(t);const a=`/GANSU-Lite/basis/${t}.gbs`,n=await fetch(a);if(!n.ok)throw new Error(`Failed to load basis ${t}: ${n.status}`);const e=L.fromGBS(await n.text());return x.set(t,e),e}async function T(t){const a=await w("6-31g(d,p)"),n=g(O),e=new k(n,a,0);await $(e,"RHF").solve({eriBackend:t})}async function D(t){const a=await w("sto-3g"),n=g(z);let e=new Float64Array(n.length*3);for(let r=0;r<n.length;r++)e[3*r]=n[r].coordinate.x,e[3*r+1]=n[r].coordinate.y,e[3*r+2]=n[r].coordinate.z;const o=.4,c=25;for(let r=0;r<c;r++){const i=n.map((s,p)=>({...s,coordinate:{x:e[3*p],y:e[3*p+1],z:e[3*p+2]}})),d=new k(i,a,0),l=$(d,"RHF");await l.solve({eriBackend:t});const m=j(d.primitiveShells,d.atoms,d.cgtoNormalizationFactors,d.numBasis,d.numAlphaSpins,l.density,l.coefficients,l.orbitalEnergies).total;let u=0;for(let s=0;s<m.length;s++)Math.abs(m[s])>u&&(u=Math.abs(m[s]));if(u<5e-4)break;for(let s=0;s<m.length;s++)e[s]-=o*m[s]}}async function U(t){const a=await w("sto-3g"),n=g(G),e=new Float64Array(n.length*3);for(let o=0;o<n.length;o++)e[3*o]=n[o].coordinate.x*v,e[3*o+1]=n[o].coordinate.y*v,e[3*o+2]=n[o].coordinate.z*v;await N(n.map(o=>o.atomicNumber),e,a,0,5e-4,void 0,void 0)}async function W(t){const a=await w("6-31g(d,p)"),n=g(O),e=new k(n,a,0);await $(e,"RHF",{functional:"B3LYP",gridLevel:"medium"}).solve({eriBackend:t})}const M=[{id:"h2o_hf_631gdp",label:"H2O HF/6-31G(d,p)",run:T},{id:"h2o_geomopt_sto3g",label:"H2O geom opt HF/STO-3G",run:D},{id:"co2_hess_sto3g",label:"CO2 Hessian HF/STO-3G",run:U},{id:"h2o_b3lyp_631gdp",label:"H2O B3LYP/6-31G(d,p)",run:W}];function J(t){const a=[...t].sort((e,o)=>e-o),n=a.length;return n===0?NaN:n%2===1?a[(n-1)/2]:(a[n/2-1]+a[n/2])/2}function F(){const t=navigator;return{ua:t.userAgent,platform:t.platform??"unknown",deviceMemoryGB:t.deviceMemory??null,hardwareConcurrency:t.hardwareConcurrency??0,screen:`${screen.width}×${screen.height}`,language:t.language||"unknown",timestamp:new Date().toISOString()}}function C(){try{const t=performance.getEntriesByType("navigation")[0];return t?Math.round(t.loadEventEnd-t.fetchStart):null}catch{return null}}async function E(){return await _(),R()?I()==="wasm-simd"?{backend:"wasm",label:"WASM-SIMD"}:{backend:"wasm",label:"WASM"}:{backend:"js",label:"JS"}}const Y=document.getElementById("app");let b=null;function P(){const t=F(),a=C();Y.innerHTML=`
    <div class="cb-page">
      <header class="cb-header">
        <h1>GANSU Lite — Cross-device Benchmark</h1>
        <p class="cb-subtitle">Reproducible end-to-end timing for the paper. Median of 3 timed runs (+ 1 warmup discarded).</p>
      </header>

      <section class="cb-panel">
        <h2>Device & Backend</h2>
        <table class="cb-info">
          <tr><th>User-Agent</th><td>${t.ua}</td></tr>
          <tr><th>Platform</th><td>${t.platform}</td></tr>
          <tr><th>Device memory</th><td>${t.deviceMemoryGB?`${t.deviceMemoryGB} GB`:"unknown"}</td></tr>
          <tr><th>Cores (navigator)</th><td>${t.hardwareConcurrency}</td></tr>
          <tr><th>Screen</th><td>${t.screen}</td></tr>
          <tr><th>Initial load</th><td>${a!=null?`${a} ms`:"unavailable"}</td></tr>
          <tr><th>Backend</th><td id="cb-backend">detecting…</td></tr>
        </table>
      </section>

      <section class="cb-panel">
        <h2>Workloads</h2>
        <p class="cb-note">All fixed: SCF tol = 1e-7, DIIS, default initial guess. Single Web Worker thread (no parallelism — single-core perf + memory bandwidth).</p>
        <ul class="cb-workloads">
          ${M.map(n=>`<li id="cb-row-${n.id}"><span class="cb-wl-label">${n.label}</span><span class="cb-wl-status">pending</span></li>`).join("")}
        </ul>
        <button id="cb-run">Run all workloads</button>
        <button id="cb-stop" disabled>Stop</button>
      </section>

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
  `,K(),document.getElementById("cb-run").addEventListener("click",X),document.getElementById("cb-stop").addEventListener("click",()=>{h=!0}),document.getElementById("cb-copy-json")?.addEventListener("click",q),document.getElementById("cb-download-json")?.addEventListener("click",()=>S("benchmark.json",JSON.stringify(b,null,2),"application/json")),document.getElementById("cb-download-csv")?.addEventListener("click",()=>S("benchmark.csv",V(b),"text/csv")),E().then(({label:n})=>{const e=document.getElementById("cb-backend");e&&(e.textContent=n)})}let h=!1;function f(t,a){const n=document.getElementById(`cb-row-${t}`);if(n){const e=n.querySelector(".cb-wl-status");e&&(e.textContent=a)}}async function X(){h=!1;const t=document.getElementById("cb-run"),a=document.getElementById("cb-stop");t.disabled=!0,a.disabled=!1;const{backend:n,label:e}=await E(),o=F(),c=C(),r=[];for(const i of M){if(h){f(i.id,"stopped");continue}const d=[];let l;try{for(let s=0;s<y&&!h;s++){f(i.id,`run ${s+1}/${y}…`);const p=performance.now();await i.run(n);const B=performance.now()-p;d.push(B),f(i.id,`run ${s+1}/${y}: ${B.toFixed(0)} ms`),await new Promise(H=>setTimeout(H,50))}}catch(s){l=s instanceof Error?s.message:String(s),f(i.id,`ERROR: ${l}`)}const m=d.slice(A),u=m.length>0?J(m):NaN;r.push({id:i.id,label:i.label,warmupMs:d[0]??NaN,runs:m,medianMs:u,...l?{error:l}:{}}),!l&&m.length>0&&f(i.id,`median ${u.toFixed(0)} ms`)}b={device:o,backend:e,initialLoadMs:c,workloads:r,toolUrl:location.href},Z(b),t.disabled=!1,a.disabled=!0}function Z(t){const a=document.getElementById("cb-results");a.style.display="";const n=document.getElementById("cb-results-body");n.innerHTML="";for(const e of t.workloads){const o=document.createElement("tr");o.innerHTML=`
      <td>${e.label}${e.error?` <span class="cb-err">(${e.error})</span>`:""}</td>
      <td class="num">${isFinite(e.warmupMs)?e.warmupMs.toFixed(0):"—"}</td>
      <td class="num">${e.runs[0]!=null?e.runs[0].toFixed(0):"—"}</td>
      <td class="num">${e.runs[1]!=null?e.runs[1].toFixed(0):"—"}</td>
      <td class="num">${e.runs[2]!=null?e.runs[2].toFixed(0):"—"}</td>
      <td class="num"><b>${isFinite(e.medianMs)?e.medianMs.toFixed(0):"—"}</b></td>
    `,n.appendChild(o)}document.getElementById("cb-json-pre").textContent=JSON.stringify(t,null,2)}function q(){b&&navigator.clipboard.writeText(JSON.stringify(b,null,2)).then(()=>alert("JSON copied to clipboard")).catch(()=>alert("Copy failed — open browser console"))}function S(t,a,n){const e=new Blob([a],{type:n}),o=URL.createObjectURL(e),c=document.createElement("a");c.href=o,c.download=t,c.click(),setTimeout(()=>URL.revokeObjectURL(o),1e3)}function V(t){if(!t)return"";const a=["workload_id","workload_label","backend","warmup_ms","run1_ms","run2_ms","run3_ms","median_ms","error"],n=t.workloads.map(o=>[o.id,o.label,t.backend,isFinite(o.warmupMs)?o.warmupMs.toFixed(2):"",o.runs[0]?.toFixed(2)??"",o.runs[1]?.toFixed(2)??"",o.runs[2]?.toFixed(2)??"",isFinite(o.medianMs)?o.medianMs.toFixed(2):"",o.error??""].map(c=>`"${String(c).replace(/"/g,'""')}"`).join(","));return[`# device.ua,"${t.device.ua.replace(/"/g,'""')}"`,`# device.platform,"${t.device.platform}"`,`# device.deviceMemoryGB,"${t.device.deviceMemoryGB??""}"`,`# device.hardwareConcurrency,"${t.device.hardwareConcurrency}"`,`# device.screen,"${t.device.screen}"`,`# backend,"${t.backend}"`,`# initialLoadMs,"${t.initialLoadMs??""}"`,`# timestamp,"${t.device.timestamp}"`].join(`
`)+`
`+a.join(",")+`
`+n.join(`
`)+`
`}function K(){if(document.getElementById("cb-styles"))return;const t=document.createElement("style");t.id="cb-styles",t.textContent=`
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
  `,document.head.appendChild(t)}P();
