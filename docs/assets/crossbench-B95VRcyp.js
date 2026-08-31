import"./styles-B8WF9G-F.js";import{p as E}from"./parseXYZ-ByfvMHmk.js";import{g as j,h as O,z as et,o as nt,f as st,n as ot,C as at,F as rt,M as it}from"./ri-DT3jfd0v.js";import{c as ct}from"./gradient-Bmg7fDyY.js";import{c as lt}from"./hessian-DpHI2yIv.js";import"./properties-Odo1T7fO.js";const dt="modulepreload",mt=function(t){return"/GANSU-Lite/"+t},W={},ut=function(r,s,n){let a=Promise.resolve();if(s&&s.length>0){let m=function(d){return Promise.all(d.map(c=>Promise.resolve(c).then(p=>({status:"fulfilled",value:p}),p=>({status:"rejected",reason:p}))))};document.getElementsByTagName("link");const o=document.querySelector("meta[property=csp-nonce]"),l=o?.nonce||o?.getAttribute("nonce");a=m(s.map(d=>{if(d=mt(d),d in W)return;W[d]=!0;const c=d.endsWith(".css"),p=c?'[rel="stylesheet"]':"";if(document.querySelector(`link[href="${d}"]${p}`))return;const i=document.createElement("link");if(i.rel=c?"stylesheet":dt,c||(i.as="script"),i.crossOrigin="",i.href=d,l&&i.setAttribute("nonce",l),document.head.appendChild(i),c)return new Promise((u,b)=>{i.addEventListener("load",u),i.addEventListener("error",()=>b(new Error(`Unable to preload CSS for ${d}`)))})}))}function e(o){const l=new Event("vite:preloadError",{cancelable:!0});if(l.payload=o,window.dispatchEvent(l),!l.defaultPrevented)throw o}return a.then(o=>{for(const l of o||[])l.status==="rejected"&&e(l.reason);return r().catch(e)})},I=4,bt=1,A=1.8897259886,U=`3
H2O standard
O   0.000000   0.000000   0.117790
H   0.000000   0.755450  -0.471160
H   0.000000  -0.755450  -0.471160`,pt=`3
H2O displaced
O   0.000000   0.000000   0.000000
H   0.000000   0.850000   0.620000
H   0.000000  -0.850000   0.620000`,ft=`3
CO2 linear
C   0.000000   0.000000   0.000000
O   0.000000   0.000000   1.160000
O   0.000000   0.000000  -1.160000`,L=new Map;async function $(t){if(L.has(t))return L.get(t);const r=`/GANSU-Lite/basis/${t}.gbs`,s=await fetch(r);if(!s.ok)throw new Error(`Failed to load basis ${t}: ${s.status}`);const n=st.fromGBS(await s.text());return L.set(t,n),n}async function ht(t){const r=await $("6-31g(d,p)"),s=E(U),n=new j(s,r,0);await O(n,"RHF").solve({eriBackend:t})}async function wt(t){const r=await $("sto-3g"),s=E(pt);let n=new Float64Array(s.length*3);for(let o=0;o<s.length;o++)n[3*o]=s[o].coordinate.x,n[3*o+1]=s[o].coordinate.y,n[3*o+2]=s[o].coordinate.z;const a=.4,e=25;for(let o=0;o<e;o++){const l=s.map((i,u)=>({...i,coordinate:{x:n[3*u],y:n[3*u+1],z:n[3*u+2]}})),m=new j(l,r,0),d=O(m,"RHF");await d.solve({eriBackend:t});const c=ct(m.primitiveShells,m.atoms,m.cgtoNormalizationFactors,m.numBasis,m.numAlphaSpins,d.density,d.coefficients,d.orbitalEnergies).total;let p=0;for(let i=0;i<c.length;i++)Math.abs(c[i])>p&&(p=Math.abs(c[i]));if(p<5e-4)break;for(let i=0;i<c.length;i++)n[i]-=a*c[i]}}async function gt(t){const r=await $("sto-3g"),s=E(ft),n=new Float64Array(s.length*3);for(let a=0;a<s.length;a++)n[3*a]=s[a].coordinate.x*A,n[3*a+1]=s[a].coordinate.y*A,n[3*a+2]=s[a].coordinate.z*A;await lt(s.map(a=>a.atomicNumber),n,r,0,5e-4,void 0,void 0)}async function yt(t){const r=await $("6-31g(d,p)"),s=E(U),n=new j(s,r,0);await O(n,"RHF",{functional:"B3LYP",gridLevel:"medium"}).solve({eriBackend:t})}const J=[{id:"h2o_hf_631gdp",label:"H2O HF/6-31G(d,p)",run:ht},{id:"h2o_geomopt_sto3g",label:"H2O geom opt HF/STO-3G",run:wt},{id:"co2_hess_sto3g",label:"CO2 Hessian HF/STO-3G",run:gt},{id:"h2o_b3lyp_631gdp",label:"H2O B3LYP/6-31G(d,p)",run:yt}];function T(t){const r=[...t].sort((n,a)=>n-a),s=r.length;return s===0?NaN:s%2===1?r[(s-1)/2]:(r[s/2-1]+r[s/2])/2}function P(){const t=navigator;return{ua:t.userAgent,platform:t.platform??"unknown",deviceMemoryGB:t.deviceMemory??null,hardwareConcurrency:t.hardwareConcurrency??0,screen:`${screen.width}×${screen.height}`,language:t.language||"unknown",timestamp:new Date().toISOString()}}function Y(){try{const t=performance.getEntriesByType("navigation")[0];return t?Math.round(t.loadEventEnd-t.fetchStart):null}catch{return null}}async function V(t){const r=await fetch(t);if(!r.ok)throw new Error(`Cannot fetch ${t}`);const s=await r.arrayBuffer(),n=await WebAssembly.compile(s);let a;const e={"./wasm_eri_bg.js":{__wbindgen_init_externref_table:()=>{const i=a.exports.__wbindgen_externrefs,u=i.grow(4);i.set(0,void 0),i.set(u+0,void 0),i.set(u+1,null),i.set(u+2,!0),i.set(u+3,!1)}}};a=await WebAssembly.instantiate(n,e);const o=a.exports;o.__wbindgen_start&&o.__wbindgen_start();const l=o.memory,m=o.__wbindgen_malloc,d=o.__wbindgen_free,c=i=>{const u=m(i.length*8,8)>>>0;return new Float64Array(l.buffer).set(i,u/8),[u,i.length]},p=i=>{const[u,b]=i,f=new Float64Array(l.buffer).slice(u/8,u/8+b);return d(u,b*8,8),f};return{computeERIs:(i,u,b,f)=>{const[w,h]=c(i),[g,M]=c(u);return p(o.compute_eris_wasm(w,h,g,M,b,f))},computeFockRhf:(i,u,b,f)=>{const[w,h]=c(i),[g,M]=c(u),[x,B]=c(b);return p(o.compute_fock_rhf(w,h,g,M,x,B,f))}}}let S,_;async function Mt(){if(S!==void 0)return S;try{S=await V("/GANSU-Lite/wasm/wasm_eri_bg.wasm")}catch{S=null}return S}async function vt(){if(_!==void 0)return _;try{_=await V("/GANSU-Lite/wasm/wasm_eri_simd_bg.wasm")}catch{_=null}return _}async function X(){return await et("/GANSU-Lite/"),nt()?ot()==="wasm-simd"?{backend:"wasm",label:"WASM-SIMD"}:{backend:"wasm",label:"WASM"}:{backend:"js",label:"JS"}}function G(t,r,s,n){let a=t,e=r;if(a<e){const c=a;a=e,e=c}let o=s,l=n;if(o<l){const c=o;o=l,l=c}let m=a*(a+1)/2+e,d=o*(o+1)/2+l;if(m<d){const c=m;m=d,d=c}return m*(m+1)/2+d}function xt(t,r,s,n){const a=new Float64Array(n*n);for(let e=0;e<n;e++)for(let o=e;o<n;o++){let l=0;for(let d=0;d<n;d++)for(let c=0;c<n;c++){const p=r[d*n+c];Math.abs(p)<1e-10||(l+=p*(t[G(e,o,d,c)]-.5*t[G(e,d,o,c)]))}const m=s[e*n+o]+l;a[e*n+o]=m,a[o*n+e]=m}return a}async function St(t,r,s,n){const a=t.numBasis,e=new rt(r,a),o=O(t,"RHF");o.computeNuclearRepulsionEnergy(),o.computeCoreHamiltonianMatrix(),o.computeTransformMatrix(),o.eri=e,o.guessInitialFockMatrix();const l=o.coreHamiltonianMatrix.data;let m=0,d=0,c=0;const p=performance.now();for(let i=0;i<200;i++){o.computeCoefficientMatrix(),o.computeDensityMatrix();const u=o.densityMatrix.data,b=await n(r,u,l,a);if(o.fockMatrix=new it(a,a,b),c=o.computeEnergy(),i>0&&Math.abs(c-m)<1e-8){d=i+1;break}o.updateFockMatrix(),m=c}return{eriMs:s,scfMs:performance.now()-p,iters:d,energy:c}}const H=[{molecule:"H2O",xyzFile:"H2O.xyz",basis:"6-31g(d,p)"},{molecule:"CO2",xyzFile:"CO2.xyz",basis:"6-31g(d,p)"},{molecule:"Benzene",xyzFile:"Benzene.xyz",basis:"6-31g(d,p)"}],D=3;async function _t(t){const r=await $(t.basis),s=`/GANSU-Lite/xyz/${t.xyzFile}`,n=await fetch(s);if(!n.ok)throw new Error(`Cannot load ${s}`);const a=await n.text(),e=E(a);return{mol:new j(e,r,0),xyz:a}}async function kt(t,r,s,n){const{mol:a}=await _t(t),e=a.numBasis,o=a.primitiveShells,l=a.cgtoNormalizationFactors,m=(await ut(async()=>{const{packShells:b}=await import("./ri-DT3jfd0v.js").then(f=>f.W);return{packShells:b}},[])).packShells(o),d=new Float64Array(l);async function c(b,f,w){const h=[],g=[];let M=0,x;try{for(let y=0;y<D;y++){n(`  ${b} run ${y+1}/${D}…`);const Z=performance.now(),K=await f(),Q=performance.now()-Z,N=await St(a,K,Q,w);h.push(N.eriMs),g.push(N.scfMs),y===0&&(M=N.iters),await new Promise(tt=>setTimeout(tt,30))}}catch(y){x=y instanceof Error?y.message:String(y)}const B=h.length?T(h):NaN,z=g.length?T(g):NaN;return{eriMedianMs:B,scfMedianMs:z,totalMedianMs:B+z,iters:M,eriRuns:h,scfRuns:g,...x?{error:x}:{}}}const p=await c("JS",()=>at(o,l,e,1e-10,void 0,"js"),(b,f,w,h)=>xt(b,f,w,h)),i=r?await c("WASM baseline",()=>r.computeERIs(m,d,e,1e-10),(b,f,w,h)=>r.computeFockRhf(b,f,w,h)):{eriMedianMs:NaN,scfMedianMs:NaN,totalMedianMs:NaN,iters:0,eriRuns:[],scfRuns:[],error:"WASM baseline not available"},u=s?await c("WASM+SIMD",()=>s.computeERIs(m,d,e,1e-10),(b,f,w,h)=>s.computeFockRhf(b,f,w,h)):{eriMedianMs:NaN,scfMedianMs:NaN,totalMedianMs:NaN,iters:0,eriRuns:[],scfRuns:[],error:"SIMD not available"};return{molecule:t.molecule,basis:t.basis,nbasis:e,js:p,wasmBase:i,wasmSimd:u}}let C=null;async function Et(){const t=document.getElementById("cb-run-compare"),r=document.getElementById("cb-stop-compare");t&&(t.disabled=!0),r&&(r.disabled=!1);const s=o=>{const l=document.getElementById("cb-compare-status");l&&(l.textContent=o)};s("Loading WASM binaries…");const n=await Mt(),a=await vt(),e=[];for(let o=0;o<H.length;o++){const l=H[o];s(`[${o+1}/${H.length}] ${l.molecule}/${l.basis}…`);try{e.push(await kt(l,n,a,s)),$t(e)}catch(m){const d=m instanceof Error?m.message:String(m);s(`  ERROR: ${d}`)}}C=e,s("Done."),t&&(t.disabled=!1),r&&(r.disabled=!0)}function $t(t){const r=document.getElementById("cb-compare-section");if(!r)return;r.style.display="";const s=document.getElementById("cb-compare-body");if(!s)return;s.innerHTML="";const n=e=>isFinite(e)?e.toFixed(0):"—",a=(e,o)=>isFinite(e)&&isFinite(o)&&o>0?`${(e/o).toFixed(1)}×`:"—";for(const e of t){const o=document.createElement("tr");o.innerHTML=`
      <td>${e.molecule}</td><td>${e.basis}</td><td class="num">${e.nbasis}</td><td class="num">${e.js.iters||"—"}</td>
      <td class="num">${n(e.js.eriMedianMs)}</td>
      <td class="num">${n(e.js.scfMedianMs)}</td>
      <td class="num"><b>${n(e.js.totalMedianMs)}</b></td>
      <td class="num">${n(e.wasmBase.eriMedianMs)}</td>
      <td class="num">${n(e.wasmBase.scfMedianMs)}</td>
      <td class="num"><b>${n(e.wasmBase.totalMedianMs)}</b></td>
      <td class="num">${a(e.js.totalMedianMs,e.wasmBase.totalMedianMs)}</td>
      <td class="num">${n(e.wasmSimd.eriMedianMs)}</td>
      <td class="num">${n(e.wasmSimd.scfMedianMs)}</td>
      <td class="num"><b>${n(e.wasmSimd.totalMedianMs)}</b></td>
      <td class="num">${a(e.js.totalMedianMs,e.wasmSimd.totalMedianMs)}</td>
    `,s.appendChild(o)}}const Bt=document.getElementById("app");let v=null;function Ft(){const t=P(),r=Y();Bt.innerHTML=`
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
          <tr><th>Initial load</th><td>${r!=null?`${r} ms`:"unavailable"}</td></tr>
          <tr><th>Backend</th><td id="cb-backend">detecting…</td></tr>
        </table>
      </section>

      <section class="cb-panel">
        <h2>Workloads</h2>
        <p class="cb-note">All fixed: SCF tol = 1e-7, DIIS, default initial guess. Single Web Worker thread (no parallelism — single-core perf + memory bandwidth).</p>
        <ul class="cb-workloads">
          ${J.map(s=>`<li id="cb-row-${s.id}"><span class="cb-wl-label">${s.label}</span><span class="cb-wl-status">pending</span></li>`).join("")}
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
  `,It(),document.getElementById("cb-run").addEventListener("click",Ct),document.getElementById("cb-stop").addEventListener("click",()=>{R=!0}),document.getElementById("cb-run-compare")?.addEventListener("click",Et),document.getElementById("cb-copy-compare-json")?.addEventListener("click",()=>q(JSON.stringify(C,null,2))),document.getElementById("cb-download-compare-json")?.addEventListener("click",()=>F("backend-comparison.json",JSON.stringify(C,null,2),"application/json")),document.getElementById("cb-download-compare-csv")?.addEventListener("click",()=>F("backend-comparison.csv",Ot(C),"text/csv")),document.getElementById("cb-copy-json")?.addEventListener("click",jt),document.getElementById("cb-download-json")?.addEventListener("click",()=>F("benchmark.json",JSON.stringify(v,null,2),"application/json")),document.getElementById("cb-download-csv")?.addEventListener("click",()=>F("benchmark.csv",Nt(v),"text/csv")),X().then(({label:s})=>{const n=document.getElementById("cb-backend");n&&(n.textContent=s)})}let R=!1;function k(t,r){const s=document.getElementById(`cb-row-${t}`);if(s){const n=s.querySelector(".cb-wl-status");n&&(n.textContent=r)}}async function Ct(){R=!1;const t=document.getElementById("cb-run"),r=document.getElementById("cb-stop");t.disabled=!0,r.disabled=!1;const{backend:s,label:n}=await X(),a=P(),e=Y(),o=[];for(const l of J){if(R){k(l.id,"stopped");continue}const m=[];let d;try{for(let i=0;i<I&&!R;i++){k(l.id,`run ${i+1}/${I}…`);const u=performance.now();await l.run(s);const b=performance.now()-u;m.push(b),k(l.id,`run ${i+1}/${I}: ${b.toFixed(0)} ms`),await new Promise(f=>setTimeout(f,50))}}catch(i){d=i instanceof Error?i.message:String(i),k(l.id,`ERROR: ${d}`)}const c=m.slice(bt),p=c.length>0?T(c):NaN;o.push({id:l.id,label:l.label,warmupMs:m[0]??NaN,runs:c,medianMs:p,...d?{error:d}:{}}),!d&&c.length>0&&k(l.id,`median ${p.toFixed(0)} ms`)}v={device:a,backend:n,initialLoadMs:e,workloads:o,toolUrl:location.href},Rt(v),t.disabled=!1,r.disabled=!0}function Rt(t){const r=document.getElementById("cb-results");r.style.display="";const s=document.getElementById("cb-results-body");s.innerHTML="";for(const n of t.workloads){const a=document.createElement("tr");a.innerHTML=`
      <td>${n.label}${n.error?` <span class="cb-err">(${n.error})</span>`:""}</td>
      <td class="num">${isFinite(n.warmupMs)?n.warmupMs.toFixed(0):"—"}</td>
      <td class="num">${n.runs[0]!=null?n.runs[0].toFixed(0):"—"}</td>
      <td class="num">${n.runs[1]!=null?n.runs[1].toFixed(0):"—"}</td>
      <td class="num">${n.runs[2]!=null?n.runs[2].toFixed(0):"—"}</td>
      <td class="num"><b>${isFinite(n.medianMs)?n.medianMs.toFixed(0):"—"}</b></td>
    `,s.appendChild(a)}document.getElementById("cb-json-pre").textContent=JSON.stringify(t,null,2)}function jt(){v&&q(JSON.stringify(v,null,2))}function q(t){navigator.clipboard.writeText(t).then(()=>alert("Copied to clipboard")).catch(()=>alert("Copy failed — open browser console"))}function Ot(t){if(!t||t.length===0)return"";const r=["molecule","basis","nbasis","iters","js_eri_ms","js_scf_ms","js_total_ms","wasm_eri_ms","wasm_scf_ms","wasm_total_ms","speedup_wasm_vs_js","simd_eri_ms","simd_scf_ms","simd_total_ms","speedup_simd_vs_js"],s=e=>isFinite(e)?e.toFixed(2):"",n=(e,o)=>isFinite(e)&&isFinite(o)&&o>0?(e/o).toFixed(2):"",a=t.map(e=>[e.molecule,e.basis,String(e.nbasis),String(e.js.iters),s(e.js.eriMedianMs),s(e.js.scfMedianMs),s(e.js.totalMedianMs),s(e.wasmBase.eriMedianMs),s(e.wasmBase.scfMedianMs),s(e.wasmBase.totalMedianMs),n(e.js.totalMedianMs,e.wasmBase.totalMedianMs),s(e.wasmSimd.eriMedianMs),s(e.wasmSimd.scfMedianMs),s(e.wasmSimd.totalMedianMs),n(e.js.totalMedianMs,e.wasmSimd.totalMedianMs)].map(o=>`"${o}"`).join(","));return r.join(",")+`
`+a.join(`
`)+`
`}function F(t,r,s){const n=new Blob([r],{type:s}),a=URL.createObjectURL(n),e=document.createElement("a");e.href=a,e.download=t,e.click(),setTimeout(()=>URL.revokeObjectURL(a),1e3)}function Nt(t){if(!t)return"";const r=["workload_id","workload_label","backend","warmup_ms","run1_ms","run2_ms","run3_ms","median_ms","error"],s=t.workloads.map(a=>[a.id,a.label,t.backend,isFinite(a.warmupMs)?a.warmupMs.toFixed(2):"",a.runs[0]?.toFixed(2)??"",a.runs[1]?.toFixed(2)??"",a.runs[2]?.toFixed(2)??"",isFinite(a.medianMs)?a.medianMs.toFixed(2):"",a.error??""].map(e=>`"${String(e).replace(/"/g,'""')}"`).join(","));return[`# device.ua,"${t.device.ua.replace(/"/g,'""')}"`,`# device.platform,"${t.device.platform}"`,`# device.deviceMemoryGB,"${t.device.deviceMemoryGB??""}"`,`# device.hardwareConcurrency,"${t.device.hardwareConcurrency}"`,`# device.screen,"${t.device.screen}"`,`# backend,"${t.backend}"`,`# initialLoadMs,"${t.initialLoadMs??""}"`,`# timestamp,"${t.device.timestamp}"`].join(`
`)+`
`+r.join(",")+`
`+s.join(`
`)+`
`}function It(){if(document.getElementById("cb-styles"))return;const t=document.createElement("style");t.id="cb-styles",t.textContent=`
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
  `,document.head.appendChild(t)}Ft();
