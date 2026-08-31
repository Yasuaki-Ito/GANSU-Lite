import"./styles-B8WF9G-F.js";import{p as L}from"./parseXYZ-ByfvMHmk.js";import{g as G,h as U,z as Ce,o as Re,f as Le,n as Ie,C as Ae,F as je,M as Ne}from"./ri-DT3jfd0v.js";import{c as Oe}from"./gradient-Bmg7fDyY.js";import{c as Te}from"./hessian-DpHI2yIv.js";import"./properties-Odo1T7fO.js";const He="modulepreload",ze=function(e){return"/GANSU-Lite/"+e},ae={},Pe=function(s,n,t){let o=Promise.resolve();if(n&&n.length>0){let u=function(l){return Promise.all(l.map(c=>Promise.resolve(c).then(b=>({status:"fulfilled",value:b}),b=>({status:"rejected",reason:b}))))};document.getElementsByTagName("link");const a=document.querySelector("meta[property=csp-nonce]"),i=a?.nonce||a?.getAttribute("nonce");o=u(n.map(l=>{if(l=ze(l),l in ae)return;ae[l]=!0;const c=l.endsWith(".css"),b=c?'[rel="stylesheet"]':"";if(document.querySelector(`link[href="${l}"]${b}`))return;const d=document.createElement("link");if(d.rel=c?"stylesheet":He,c||(d.as="script"),d.crossOrigin="",d.href=l,i&&d.setAttribute("nonce",i),document.head.appendChild(d),c)return new Promise((m,y)=>{d.addEventListener("load",m),d.addEventListener("error",()=>y(new Error(`Unable to preload CSS for ${l}`)))})}))}function r(a){const i=new Event("vite:preloadError",{cancelable:!0});if(i.payload=a,window.dispatchEvent(i),!i.defaultPrevented)throw a}return o.then(a=>{for(const i of a||[])i.status==="rejected"&&r(i.reason);return s().catch(r)})};function Q(){const e=navigator;return{ua:e.userAgent,platform:e.platform??"unknown",deviceMemoryGB:e.deviceMemory??null,hardwareConcurrency:e.hardwareConcurrency??0,screen:`${screen.width}×${screen.height}`,language:e.language||"unknown",timestamp:new Date().toISOString()}}function pe(){try{const e=performance.getEntriesByType("navigation")[0];return e?Math.round(e.loadEventEnd-e.fetchStart):null}catch{return null}}const j=.9572,N=104.52*Math.PI/180,O=3,We=[["O",0,0,0],["H",j*Math.sin(N/2),0,j*Math.cos(N/2)],["H",-j*Math.sin(N/2),0,j*Math.cos(N/2)]];function Ge(e){const s=Math.ceil(Math.cbrt(e))+1,n=[];for(let t=-s;t<=s;t++)for(let o=-s;o<=s;o++)for(let r=-s;r<=s;r++)n.push([t,o,r]);return n.sort((t,o)=>{const r=t[0]*t[0]+t[1]*t[1]+t[2]*t[2],a=o[0]*o[0]+o[1]*o[1]+o[2]*o[2];return r!==a?r-a:t[2]!==o[2]?t[2]-o[2]:t[1]!==o[1]?t[1]-o[1]:t[0]-o[0]}),n.slice(0,e)}function Ue(e){const s=Math.PI*(3-Math.sqrt(5));return[e*s,e*1.2345678,e*.7654321]}function De(e,[s,n,t]){const o=(a,i)=>[a[0]*Math.cos(i)-a[1]*Math.sin(i),a[0]*Math.sin(i)+a[1]*Math.cos(i),a[2]];return o(((a,i)=>[a[0]*Math.cos(i)+a[2]*Math.sin(i),a[1],-a[0]*Math.sin(i)+a[2]*Math.cos(i)])(o([e[0],e[1],e[2]],t),n),s)}function Je(e){const s=Ge(e),n=[];for(let t=0;t<e;t++){const[o,r,a]=s[t],i=Ue(t),u=o*O,l=r*O,c=a*O;for(const[b,d,m,y]of We){const[f,g,v]=De([d,m,y],i);n.push(`${b.padEnd(2)} ${(u+f).toFixed(6).padStart(12)} ${(l+g).toFixed(6).padStart(12)} ${(c+v).toFixed(6).padStart(12)}`)}}return`${e*3}
(H2O)${e} — cubic lattice a=${O} A, golden-angle orientations
${n.join(`
`)}`}const re=1.526,k=1.09,qe=109.47*Math.PI/180;function Ye(e){const s=qe/2,n=re*Math.sin(s),t=re*Math.cos(s),o=[];for(let c=0;c<e;c++)o.push([c*n,0,c%2===0?0:t]);const r=[];for(const c of o)r.push(F("C",c));for(let c=0;c<e;c++){const[b,d,m]=o[c],y=c%2===0?1:-1,f=m+y*k*Math.cos(s),g=k*Math.sin(s);r.push(F("H",[b,d+g,f])),r.push(F("H",[b,d-g,f]))}const a=o[0],i=o[e-1],u=[-Math.sin(s),0,e>1&&o[1][2]>a[2]?-Math.cos(s):Math.cos(s)],l=[Math.sin(s),0,e>1&&o[e-2][2]>i[2]?-Math.cos(s):Math.cos(s)];return r.push(F("H",[a[0]+k*u[0],a[1],a[2]+k*u[2]])),r.push(F("H",[i[0]+k*l[0],i[1],i[2]+k*l[2]])),`${e+2*e+2}
C${e}H${2*e+2} all-anti
${r.join(`
`)}`}function F(e,[s,n,t]){return`${e.padEnd(2)} ${s.toFixed(6).padStart(12)} ${n.toFixed(6).padStart(12)} ${t.toFixed(6).padStart(12)}`}const B=[{id:"water",label:"Water cluster (H₂O)ₙ",name:e=>`(H2O)${e}`,natoms:e=>3*e,xyz:Je,ladder:[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32]},{id:"alkane",label:"n-Alkane CₙH₂ₙ₊₂",name:e=>`C${e}H${2*e+2}`,natoms:e=>3*e+2,xyz:Ye,ladder:[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20]}],ee="gansu-stress-v1",te="gansu-stress-probe-v1",be=180,V="6-31g(d,p)";function Xe(){try{const e=localStorage.getItem(ee);return e?JSON.parse(e):null}catch{return null}}function Z(e){try{localStorage.setItem(ee,JSON.stringify(e))}catch{}}function Ve(){try{localStorage.removeItem(ee)}catch{}}function P(e){return e==null?"—":e<1024**2?`${(e/1024).toFixed(0)} KB`:e<1024**3?`${(e/1024**2).toFixed(0)} MB`:`${(e/1024**3).toFixed(2)} GB`}function fe(e){return e==null?"—":(e/1e3).toFixed(2)}const se={ok:"ok",timeout:"timeout (exceeded budget)",error:"error / OOM","worker-died":"worker killed (OOM)","tab-crash":"tab crash"};function ie(e,s,n,t,o,r){return new Promise(a=>{const i=new Worker(new URL("/GANSU-Lite/assets/stressWorker-C3Q-jLBX.js",import.meta.url),{type:"module"});r(i);let u=!1;const l=performance.now();let c=l,b="no progress reported",d=null;const m=w=>{u||(u=!0,clearTimeout(g),i.terminate(),r(null),a(w))},y=w=>{const h=w*(w+1)/2;return h*(h+1)/2*8},f=()=>({nbasis:d,totalMs:null,scfMs:null,iterations:null,converged:null,energy:null,eriBytes:d==null?null:y(d),heapUsedBytes:null,heapLimitBytes:null}),g=setTimeout(()=>{const w=(performance.now()-c)/1e3;m({...f(),status:"timeout",totalMs:performance.now()-l,detail:`exceeded ${(n/1e3).toFixed(0)} s budget; silent for ${w.toFixed(0)} s after "${b}"`+(w>30?" — the worker had gone quiet, so it may have been killed rather than merely slow":" — the worker was still reporting progress, so this is genuine slowness, not a crash")})},n);i.onmessage=w=>{const h=w.data;if(c=performance.now(),h.type==="progress"){const S=/(\d+) basis functions/.exec(h.message);S&&(d=Number(S[1])),b=h.message,o(h.message,h.elapsedMs);return}if(h.type==="done"){m({status:"ok",nbasis:h.nbasis,totalMs:h.totalMs,scfMs:h.scfMs,iterations:h.iterations,converged:h.converged,energy:h.energy,eriBytes:h.eriBytes,heapUsedBytes:h.heapUsedBytes,heapLimitBytes:h.heapLimitBytes,...h.converged?{}:{detail:"SCF did not converge"}});return}h.type==="error"&&m({...f(),status:"error",nbasis:h.nbasis??d,totalMs:h.elapsedMs,detail:`${h.name}: ${h.message} (phase: ${h.phase})`})},i.onerror=w=>{m({...f(),status:"worker-died",totalMs:performance.now()-l,detail:w.message||"worker terminated by the browser"})};const v={type:"stress-run",xyzText:e,basisGBS:s,charge:0,eriBackend:"auto",baseUrl:"/GANSU-Lite/",...t?{dftConfig:{functional:"B3LYP",gridLevel:"medium"}}:{}};i.postMessage(v)})}function Ze(){const e=new Set;for(const s of B[0].ladder)e.add(s*25);return[...e].filter(s=>s<=300).sort((s,n)=>s-n)}function Ke(e){return new Promise(s=>{const n=new Worker(new URL("/GANSU-Lite/assets/stressWorker-C3Q-jLBX.js",import.meta.url),{type:"module"}),t=[],o=()=>{n.terminate(),s(t)};n.onmessage=a=>{const i=a.data;if(i.type==="probe-result"){const u={nbasis:i.nbasis,bytes:i.bytes,ok:i.ok,ms:i.ms,verdict:i.verdict,...i.detail?{detail:i.detail}:{}};t.push(u),e(u);return}i.type==="probe-done"&&o()},n.onerror=()=>o();const r={type:"memory-probe",nbasisLadder:Ze(),baseUrl:"/GANSU-Lite/"};n.postMessage(r)})}function he(e,s){const n=s.verdict==="ok",t={ok:"allocated in WASM memory",failed:"ALLOCATION FAILED",unavailable:"WASM unavailable"}[s.verdict],o=document.createElement("tr");o.innerHTML=`<td class="num">${s.nbasis||"—"}</td><td class="num">${s.bytes?P(s.bytes):"—"}</td><td class="num">${s.ms.toFixed(0)}</td><td${n?"":' class="cb-err"'}>${t}${s.detail?`<br><span class="cb-err">${ge(s.detail)}</span>`:""}</td>`,e.appendChild(o)}let x=[];function Qe(e){try{localStorage.setItem(te,JSON.stringify(e))}catch{}}function et(){try{const e=localStorage.getItem(te);return e?JSON.parse(e):[]}catch{return[]}}async function tt(){const e=document.getElementById("cb-stress-probe"),s=document.getElementById("cb-probe-body"),n=document.getElementById("cb-probe-section");e&&(e.disabled=!0),n&&(n.style.display=""),s&&(s.innerHTML=""),x=[],M("probing memory ceiling…"),await Ke(a=>{x.push(a),Qe(x),s&&he(s,a)});const t=x.filter(a=>a.verdict==="ok"),o=x.find(a=>a.verdict!=="ok"),r=document.getElementById("cb-probe-summary");if(r){const a=t.length?t[t.length-1]:null,i=x.find(u=>u.verdict==="unavailable");if(i)r.textContent=`No WASM on this device, so there is no wasm32 ceiling to measure — it runs the JS backend, whose limit is device RAM and shows up in the ladder above. ${i.detail??""}`;else if(!a)r.textContent="Could not allocate even the smallest probe array — something is wrong.";else{const u=o?{failed:`the allocator refused ${o.nbasis} (${P(o.bytes)})`,unavailable:"",ok:""}[o.verdict]:"the probe ladder topped out at 300 basis functions without hitting a ceiling";r.textContent=`WASM memory ceiling: ${o?"":"≥ "}${a.nbasis} basis functions — a ${P(a.bytes)} ERI array allocates inside the module's linear memory; ${u}.`}}M("memory probe done"),p.output&&ne(),e&&(e.disabled=!1)}function ye(e,s){let n=null;for(const t of e)if(!(t.status!=="ok"||t.totalMs==null)){if(t.totalMs/1e3>s)break;n=t}return n}function st(e,s){if(e.length<2||s.nbasis==null)return null;const n=e[e.length-2],t=e[e.length-1];if(n.nbasis==null||t.nbasis==null||n.totalMs==null||t.totalMs==null||n.nbasis>=t.nbasis||t.nbasis>=s.nbasis||n.totalMs<=0)return null;const o=Math.log(t.totalMs/n.totalMs)/Math.log(t.nbasis/n.nbasis);return!Number.isFinite(o)||o<=0?null:{seconds:t.totalMs*Math.pow(s.nbasis/t.nbasis,o)/1e3,exponent:o}}function nt(e,s){const n=e.filter(a=>a.status==="ok"),t=n.length?n[n.length-1]:null,o=e.find(a=>a.status!=="ok")??null;let r=null;if(o&&o.status==="timeout"){const a=st(n,o);a&&(r={estimatedSeconds:a.seconds,budgetSeconds:s,exponent:a.exponent,fragile:a.seconds<1.25*s})}return{borderline:r,largestCompleted:t&&t.nbasis!=null&&t.totalMs!=null?{name:t.name,natoms:t.natoms,nbasis:t.nbasis,seconds:+(t.totalMs/1e3).toFixed(1)}:null,firstFailure:o?{name:o.name,natoms:o.natoms,nbasis:o.nbasis,mode:se[o.status]}:null}}let p={series:B[0],points:[],running:!1,stopRequested:!1,activeWorker:null,output:null};function ot(){return`
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
            ${B.map(e=>`<option value="${e.id}">${e.label}</option>`).join("")}
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
          <input id="cb-stress-ladder" type="text" value="${B[0].ladder.join(", ")}">
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
        The SCF's dominant allocation is the unique-ERI array (N⁴/8 doubles), and it lives inside
        the WASM module's linear memory. The probe asks the module's own allocator for exactly that
        array, at each size, and writes across it — tens of milliseconds per size instead of the
        tens of minutes an SCF at that size would take. Synthetic substitutes do not work: a JS
        <code>Float64Array</code> and a bare <code>WebAssembly.Memory</code> both report success at
        225 basis functions, where the real calculation traps. The ceiling is wasm32's 4 GiB address
        space minus what the module already holds, so it is a property of the engine rather than of
        the machine's RAM.
      </p>
      <button id="cb-stress-probe">Probe memory ceiling</button>
      <div id="cb-probe-section" style="display:none;margin-top:10px">
        <table class="cb-results-table" style="font-size:0.8rem;max-width:520px">
          <thead><tr><th class="num">Basis fns</th><th class="num">ERI array</th><th class="num">ms</th><th>Result</th></tr></thead>
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
  `}function at(e){const s=r=>document.getElementById(r),n=s("cb-stress-series"),t=s("cb-stress-ladder");n?.addEventListener("change",()=>{const r=B.find(a=>a.id===n.value);r&&t&&(p.series=r,t.value=r.ladder.join(", "))}),s("cb-stress-run")?.addEventListener("click",()=>{lt(e)}),s("cb-stress-stop")?.addEventListener("click",()=>{p.stopRequested=!0,p.activeWorker?.terminate(),M("stopping…")}),s("cb-stress-reset")?.addEventListener("click",()=>{Ve();try{localStorage.removeItem(te)}catch{}x=[],p.points=[],p.output=null,W();const r=s("cb-stress-recovered");r&&(r.style.display="none");const a=s("cb-stress-summary");a&&(a.style.display="none"),M("saved progress cleared")}),s("cb-stress-probe")?.addEventListener("click",()=>{tt()}),s("cb-stress-threshold")?.addEventListener("input",ne),s("cb-stress-copy-row")?.addEventListener("click",()=>ce(ve())),s("cb-stress-copy-json")?.addEventListener("click",()=>ce(JSON.stringify(p.output,null,2))),s("cb-stress-download-json")?.addEventListener("click",()=>le("stress-test.json",JSON.stringify(p.output,null,2),"application/json")),s("cb-stress-download-csv")?.addEventListener("click",()=>le("stress-test.csv",ut(p.output),"text/csv")),rt();const o=et();if(o.length){x=o;const r=document.getElementById("cb-probe-section"),a=document.getElementById("cb-probe-body");if(r&&a){r.style.display="";for(const i of o)he(a,i)}}}function rt(){const e=Xe();if(!e)return;p.points=e.points??[];const s=B.find(n=>n.id===e.seriesId);if(s&&(p.series=s),e.inFlight){const n=e.inFlight;p.points.some(o=>o.n===n.n)||p.points.push({n:n.n,name:n.name,natoms:n.natoms,nbasis:null,status:"tab-crash",totalMs:null,scfMs:null,iterations:null,converged:null,energy:null,eriBytes:null,heapUsedBytes:null,heapLimitBytes:null,detail:"tab died while this point was running"});const t=document.getElementById("cb-stress-recovered");t&&(t.style.display="",t.textContent=`Recovered from a previous session: the tab died while running ${n.name}. Results up to that point are restored below and ${n.name} is recorded as a tab crash.`),Z({...e,points:p.points,inFlight:null})}p.points.length&&W()}function M(e){const s=document.getElementById("cb-stress-status");s&&(s.textContent=e)}async function it(){const e=await fetch(`/GANSU-Lite/basis/${V}.gbs`);if(!e.ok)throw new Error(`Cannot load basis ${V}: HTTP ${e.status}`);return e.text()}function ct(e,s){const n=e.split(/[,\s]+/).map(t=>parseInt(t,10)).filter(t=>Number.isFinite(t)&&t>0);return n.length?n:s}async function lt(e){if(p.running)return;const s=document.getElementById("cb-stress-run"),n=document.getElementById("cb-stress-stop"),t=Math.max(10,Number(document.getElementById("cb-stress-budget")?.value)||be),o=document.getElementById("cb-stress-theory")?.value==="b3lyp",r=ct(document.getElementById("cb-stress-ladder")?.value??"",p.series.ladder);p.running=!0,p.stopRequested=!1,p.points=[],p.output=null,s&&(s.disabled=!0),n&&(n.disabled=!1),W();let a;try{M("loading basis set…"),a=await it()}catch(c){M(`FAILED: ${c instanceof Error?c.message:String(c)}`),p.running=!1,s&&(s.disabled=!1),n&&(n.disabled=!0);return}const i=p.series,u=o?"B3LYP":"RHF",l=t*1e3;try{M("warm-up…"),await ie(i.xyz(r[0]),a,l,o,()=>{},c=>{p.activeWorker=c})}catch{}for(const c of r){if(p.stopRequested){M("stopped");break}const b=i.name(c),d=i.natoms(c);Z({seriesId:i.id,method:u,budgetSeconds:t,points:p.points,inFlight:{n:c,name:b,natoms:d,startedAt:Date.now()}}),M(`${b} — starting…`);const m=await ie(i.xyz(c),a,l,o,(f,g)=>M(`${b} — ${f} (${(g/1e3).toFixed(1)} s)`),f=>{p.activeWorker=f}),y={n:c,name:b,natoms:d,...m};if(p.points.push(y),W(),Z({seriesId:i.id,method:u,budgetSeconds:t,points:p.points,inFlight:null}),m.status!=="ok"){M(`${b} failed: ${se[m.status]} — stopping ladder`);break}M(`${b} ok — ${fe(m.totalMs)} s`),await new Promise(f=>setTimeout(f,200))}p.output={device:Q(),backend:e(),seriesId:i.id,seriesLabel:i.label,method:u,basis:V,budgetSeconds:t,points:p.points,summary:nt(p.points,t),memoryProbe:x,toolUrl:location.href},ne(),p.running=!1,p.activeWorker=null,s&&(s.disabled=!1),n&&(n.disabled=!0),p.stopRequested||M("done")}function W(){const e=document.getElementById("cb-stress-body");if(e){e.innerHTML="";for(const s of p.points){const n=document.createElement("tr"),t=s.status!=="ok";n.innerHTML=`
      <td>${s.name}</td>
      <td class="num">${s.natoms}</td>
      <td class="num">${s.nbasis??"—"}</td>
      <td class="num">${P(s.eriBytes)}</td>
      <td class="num">${s.iterations??"—"}</td>
      <td class="num">${fe(s.totalMs)}</td>
      <td${t?' class="cb-err"':""}>${se[s.status]}${s.detail?`<br><span class="cb-err">${ge(s.detail)}</span>`:""}</td>
    `,e.appendChild(n)}}}function ge(e){return e.replace(/[&<>"]/g,s=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"})[s])}function we(){const e=document.getElementById("cb-stress-threshold"),s=Number(e?.value);return Number.isFinite(s)&&s>0?s:180}function ve(){const e=p.output;if(!e)return"";const s=we(),n=ye(e.points,s),t=e.summary.largestCompleted,o=e.summary.firstFailure,r=x.filter(l=>l.verdict==="ok"),a=r.length?r[r.length-1].nbasis:null,i=x.find(l=>l.verdict==="failed");return`| ${[dt(e.device),n?n.name:"—",n?String(n.natoms):"—",n?.nbasis!=null?String(n.nbasis):"—",n&&n.totalMs!=null?(n.totalMs/1e3).toFixed(1):"—",t&&(n?.nbasis==null||t.nbasis>n.nbasis)?`${t.name} (${t.seconds}s)`:"—",o?`${o.name}: ${o.mode}`:"(no failure in ladder)",a==null?"(probe not run)":`${i?"":"≥"}${a}`].join(" | ")} |`}function dt(e){const s=e.ua;return/iPhone/.test(s)?"iPhone":/iPad/.test(s)?"iPad":/Mac OS X|Macintosh/.test(s)?"Apple-silicon laptop":/Android/.test(s)?"Android":/Windows/.test(s)?"Windows":e.platform}function ne(){const e=document.getElementById("cb-stress-summary"),s=document.getElementById("cb-stress-row");if(!e||!s||!p.output)return;e.style.display="";const n=p.output.summary,t=we(),o=`| Device | Largest within ${t}s | Atoms | Basis fns | Time (s) | Ran further (over threshold) | Ladder ended | WASM ceiling |`,r="|---|---|--:|--:|--:|---|---|--:|",a=[];if(ye(p.output.points,t)||a.push(`# Nothing completed within ${t} s — lower the threshold or read the timings directly.`),n.firstFailure||a.push("# The ladder finished without a failure — extend it to find where this device stops."),n.borderline){const l=n.borderline;a.push(`# ${n.firstFailure?.name} is extrapolated to need ~${l.estimatedSeconds.toFixed(0)} s (fitted N^${l.exponent.toFixed(1)}) against the ${l.budgetSeconds} s run budget.`),l.fragile&&a.push("# FRAGILE: that is within 25% of the budget, so run-to-run noise alone could flip it.",`# The "${t}s" column is unaffected — it is read off the recorded timings — but the`,`# ladder stopped early, so re-run at ${Math.ceil(l.estimatedSeconds*1.5/60)*60} s if you need rungs beyond this one.`)}const u=p.points.filter(l=>l.status==="ok"&&l.converged===!1);u.length&&a.push(`# SCF did not converge for: ${u.map(l=>l.name).join(", ")}`),x.length||a.push('# Memory ceiling not measured — press "Probe memory ceiling".'),s.textContent=[o,r,ve(),...a].join(`
`)}function ce(e){navigator.clipboard.writeText(e).then(()=>alert("Copied to clipboard")).catch(()=>alert("Copy failed — open browser console"))}function le(e,s,n){const t=new Blob([s],{type:n}),o=URL.createObjectURL(t),r=document.createElement("a");r.href=o,r.download=e,r.click(),setTimeout(()=>URL.revokeObjectURL(o),1e3)}function ut(e){if(!e)return"";const s=["n","system","atoms","nbasis","eri_bytes","iterations","converged","setup_plus_scf_ms","scf_ms","energy_hartree","status","detail"],n=e.points.map(o=>[o.n,o.name,o.natoms,o.nbasis??"",o.eriBytes??"",o.iterations??"",o.converged==null?"":String(o.converged),o.totalMs?.toFixed(2)??"",o.scfMs?.toFixed(2)??"",o.energy?.toFixed(8)??"",o.status,o.detail??""].map(r=>`"${String(r).replace(/"/g,'""')}"`).join(","));return[`# device.ua,"${e.device.ua.replace(/"/g,'""')}"`,`# device.platform,"${e.device.platform}"`,`# device.deviceMemoryGB,"${e.device.deviceMemoryGB??""}"`,`# device.hardwareConcurrency,"${e.device.hardwareConcurrency}"`,`# backend,"${e.backend}"`,`# series,"${e.seriesLabel}"`,`# method,"${e.method}/${e.basis}"`,`# budget_seconds,"${e.budgetSeconds}"`,`# timestamp,"${e.device.timestamp}"`].join(`
`)+`
`+s.join(",")+`
`+n.join(`
`)+`
`}const J=4,mt=1,q=1.8897259886,Me=`3
H2O standard
O   0.000000   0.000000   0.117790
H   0.000000   0.755450  -0.471160
H   0.000000  -0.755450  -0.471160`,pt=`3
H2O displaced
O   0.000000   0.000000   0.000000
H   0.000000   0.850000   0.620000
H   0.000000  -0.850000   0.620000`,bt=`3
CO2 linear
C   0.000000   0.000000   0.000000
O   0.000000   0.000000   1.160000
O   0.000000   0.000000  -1.160000`,Y=new Map;async function I(e){if(Y.has(e))return Y.get(e);const s=`/GANSU-Lite/basis/${e}.gbs`,n=await fetch(s);if(!n.ok)throw new Error(`Failed to load basis ${e}: ${n.status}`);const t=Le.fromGBS(await n.text());return Y.set(e,t),t}async function ft(e){const s=await I("6-31g(d,p)"),n=L(Me),t=new G(n,s,0);await U(t,"RHF").solve({eriBackend:e})}async function ht(e){const s=await I("sto-3g"),n=L(pt);let t=new Float64Array(n.length*3);for(let a=0;a<n.length;a++)t[3*a]=n[a].coordinate.x,t[3*a+1]=n[a].coordinate.y,t[3*a+2]=n[a].coordinate.z;const o=.4,r=25;for(let a=0;a<r;a++){const i=n.map((d,m)=>({...d,coordinate:{x:t[3*m],y:t[3*m+1],z:t[3*m+2]}})),u=new G(i,s,0),l=U(u,"RHF");await l.solve({eriBackend:e});const c=Oe(u.primitiveShells,u.atoms,u.cgtoNormalizationFactors,u.numBasis,u.numAlphaSpins,l.density,l.coefficients,l.orbitalEnergies).total;let b=0;for(let d=0;d<c.length;d++)Math.abs(c[d])>b&&(b=Math.abs(c[d]));if(b<5e-4)break;for(let d=0;d<c.length;d++)t[d]-=o*c[d]}}async function yt(e){const s=await I("sto-3g"),n=L(bt),t=new Float64Array(n.length*3);for(let o=0;o<n.length;o++)t[3*o]=n[o].coordinate.x*q,t[3*o+1]=n[o].coordinate.y*q,t[3*o+2]=n[o].coordinate.z*q;await Te(n.map(o=>o.atomicNumber),t,s,0,5e-4,void 0,void 0)}async function gt(e){const s=await I("6-31g(d,p)"),n=L(Me),t=new G(n,s,0);await U(t,"RHF",{functional:"B3LYP",gridLevel:"medium"}).solve({eriBackend:e})}const xe=[{id:"h2o_hf_631gdp",label:"H2O HF/6-31G(d,p)",run:ft},{id:"h2o_geomopt_sto3g",label:"H2O geom opt HF/STO-3G",run:ht},{id:"co2_hess_sto3g",label:"CO2 Hessian HF/STO-3G",run:yt},{id:"h2o_b3lyp_631gdp",label:"H2O B3LYP/6-31G(d,p)",run:gt}];function K(e){const s=[...e].sort((t,o)=>t-o),n=s.length;return n===0?NaN:n%2===1?s[(n-1)/2]:(s[n/2-1]+s[n/2])/2}async function Se(e){const s=await fetch(e);if(!s.ok)throw new Error(`Cannot fetch ${e}`);const n=await s.arrayBuffer(),t=await WebAssembly.compile(n);let o;const r={"./wasm_eri_bg.js":{__wbindgen_init_externref_table:()=>{const d=o.exports.__wbindgen_externrefs,m=d.grow(4);d.set(0,void 0),d.set(m+0,void 0),d.set(m+1,null),d.set(m+2,!0),d.set(m+3,!1)}}};o=await WebAssembly.instantiate(t,r);const a=o.exports;a.__wbindgen_start&&a.__wbindgen_start();const i=a.memory,u=a.__wbindgen_malloc,l=a.__wbindgen_free,c=d=>{const m=u(d.length*8,8)>>>0;return new Float64Array(i.buffer).set(d,m/8),[m,d.length]},b=d=>{const[m,y]=d,f=new Float64Array(i.buffer).slice(m/8,m/8+y);return l(m,y*8,8),f};return{computeERIs:(d,m,y,f)=>{const[g,v]=c(d),[w,h]=c(m);return b(a.compute_eris_wasm(g,v,w,h,y,f))},computeFockRhf:(d,m,y,f)=>{const[g,v]=c(d),[w,h]=c(m),[S,A]=c(y);return b(a.compute_fock_rhf(g,v,w,h,S,A,f))}}}let _,C;async function wt(){if(_!==void 0)return _;try{_=await Se("/GANSU-Lite/wasm/wasm_eri_bg.wasm")}catch{_=null}return _}async function vt(){if(C!==void 0)return C;try{C=await Se("/GANSU-Lite/wasm/wasm_eri_simd_bg.wasm")}catch{C=null}return C}async function $e(){return await Ce("/GANSU-Lite/"),Re()?Ie()==="wasm-simd"?{backend:"wasm",label:"WASM-SIMD"}:{backend:"wasm",label:"WASM"}:{backend:"js",label:"JS"}}function de(e,s,n,t){let o=e,r=s;if(o<r){const c=o;o=r,r=c}let a=n,i=t;if(a<i){const c=a;a=i,i=c}let u=o*(o+1)/2+r,l=a*(a+1)/2+i;if(u<l){const c=u;u=l,l=c}return u*(u+1)/2+l}function Mt(e,s,n,t){const o=new Float64Array(t*t);for(let r=0;r<t;r++)for(let a=r;a<t;a++){let i=0;for(let l=0;l<t;l++)for(let c=0;c<t;c++){const b=s[l*t+c];Math.abs(b)<1e-10||(i+=b*(e[de(r,a,l,c)]-.5*e[de(r,l,a,c)]))}const u=n[r*t+a]+i;o[r*t+a]=u,o[a*t+r]=u}return o}async function xt(e,s,n,t){const o=e.numBasis,r=new je(s,o),a=U(e,"RHF");a.computeNuclearRepulsionEnergy(),a.computeCoreHamiltonianMatrix(),a.computeTransformMatrix(),a.eri=r,a.guessInitialFockMatrix();const i=a.coreHamiltonianMatrix.data;let u=0,l=0,c=0;const b=performance.now();for(let d=0;d<200;d++){a.computeCoefficientMatrix(),a.computeDensityMatrix();const m=a.densityMatrix.data,y=await t(s,m,i,o);if(a.fockMatrix=new Ne(o,o,y),c=a.computeEnergy(),d>0&&Math.abs(c-u)<1e-8){l=d+1;break}a.updateFockMatrix(),u=c}return{eriMs:n,scfMs:performance.now()-b,iters:l,energy:c}}const X=[{molecule:"H2O",xyzFile:"H2O.xyz",basis:"6-31g(d,p)"},{molecule:"CO2",xyzFile:"CO2.xyz",basis:"6-31g(d,p)"},{molecule:"Benzene",xyzFile:"Benzene.xyz",basis:"6-31g(d,p)"}],ue=3;async function St(e){const s=await I(e.basis),n=`/GANSU-Lite/xyz/${e.xyzFile}`,t=await fetch(n);if(!t.ok)throw new Error(`Cannot load ${n}`);const o=await t.text(),r=L(o);return{mol:new G(r,s,0),xyz:o}}async function $t(e,s,n,t){const{mol:o}=await St(e),r=o.numBasis,a=o.primitiveShells,i=o.cgtoNormalizationFactors,u=(await Pe(async()=>{const{packShells:y}=await import("./ri-DT3jfd0v.js").then(f=>f.W);return{packShells:y}},[])).packShells(a),l=new Float64Array(i);async function c(y,f,g){const v=[],w=[];let h=0,S;try{for(let $=0;$<ue;$++){t(`  ${y} run ${$+1}/${ue}…`);const Be=performance.now(),Ee=await f(),Fe=performance.now()-Be,D=await xt(o,Ee,Fe,g);v.push(D.eriMs),w.push(D.scfMs),$===0&&(h=D.iters),await new Promise(_e=>setTimeout(_e,30))}}catch($){S=$ instanceof Error?$.message:String($)}const A=v.length?K(v):NaN,oe=w.length?K(w):NaN;return{eriMedianMs:A,scfMedianMs:oe,totalMedianMs:A+oe,iters:h,eriRuns:v,scfRuns:w,...S?{error:S}:{}}}const b=await c("JS",()=>Ae(a,i,r,1e-10,void 0,"js"),(y,f,g,v)=>Mt(y,f,g,v)),d=s?await c("WASM baseline",()=>s.computeERIs(u,l,r,1e-10),(y,f,g,v)=>s.computeFockRhf(y,f,g,v)):{eriMedianMs:NaN,scfMedianMs:NaN,totalMedianMs:NaN,iters:0,eriRuns:[],scfRuns:[],error:"WASM baseline not available"},m=n?await c("WASM+SIMD",()=>n.computeERIs(u,l,r,1e-10),(y,f,g,v)=>n.computeFockRhf(y,f,g,v)):{eriMedianMs:NaN,scfMedianMs:NaN,totalMedianMs:NaN,iters:0,eriRuns:[],scfRuns:[],error:"SIMD not available"};return{molecule:e.molecule,basis:e.basis,nbasis:r,js:b,wasmBase:d,wasmSimd:m}}let H=null;async function kt(){const e=document.getElementById("cb-run-compare"),s=document.getElementById("cb-stop-compare");e&&(e.disabled=!0),s&&(s.disabled=!1);const n=a=>{const i=document.getElementById("cb-compare-status");i&&(i.textContent=a)};n("Loading WASM binaries…");const t=await wt(),o=await vt(),r=[];for(let a=0;a<X.length;a++){const i=X[a];n(`[${a+1}/${X.length}] ${i.molecule}/${i.basis}…`);try{r.push(await $t(i,t,o,n)),Bt(r)}catch(u){const l=u instanceof Error?u.message:String(u);n(`  ERROR: ${l}`)}}H=r,n("Done."),e&&(e.disabled=!1),s&&(s.disabled=!0)}function Bt(e){const s=document.getElementById("cb-compare-section");if(!s)return;s.style.display="";const n=document.getElementById("cb-compare-body");if(!n)return;n.innerHTML="";const t=r=>isFinite(r)?r.toFixed(0):"—",o=(r,a)=>isFinite(r)&&isFinite(a)&&a>0?`${(r/a).toFixed(1)}×`:"—";for(const r of e){const a=document.createElement("tr");a.innerHTML=`
      <td>${r.molecule}</td><td>${r.basis}</td><td class="num">${r.nbasis}</td><td class="num">${r.js.iters||"—"}</td>
      <td class="num">${t(r.js.eriMedianMs)}</td>
      <td class="num">${t(r.js.scfMedianMs)}</td>
      <td class="num"><b>${t(r.js.totalMedianMs)}</b></td>
      <td class="num">${t(r.wasmBase.eriMedianMs)}</td>
      <td class="num">${t(r.wasmBase.scfMedianMs)}</td>
      <td class="num"><b>${t(r.wasmBase.totalMedianMs)}</b></td>
      <td class="num">${o(r.js.totalMedianMs,r.wasmBase.totalMedianMs)}</td>
      <td class="num">${t(r.wasmSimd.eriMedianMs)}</td>
      <td class="num">${t(r.wasmSimd.scfMedianMs)}</td>
      <td class="num"><b>${t(r.wasmSimd.totalMedianMs)}</b></td>
      <td class="num">${o(r.js.totalMedianMs,r.wasmSimd.totalMedianMs)}</td>
    `,n.appendChild(a)}}const Et=document.getElementById("app");let E=null,me="detecting…";function Ft(){const e=Q(),s=pe();Et.innerHTML=`
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
          <tr><th>Initial load</th><td>${s!=null?`${s} ms`:"unavailable"}</td></tr>
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

      ${ot()}

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
  `,At(),document.getElementById("cb-run").addEventListener("click",_t),document.getElementById("cb-stop").addEventListener("click",()=>{z=!0}),document.getElementById("cb-run-compare")?.addEventListener("click",kt),document.getElementById("cb-copy-compare-json")?.addEventListener("click",()=>ke(JSON.stringify(H,null,2))),document.getElementById("cb-download-compare-json")?.addEventListener("click",()=>T("backend-comparison.json",JSON.stringify(H,null,2),"application/json")),document.getElementById("cb-download-compare-csv")?.addEventListener("click",()=>T("backend-comparison.csv",Lt(H),"text/csv")),document.getElementById("cb-copy-json")?.addEventListener("click",Rt),document.getElementById("cb-download-json")?.addEventListener("click",()=>T("benchmark.json",JSON.stringify(E,null,2),"application/json")),document.getElementById("cb-download-csv")?.addEventListener("click",()=>T("benchmark.csv",It(E),"text/csv")),at(()=>me),$e().then(({label:n})=>{me=n;const t=document.getElementById("cb-backend");t&&(t.textContent=n)})}let z=!1;function R(e,s){const n=document.getElementById(`cb-row-${e}`);if(n){const t=n.querySelector(".cb-wl-status");t&&(t.textContent=s)}}async function _t(){z=!1;const e=document.getElementById("cb-run"),s=document.getElementById("cb-stop");e.disabled=!0,s.disabled=!1;const{backend:n,label:t}=await $e(),o=Q(),r=pe(),a=[];for(const i of xe){if(z){R(i.id,"stopped");continue}const u=[];let l;try{for(let d=0;d<J&&!z;d++){R(i.id,`run ${d+1}/${J}…`);const m=performance.now();await i.run(n);const y=performance.now()-m;u.push(y),R(i.id,`run ${d+1}/${J}: ${y.toFixed(0)} ms`),await new Promise(f=>setTimeout(f,50))}}catch(d){l=d instanceof Error?d.message:String(d),R(i.id,`ERROR: ${l}`)}const c=u.slice(mt),b=c.length>0?K(c):NaN;a.push({id:i.id,label:i.label,warmupMs:u[0]??NaN,runs:c,medianMs:b,...l?{error:l}:{}}),!l&&c.length>0&&R(i.id,`median ${b.toFixed(0)} ms`)}E={device:o,backend:t,initialLoadMs:r,workloads:a,toolUrl:location.href},Ct(E),e.disabled=!1,s.disabled=!0}function Ct(e){const s=document.getElementById("cb-results");s.style.display="";const n=document.getElementById("cb-results-body");n.innerHTML="";for(const t of e.workloads){const o=document.createElement("tr");o.innerHTML=`
      <td>${t.label}${t.error?` <span class="cb-err">(${t.error})</span>`:""}</td>
      <td class="num">${isFinite(t.warmupMs)?t.warmupMs.toFixed(0):"—"}</td>
      <td class="num">${t.runs[0]!=null?t.runs[0].toFixed(0):"—"}</td>
      <td class="num">${t.runs[1]!=null?t.runs[1].toFixed(0):"—"}</td>
      <td class="num">${t.runs[2]!=null?t.runs[2].toFixed(0):"—"}</td>
      <td class="num"><b>${isFinite(t.medianMs)?t.medianMs.toFixed(0):"—"}</b></td>
    `,n.appendChild(o)}document.getElementById("cb-json-pre").textContent=JSON.stringify(e,null,2)}function Rt(){E&&ke(JSON.stringify(E,null,2))}function ke(e){navigator.clipboard.writeText(e).then(()=>alert("Copied to clipboard")).catch(()=>alert("Copy failed — open browser console"))}function Lt(e){if(!e||e.length===0)return"";const s=["molecule","basis","nbasis","iters","js_eri_ms","js_scf_ms","js_total_ms","wasm_eri_ms","wasm_scf_ms","wasm_total_ms","speedup_wasm_vs_js","simd_eri_ms","simd_scf_ms","simd_total_ms","speedup_simd_vs_js"],n=r=>isFinite(r)?r.toFixed(2):"",t=(r,a)=>isFinite(r)&&isFinite(a)&&a>0?(r/a).toFixed(2):"",o=e.map(r=>[r.molecule,r.basis,String(r.nbasis),String(r.js.iters),n(r.js.eriMedianMs),n(r.js.scfMedianMs),n(r.js.totalMedianMs),n(r.wasmBase.eriMedianMs),n(r.wasmBase.scfMedianMs),n(r.wasmBase.totalMedianMs),t(r.js.totalMedianMs,r.wasmBase.totalMedianMs),n(r.wasmSimd.eriMedianMs),n(r.wasmSimd.scfMedianMs),n(r.wasmSimd.totalMedianMs),t(r.js.totalMedianMs,r.wasmSimd.totalMedianMs)].map(a=>`"${a}"`).join(","));return s.join(",")+`
`+o.join(`
`)+`
`}function T(e,s,n){const t=new Blob([s],{type:n}),o=URL.createObjectURL(t),r=document.createElement("a");r.href=o,r.download=e,r.click(),setTimeout(()=>URL.revokeObjectURL(o),1e3)}function It(e){if(!e)return"";const s=["workload_id","workload_label","backend","warmup_ms","run1_ms","run2_ms","run3_ms","median_ms","error"],n=e.workloads.map(o=>[o.id,o.label,e.backend,isFinite(o.warmupMs)?o.warmupMs.toFixed(2):"",o.runs[0]?.toFixed(2)??"",o.runs[1]?.toFixed(2)??"",o.runs[2]?.toFixed(2)??"",isFinite(o.medianMs)?o.medianMs.toFixed(2):"",o.error??""].map(r=>`"${String(r).replace(/"/g,'""')}"`).join(","));return[`# device.ua,"${e.device.ua.replace(/"/g,'""')}"`,`# device.platform,"${e.device.platform}"`,`# device.deviceMemoryGB,"${e.device.deviceMemoryGB??""}"`,`# device.hardwareConcurrency,"${e.device.hardwareConcurrency}"`,`# device.screen,"${e.device.screen}"`,`# backend,"${e.backend}"`,`# initialLoadMs,"${e.initialLoadMs??""}"`,`# timestamp,"${e.device.timestamp}"`].join(`
`)+`
`+s.join(",")+`
`+n.join(`
`)+`
`}function At(){if(document.getElementById("cb-styles"))return;const e=document.createElement("style");e.id="cb-styles",e.textContent=`
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
  `,document.head.appendChild(e)}Ft();
