import"./styles-B8WF9G-F.js";import{p as L}from"./parseXYZ-ByfvMHmk.js";import{g as G,h as U,z as Ee,o as Fe,f as _e,n as Ce,C as Re,F as Le,M as Ie}from"./ri-DT3jfd0v.js";import{c as Ae}from"./gradient-Bmg7fDyY.js";import{c as je}from"./hessian-DpHI2yIv.js";import"./properties-Odo1T7fO.js";const Oe="modulepreload",Ne=function(e){return"/GANSU-Lite/"+e},oe={},He=function(t,n,s){let o=Promise.resolve();if(n&&n.length>0){let u=function(d){return Promise.all(d.map(c=>Promise.resolve(c).then(b=>({status:"fulfilled",value:b}),b=>({status:"rejected",reason:b}))))};document.getElementsByTagName("link");const r=document.querySelector("meta[property=csp-nonce]"),i=r?.nonce||r?.getAttribute("nonce");o=u(n.map(d=>{if(d=Ne(d),d in oe)return;oe[d]=!0;const c=d.endsWith(".css"),b=c?'[rel="stylesheet"]':"";if(document.querySelector(`link[href="${d}"]${b}`))return;const l=document.createElement("link");if(l.rel=c?"stylesheet":Oe,c||(l.as="script"),l.crossOrigin="",l.href=d,i&&l.setAttribute("nonce",i),document.head.appendChild(l),c)return new Promise((m,y)=>{l.addEventListener("load",m),l.addEventListener("error",()=>y(new Error(`Unable to preload CSS for ${d}`)))})}))}function a(r){const i=new Event("vite:preloadError",{cancelable:!0});if(i.payload=r,window.dispatchEvent(i),!i.defaultPrevented)throw r}return o.then(r=>{for(const i of r||[])i.status==="rejected"&&a(i.reason);return t().catch(a)})};function Q(){const e=navigator;return{ua:e.userAgent,platform:e.platform??"unknown",deviceMemoryGB:e.deviceMemory??null,hardwareConcurrency:e.hardwareConcurrency??0,screen:`${screen.width}×${screen.height}`,language:e.language||"unknown",timestamp:new Date().toISOString()}}function me(){try{const e=performance.getEntriesByType("navigation")[0];return e?Math.round(e.loadEventEnd-e.fetchStart):null}catch{return null}}const j=.9572,O=104.52*Math.PI/180,N=3,Te=[["O",0,0,0],["H",j*Math.sin(O/2),0,j*Math.cos(O/2)],["H",-j*Math.sin(O/2),0,j*Math.cos(O/2)]];function ze(e){const t=Math.ceil(Math.cbrt(e))+1,n=[];for(let s=-t;s<=t;s++)for(let o=-t;o<=t;o++)for(let a=-t;a<=t;a++)n.push([s,o,a]);return n.sort((s,o)=>{const a=s[0]*s[0]+s[1]*s[1]+s[2]*s[2],r=o[0]*o[0]+o[1]*o[1]+o[2]*o[2];return a!==r?a-r:s[2]!==o[2]?s[2]-o[2]:s[1]!==o[1]?s[1]-o[1]:s[0]-o[0]}),n.slice(0,e)}function We(e){const t=Math.PI*(3-Math.sqrt(5));return[e*t,e*1.2345678,e*.7654321]}function Pe(e,[t,n,s]){const o=(r,i)=>[r[0]*Math.cos(i)-r[1]*Math.sin(i),r[0]*Math.sin(i)+r[1]*Math.cos(i),r[2]];return o(((r,i)=>[r[0]*Math.cos(i)+r[2]*Math.sin(i),r[1],-r[0]*Math.sin(i)+r[2]*Math.cos(i)])(o([e[0],e[1],e[2]],s),n),t)}function Ge(e){const t=ze(e),n=[];for(let s=0;s<e;s++){const[o,a,r]=t[s],i=We(s),u=o*N,d=a*N,c=r*N;for(const[b,l,m,y]of Te){const[f,g,v]=Pe([l,m,y],i);n.push(`${b.padEnd(2)} ${(u+f).toFixed(6).padStart(12)} ${(d+g).toFixed(6).padStart(12)} ${(c+v).toFixed(6).padStart(12)}`)}}return`${e*3}
(H2O)${e} — cubic lattice a=${N} A, golden-angle orientations
${n.join(`
`)}`}const ae=1.526,k=1.09,Ue=109.47*Math.PI/180;function De(e){const t=Ue/2,n=ae*Math.sin(t),s=ae*Math.cos(t),o=[];for(let c=0;c<e;c++)o.push([c*n,0,c%2===0?0:s]);const a=[];for(const c of o)a.push(F("C",c));for(let c=0;c<e;c++){const[b,l,m]=o[c],y=c%2===0?1:-1,f=m+y*k*Math.cos(t),g=k*Math.sin(t);a.push(F("H",[b,l+g,f])),a.push(F("H",[b,l-g,f]))}const r=o[0],i=o[e-1],u=[-Math.sin(t),0,e>1&&o[1][2]>r[2]?-Math.cos(t):Math.cos(t)],d=[Math.sin(t),0,e>1&&o[e-2][2]>i[2]?-Math.cos(t):Math.cos(t)];return a.push(F("H",[r[0]+k*u[0],r[1],r[2]+k*u[2]])),a.push(F("H",[i[0]+k*d[0],i[1],i[2]+k*d[2]])),`${e+2*e+2}
C${e}H${2*e+2} all-anti
${a.join(`
`)}`}function F(e,[t,n,s]){return`${e.padEnd(2)} ${t.toFixed(6).padStart(12)} ${n.toFixed(6).padStart(12)} ${s.toFixed(6).padStart(12)}`}const B=[{id:"water",label:"Water cluster (H₂O)ₙ",name:e=>`(H2O)${e}`,natoms:e=>3*e,xyz:Ge,ladder:[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32]},{id:"alkane",label:"n-Alkane CₙH₂ₙ₊₂",name:e=>`C${e}H${2*e+2}`,natoms:e=>3*e+2,xyz:De,ladder:[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20]}],ee="gansu-stress-v1",te="gansu-stress-probe-v1",pe=180,V="6-31g(d,p)";function Je(){try{const e=localStorage.getItem(ee);return e?JSON.parse(e):null}catch{return null}}function Z(e){try{localStorage.setItem(ee,JSON.stringify(e))}catch{}}function qe(){try{localStorage.removeItem(ee)}catch{}}function W(e){return e==null?"—":e<1024**2?`${(e/1024).toFixed(0)} KB`:e<1024**3?`${(e/1024**2).toFixed(0)} MB`:`${(e/1024**3).toFixed(2)} GB`}function be(e){return e==null?"—":(e/1e3).toFixed(2)}const se={ok:"ok",timeout:"timeout (exceeded budget)",error:"error / OOM","worker-died":"worker killed (OOM)","tab-crash":"tab crash"};function re(e,t,n,s,o,a){return new Promise(r=>{const i=new Worker(new URL("/GANSU-Lite/assets/stressWorker-C3Q-jLBX.js",import.meta.url),{type:"module"});a(i);let u=!1;const d=performance.now();let c=d,b="no progress reported",l=null;const m=w=>{u||(u=!0,clearTimeout(g),i.terminate(),a(null),r(w))},y=w=>{const h=w*(w+1)/2;return h*(h+1)/2*8},f=()=>({nbasis:l,totalMs:null,scfMs:null,iterations:null,converged:null,energy:null,eriBytes:l==null?null:y(l),heapUsedBytes:null,heapLimitBytes:null}),g=setTimeout(()=>{const w=(performance.now()-c)/1e3;m({...f(),status:"timeout",totalMs:performance.now()-d,detail:`exceeded ${(n/1e3).toFixed(0)} s budget; silent for ${w.toFixed(0)} s after "${b}"`+(w>30?" — the worker had gone quiet, so it may have been killed rather than merely slow":" — the worker was still reporting progress, so this is genuine slowness, not a crash")})},n);i.onmessage=w=>{const h=w.data;if(c=performance.now(),h.type==="progress"){const S=/(\d+) basis functions/.exec(h.message);S&&(l=Number(S[1])),b=h.message,o(h.message,h.elapsedMs);return}if(h.type==="done"){m({status:"ok",nbasis:h.nbasis,totalMs:h.totalMs,scfMs:h.scfMs,iterations:h.iterations,converged:h.converged,energy:h.energy,eriBytes:h.eriBytes,heapUsedBytes:h.heapUsedBytes,heapLimitBytes:h.heapLimitBytes,...h.converged?{}:{detail:"SCF did not converge"}});return}h.type==="error"&&m({...f(),status:"error",nbasis:h.nbasis??l,totalMs:h.elapsedMs,detail:`${h.name}: ${h.message} (phase: ${h.phase})`})},i.onerror=w=>{m({...f(),status:"worker-died",totalMs:performance.now()-d,detail:w.message||"worker terminated by the browser"})};const v={type:"stress-run",xyzText:e,basisGBS:t,charge:0,eriBackend:"auto",baseUrl:"/GANSU-Lite/",...s?{dftConfig:{functional:"B3LYP",gridLevel:"medium"}}:{}};i.postMessage(v)})}function Ye(){const e=new Set;for(const t of B[0].ladder)e.add(t*25);return[...e].filter(t=>t<=300).sort((t,n)=>t-n)}function Xe(e){return new Promise(t=>{const n=new Worker(new URL("/GANSU-Lite/assets/stressWorker-C3Q-jLBX.js",import.meta.url),{type:"module"}),s=[],o=()=>{n.terminate(),t(s)};n.onmessage=r=>{const i=r.data;if(i.type==="probe-result"){const u={nbasis:i.nbasis,bytes:i.bytes,ok:i.ok,ms:i.ms,verdict:i.verdict,...i.detail?{detail:i.detail}:{}};s.push(u),e(u);return}i.type==="probe-done"&&o()},n.onerror=()=>o();const a={type:"memory-probe",nbasisLadder:Ye(),baseUrl:"/GANSU-Lite/"};n.postMessage(a)})}function fe(e,t){const n=t.verdict==="ok",s={ok:"allocated in WASM memory",failed:"ALLOCATION FAILED",unavailable:"WASM unavailable"}[t.verdict],o=document.createElement("tr");o.innerHTML=`<td class="num">${t.nbasis||"—"}</td><td class="num">${t.bytes?W(t.bytes):"—"}</td><td class="num">${t.ms.toFixed(0)}</td><td${n?"":' class="cb-err"'}>${s}${t.detail?`<br><span class="cb-err">${he(t.detail)}</span>`:""}</td>`,e.appendChild(o)}let x=[];function Ve(e){try{localStorage.setItem(te,JSON.stringify(e))}catch{}}function Ze(){try{const e=localStorage.getItem(te);return e?JSON.parse(e):[]}catch{return[]}}async function Ke(){const e=document.getElementById("cb-stress-probe"),t=document.getElementById("cb-probe-body"),n=document.getElementById("cb-probe-section");e&&(e.disabled=!0),n&&(n.style.display=""),t&&(t.innerHTML=""),x=[],M("probing memory ceiling…"),await Xe(r=>{x.push(r),Ve(x),t&&fe(t,r)});const s=x.filter(r=>r.verdict==="ok"),o=x.find(r=>r.verdict!=="ok"),a=document.getElementById("cb-probe-summary");if(a){const r=s.length?s[s.length-1]:null,i=x.find(u=>u.verdict==="unavailable");if(i)a.textContent=`No WASM on this device, so there is no wasm32 ceiling to measure — it runs the JS backend, whose limit is device RAM and shows up in the ladder above. ${i.detail??""}`;else if(!r)a.textContent="Could not allocate even the smallest probe array — something is wrong.";else{const u=o?{failed:`the allocator refused ${o.nbasis} (${W(o.bytes)})`,unavailable:"",ok:""}[o.verdict]:"the probe ladder topped out at 300 basis functions without hitting a ceiling";a.textContent=`WASM memory ceiling: ${o?"":"≥ "}${r.nbasis} basis functions — a ${W(r.bytes)} ERI array allocates inside the module's linear memory; ${u}.`}}M("memory probe done"),e&&(e.disabled=!1)}function Qe(e,t){if(e.length<2||t.nbasis==null)return null;const n=e[e.length-2],s=e[e.length-1];if(n.nbasis==null||s.nbasis==null||n.totalMs==null||s.totalMs==null||n.nbasis>=s.nbasis||s.nbasis>=t.nbasis||n.totalMs<=0)return null;const o=Math.log(s.totalMs/n.totalMs)/Math.log(s.nbasis/n.nbasis);return!Number.isFinite(o)||o<=0?null:{seconds:s.totalMs*Math.pow(t.nbasis/s.nbasis,o)/1e3,exponent:o}}function et(e,t){const n=e.filter(r=>r.status==="ok"),s=n.length?n[n.length-1]:null,o=e.find(r=>r.status!=="ok")??null;let a=null;if(o&&o.status==="timeout"){const r=Qe(n,o);r&&(a={estimatedSeconds:r.seconds,budgetSeconds:t,exponent:r.exponent,fragile:r.seconds<1.25*t})}return{borderline:a,largestCompleted:s&&s.nbasis!=null&&s.totalMs!=null?{name:s.name,natoms:s.natoms,nbasis:s.nbasis,seconds:+(s.totalMs/1e3).toFixed(1)}:null,firstFailure:o?{name:o.name,natoms:o.natoms,nbasis:o.nbasis,mode:se[o.status]}:null}}let p={series:B[0],points:[],running:!1,stopRequested:!1,activeWorker:null,output:null};function tt(){return`
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
          <input id="cb-stress-budget" type="number" min="10" max="3600" step="10" value="${pe}">
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
        <pre id="cb-stress-row" class="cb-stress-row"></pre>
        <button id="cb-stress-copy-row">Copy row</button>
        <button id="cb-stress-copy-json">Copy JSON</button>
        <button id="cb-stress-download-json">Download JSON</button>
        <button id="cb-stress-download-csv">Download CSV</button>
      </div>
    </section>
  `}function st(e){const t=a=>document.getElementById(a),n=t("cb-stress-series"),s=t("cb-stress-ladder");n?.addEventListener("change",()=>{const a=B.find(r=>r.id===n.value);a&&s&&(p.series=a,s.value=a.ladder.join(", "))}),t("cb-stress-run")?.addEventListener("click",()=>{rt(e)}),t("cb-stress-stop")?.addEventListener("click",()=>{p.stopRequested=!0,p.activeWorker?.terminate(),M("stopping…")}),t("cb-stress-reset")?.addEventListener("click",()=>{qe();try{localStorage.removeItem(te)}catch{}x=[],p.points=[],p.output=null,P();const a=t("cb-stress-recovered");a&&(a.style.display="none");const r=t("cb-stress-summary");r&&(r.style.display="none"),M("saved progress cleared")}),t("cb-stress-probe")?.addEventListener("click",()=>{Ke()}),t("cb-stress-copy-row")?.addEventListener("click",()=>ie(ye())),t("cb-stress-copy-json")?.addEventListener("click",()=>ie(JSON.stringify(p.output,null,2))),t("cb-stress-download-json")?.addEventListener("click",()=>ce("stress-test.json",JSON.stringify(p.output,null,2),"application/json")),t("cb-stress-download-csv")?.addEventListener("click",()=>ce("stress-test.csv",lt(p.output),"text/csv")),nt();const o=Ze();if(o.length){x=o;const a=document.getElementById("cb-probe-section"),r=document.getElementById("cb-probe-body");if(a&&r){a.style.display="";for(const i of o)fe(r,i)}}}function nt(){const e=Je();if(!e)return;p.points=e.points??[];const t=B.find(n=>n.id===e.seriesId);if(t&&(p.series=t),e.inFlight){const n=e.inFlight;p.points.some(o=>o.n===n.n)||p.points.push({n:n.n,name:n.name,natoms:n.natoms,nbasis:null,status:"tab-crash",totalMs:null,scfMs:null,iterations:null,converged:null,energy:null,eriBytes:null,heapUsedBytes:null,heapLimitBytes:null,detail:"tab died while this point was running"});const s=document.getElementById("cb-stress-recovered");s&&(s.style.display="",s.textContent=`Recovered from a previous session: the tab died while running ${n.name}. Results up to that point are restored below and ${n.name} is recorded as a tab crash.`),Z({...e,points:p.points,inFlight:null})}p.points.length&&P()}function M(e){const t=document.getElementById("cb-stress-status");t&&(t.textContent=e)}async function ot(){const e=await fetch(`/GANSU-Lite/basis/${V}.gbs`);if(!e.ok)throw new Error(`Cannot load basis ${V}: HTTP ${e.status}`);return e.text()}function at(e,t){const n=e.split(/[,\s]+/).map(s=>parseInt(s,10)).filter(s=>Number.isFinite(s)&&s>0);return n.length?n:t}async function rt(e){if(p.running)return;const t=document.getElementById("cb-stress-run"),n=document.getElementById("cb-stress-stop"),s=Math.max(10,Number(document.getElementById("cb-stress-budget")?.value)||pe),o=document.getElementById("cb-stress-theory")?.value==="b3lyp",a=at(document.getElementById("cb-stress-ladder")?.value??"",p.series.ladder);p.running=!0,p.stopRequested=!1,p.points=[],p.output=null,t&&(t.disabled=!0),n&&(n.disabled=!1),P();let r;try{M("loading basis set…"),r=await ot()}catch(c){M(`FAILED: ${c instanceof Error?c.message:String(c)}`),p.running=!1,t&&(t.disabled=!1),n&&(n.disabled=!0);return}const i=p.series,u=o?"B3LYP":"RHF",d=s*1e3;try{M("warm-up…"),await re(i.xyz(a[0]),r,d,o,()=>{},c=>{p.activeWorker=c})}catch{}for(const c of a){if(p.stopRequested){M("stopped");break}const b=i.name(c),l=i.natoms(c);Z({seriesId:i.id,method:u,budgetSeconds:s,points:p.points,inFlight:{n:c,name:b,natoms:l,startedAt:Date.now()}}),M(`${b} — starting…`);const m=await re(i.xyz(c),r,d,o,(f,g)=>M(`${b} — ${f} (${(g/1e3).toFixed(1)} s)`),f=>{p.activeWorker=f}),y={n:c,name:b,natoms:l,...m};if(p.points.push(y),P(),Z({seriesId:i.id,method:u,budgetSeconds:s,points:p.points,inFlight:null}),m.status!=="ok"){M(`${b} failed: ${se[m.status]} — stopping ladder`);break}M(`${b} ok — ${be(m.totalMs)} s`),await new Promise(f=>setTimeout(f,200))}p.output={device:Q(),backend:e(),seriesId:i.id,seriesLabel:i.label,method:u,basis:V,budgetSeconds:s,points:p.points,summary:et(p.points,s),memoryProbe:x,toolUrl:location.href},ct(),p.running=!1,p.activeWorker=null,t&&(t.disabled=!1),n&&(n.disabled=!0),p.stopRequested||M("done")}function P(){const e=document.getElementById("cb-stress-body");if(e){e.innerHTML="";for(const t of p.points){const n=document.createElement("tr"),s=t.status!=="ok";n.innerHTML=`
      <td>${t.name}</td>
      <td class="num">${t.natoms}</td>
      <td class="num">${t.nbasis??"—"}</td>
      <td class="num">${W(t.eriBytes)}</td>
      <td class="num">${t.iterations??"—"}</td>
      <td class="num">${be(t.totalMs)}</td>
      <td${s?' class="cb-err"':""}>${se[t.status]}${t.detail?`<br><span class="cb-err">${he(t.detail)}</span>`:""}</td>
    `,e.appendChild(n)}}}function he(e){return e.replace(/[&<>"]/g,t=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"})[t])}function ye(){const e=p.output;if(!e)return"";const t=e.summary.largestCompleted,n=e.summary.firstFailure;return`| ${[it(e.device),t?t.name:"—",t?String(t.natoms):"—",t?String(t.nbasis):"—",t?t.seconds.toFixed(1):"—",n?n.name:"(no failure in ladder)",n?n.mode:"—"].join(" | ")} |`}function it(e){const t=e.ua;return/iPhone/.test(t)?"iPhone":/iPad/.test(t)?"iPad":/Mac OS X|Macintosh/.test(t)?"Apple-silicon laptop":/Android/.test(t)?"Android":/Windows/.test(t)?"Windows":e.platform}function ct(){const e=document.getElementById("cb-stress-summary"),t=document.getElementById("cb-stress-row");if(!e||!t||!p.output)return;e.style.display="";const n=p.output.summary,s="| Device | Largest completed | Atoms | Basis fns | Time (s) | First failing size | Failure mode |",o="|---|---|---|---|---|---|---|",a=[];if(n.firstFailure||a.push("# The ladder finished without a failure — extend it to find the limit."),n.borderline){const i=n.borderline;a.push(`# ${n.firstFailure?.name} is extrapolated to need ~${i.estimatedSeconds.toFixed(0)} s (fitted N^${i.exponent.toFixed(1)}) against a ${i.budgetSeconds} s budget.`),i.fragile&&a.push("# FRAGILE: that is within 25% of the budget, so run-to-run noise alone could flip it.",`# Treat this row as "limit at a ${i.budgetSeconds} s budget", not as the device's ceiling.`,`# For the device's own ceiling, re-run with the budget at ${Math.ceil(i.estimatedSeconds*1.5/60)*60} s or more.`)}const r=p.points.filter(i=>i.status==="ok"&&i.converged===!1);r.length&&a.push(`# SCF did not converge for: ${r.map(i=>i.name).join(", ")}`),t.textContent=[s,o,ye(),...a].join(`
`)}function ie(e){navigator.clipboard.writeText(e).then(()=>alert("Copied to clipboard")).catch(()=>alert("Copy failed — open browser console"))}function ce(e,t,n){const s=new Blob([t],{type:n}),o=URL.createObjectURL(s),a=document.createElement("a");a.href=o,a.download=e,a.click(),setTimeout(()=>URL.revokeObjectURL(o),1e3)}function lt(e){if(!e)return"";const t=["n","system","atoms","nbasis","eri_bytes","iterations","converged","setup_plus_scf_ms","scf_ms","energy_hartree","status","detail"],n=e.points.map(o=>[o.n,o.name,o.natoms,o.nbasis??"",o.eriBytes??"",o.iterations??"",o.converged==null?"":String(o.converged),o.totalMs?.toFixed(2)??"",o.scfMs?.toFixed(2)??"",o.energy?.toFixed(8)??"",o.status,o.detail??""].map(a=>`"${String(a).replace(/"/g,'""')}"`).join(","));return[`# device.ua,"${e.device.ua.replace(/"/g,'""')}"`,`# device.platform,"${e.device.platform}"`,`# device.deviceMemoryGB,"${e.device.deviceMemoryGB??""}"`,`# device.hardwareConcurrency,"${e.device.hardwareConcurrency}"`,`# backend,"${e.backend}"`,`# series,"${e.seriesLabel}"`,`# method,"${e.method}/${e.basis}"`,`# budget_seconds,"${e.budgetSeconds}"`,`# timestamp,"${e.device.timestamp}"`].join(`
`)+`
`+t.join(",")+`
`+n.join(`
`)+`
`}const J=4,dt=1,q=1.8897259886,ge=`3
H2O standard
O   0.000000   0.000000   0.117790
H   0.000000   0.755450  -0.471160
H   0.000000  -0.755450  -0.471160`,ut=`3
H2O displaced
O   0.000000   0.000000   0.000000
H   0.000000   0.850000   0.620000
H   0.000000  -0.850000   0.620000`,mt=`3
CO2 linear
C   0.000000   0.000000   0.000000
O   0.000000   0.000000   1.160000
O   0.000000   0.000000  -1.160000`,Y=new Map;async function I(e){if(Y.has(e))return Y.get(e);const t=`/GANSU-Lite/basis/${e}.gbs`,n=await fetch(t);if(!n.ok)throw new Error(`Failed to load basis ${e}: ${n.status}`);const s=_e.fromGBS(await n.text());return Y.set(e,s),s}async function pt(e){const t=await I("6-31g(d,p)"),n=L(ge),s=new G(n,t,0);await U(s,"RHF").solve({eriBackend:e})}async function bt(e){const t=await I("sto-3g"),n=L(ut);let s=new Float64Array(n.length*3);for(let r=0;r<n.length;r++)s[3*r]=n[r].coordinate.x,s[3*r+1]=n[r].coordinate.y,s[3*r+2]=n[r].coordinate.z;const o=.4,a=25;for(let r=0;r<a;r++){const i=n.map((l,m)=>({...l,coordinate:{x:s[3*m],y:s[3*m+1],z:s[3*m+2]}})),u=new G(i,t,0),d=U(u,"RHF");await d.solve({eriBackend:e});const c=Ae(u.primitiveShells,u.atoms,u.cgtoNormalizationFactors,u.numBasis,u.numAlphaSpins,d.density,d.coefficients,d.orbitalEnergies).total;let b=0;for(let l=0;l<c.length;l++)Math.abs(c[l])>b&&(b=Math.abs(c[l]));if(b<5e-4)break;for(let l=0;l<c.length;l++)s[l]-=o*c[l]}}async function ft(e){const t=await I("sto-3g"),n=L(mt),s=new Float64Array(n.length*3);for(let o=0;o<n.length;o++)s[3*o]=n[o].coordinate.x*q,s[3*o+1]=n[o].coordinate.y*q,s[3*o+2]=n[o].coordinate.z*q;await je(n.map(o=>o.atomicNumber),s,t,0,5e-4,void 0,void 0)}async function ht(e){const t=await I("6-31g(d,p)"),n=L(ge),s=new G(n,t,0);await U(s,"RHF",{functional:"B3LYP",gridLevel:"medium"}).solve({eriBackend:e})}const we=[{id:"h2o_hf_631gdp",label:"H2O HF/6-31G(d,p)",run:pt},{id:"h2o_geomopt_sto3g",label:"H2O geom opt HF/STO-3G",run:bt},{id:"co2_hess_sto3g",label:"CO2 Hessian HF/STO-3G",run:ft},{id:"h2o_b3lyp_631gdp",label:"H2O B3LYP/6-31G(d,p)",run:ht}];function K(e){const t=[...e].sort((s,o)=>s-o),n=t.length;return n===0?NaN:n%2===1?t[(n-1)/2]:(t[n/2-1]+t[n/2])/2}async function ve(e){const t=await fetch(e);if(!t.ok)throw new Error(`Cannot fetch ${e}`);const n=await t.arrayBuffer(),s=await WebAssembly.compile(n);let o;const a={"./wasm_eri_bg.js":{__wbindgen_init_externref_table:()=>{const l=o.exports.__wbindgen_externrefs,m=l.grow(4);l.set(0,void 0),l.set(m+0,void 0),l.set(m+1,null),l.set(m+2,!0),l.set(m+3,!1)}}};o=await WebAssembly.instantiate(s,a);const r=o.exports;r.__wbindgen_start&&r.__wbindgen_start();const i=r.memory,u=r.__wbindgen_malloc,d=r.__wbindgen_free,c=l=>{const m=u(l.length*8,8)>>>0;return new Float64Array(i.buffer).set(l,m/8),[m,l.length]},b=l=>{const[m,y]=l,f=new Float64Array(i.buffer).slice(m/8,m/8+y);return d(m,y*8,8),f};return{computeERIs:(l,m,y,f)=>{const[g,v]=c(l),[w,h]=c(m);return b(r.compute_eris_wasm(g,v,w,h,y,f))},computeFockRhf:(l,m,y,f)=>{const[g,v]=c(l),[w,h]=c(m),[S,A]=c(y);return b(r.compute_fock_rhf(g,v,w,h,S,A,f))}}}let _,C;async function yt(){if(_!==void 0)return _;try{_=await ve("/GANSU-Lite/wasm/wasm_eri_bg.wasm")}catch{_=null}return _}async function gt(){if(C!==void 0)return C;try{C=await ve("/GANSU-Lite/wasm/wasm_eri_simd_bg.wasm")}catch{C=null}return C}async function Me(){return await Ee("/GANSU-Lite/"),Fe()?Ce()==="wasm-simd"?{backend:"wasm",label:"WASM-SIMD"}:{backend:"wasm",label:"WASM"}:{backend:"js",label:"JS"}}function le(e,t,n,s){let o=e,a=t;if(o<a){const c=o;o=a,a=c}let r=n,i=s;if(r<i){const c=r;r=i,i=c}let u=o*(o+1)/2+a,d=r*(r+1)/2+i;if(u<d){const c=u;u=d,d=c}return u*(u+1)/2+d}function wt(e,t,n,s){const o=new Float64Array(s*s);for(let a=0;a<s;a++)for(let r=a;r<s;r++){let i=0;for(let d=0;d<s;d++)for(let c=0;c<s;c++){const b=t[d*s+c];Math.abs(b)<1e-10||(i+=b*(e[le(a,r,d,c)]-.5*e[le(a,d,r,c)]))}const u=n[a*s+r]+i;o[a*s+r]=u,o[r*s+a]=u}return o}async function vt(e,t,n,s){const o=e.numBasis,a=new Le(t,o),r=U(e,"RHF");r.computeNuclearRepulsionEnergy(),r.computeCoreHamiltonianMatrix(),r.computeTransformMatrix(),r.eri=a,r.guessInitialFockMatrix();const i=r.coreHamiltonianMatrix.data;let u=0,d=0,c=0;const b=performance.now();for(let l=0;l<200;l++){r.computeCoefficientMatrix(),r.computeDensityMatrix();const m=r.densityMatrix.data,y=await s(t,m,i,o);if(r.fockMatrix=new Ie(o,o,y),c=r.computeEnergy(),l>0&&Math.abs(c-u)<1e-8){d=l+1;break}r.updateFockMatrix(),u=c}return{eriMs:n,scfMs:performance.now()-b,iters:d,energy:c}}const X=[{molecule:"H2O",xyzFile:"H2O.xyz",basis:"6-31g(d,p)"},{molecule:"CO2",xyzFile:"CO2.xyz",basis:"6-31g(d,p)"},{molecule:"Benzene",xyzFile:"Benzene.xyz",basis:"6-31g(d,p)"}],de=3;async function Mt(e){const t=await I(e.basis),n=`/GANSU-Lite/xyz/${e.xyzFile}`,s=await fetch(n);if(!s.ok)throw new Error(`Cannot load ${n}`);const o=await s.text(),a=L(o);return{mol:new G(a,t,0),xyz:o}}async function xt(e,t,n,s){const{mol:o}=await Mt(e),a=o.numBasis,r=o.primitiveShells,i=o.cgtoNormalizationFactors,u=(await He(async()=>{const{packShells:y}=await import("./ri-DT3jfd0v.js").then(f=>f.W);return{packShells:y}},[])).packShells(r),d=new Float64Array(i);async function c(y,f,g){const v=[],w=[];let h=0,S;try{for(let $=0;$<de;$++){s(`  ${y} run ${$+1}/${de}…`);const Se=performance.now(),$e=await f(),ke=performance.now()-Se,D=await vt(o,$e,ke,g);v.push(D.eriMs),w.push(D.scfMs),$===0&&(h=D.iters),await new Promise(Be=>setTimeout(Be,30))}}catch($){S=$ instanceof Error?$.message:String($)}const A=v.length?K(v):NaN,ne=w.length?K(w):NaN;return{eriMedianMs:A,scfMedianMs:ne,totalMedianMs:A+ne,iters:h,eriRuns:v,scfRuns:w,...S?{error:S}:{}}}const b=await c("JS",()=>Re(r,i,a,1e-10,void 0,"js"),(y,f,g,v)=>wt(y,f,g,v)),l=t?await c("WASM baseline",()=>t.computeERIs(u,d,a,1e-10),(y,f,g,v)=>t.computeFockRhf(y,f,g,v)):{eriMedianMs:NaN,scfMedianMs:NaN,totalMedianMs:NaN,iters:0,eriRuns:[],scfRuns:[],error:"WASM baseline not available"},m=n?await c("WASM+SIMD",()=>n.computeERIs(u,d,a,1e-10),(y,f,g,v)=>n.computeFockRhf(y,f,g,v)):{eriMedianMs:NaN,scfMedianMs:NaN,totalMedianMs:NaN,iters:0,eriRuns:[],scfRuns:[],error:"SIMD not available"};return{molecule:e.molecule,basis:e.basis,nbasis:a,js:b,wasmBase:l,wasmSimd:m}}let T=null;async function St(){const e=document.getElementById("cb-run-compare"),t=document.getElementById("cb-stop-compare");e&&(e.disabled=!0),t&&(t.disabled=!1);const n=r=>{const i=document.getElementById("cb-compare-status");i&&(i.textContent=r)};n("Loading WASM binaries…");const s=await yt(),o=await gt(),a=[];for(let r=0;r<X.length;r++){const i=X[r];n(`[${r+1}/${X.length}] ${i.molecule}/${i.basis}…`);try{a.push(await xt(i,s,o,n)),$t(a)}catch(u){const d=u instanceof Error?u.message:String(u);n(`  ERROR: ${d}`)}}T=a,n("Done."),e&&(e.disabled=!1),t&&(t.disabled=!0)}function $t(e){const t=document.getElementById("cb-compare-section");if(!t)return;t.style.display="";const n=document.getElementById("cb-compare-body");if(!n)return;n.innerHTML="";const s=a=>isFinite(a)?a.toFixed(0):"—",o=(a,r)=>isFinite(a)&&isFinite(r)&&r>0?`${(a/r).toFixed(1)}×`:"—";for(const a of e){const r=document.createElement("tr");r.innerHTML=`
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
    `,n.appendChild(r)}}const kt=document.getElementById("app");let E=null,ue="detecting…";function Bt(){const e=Q(),t=me();kt.innerHTML=`
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
          ${we.map(n=>`<li id="cb-row-${n.id}"><span class="cb-wl-label">${n.label}</span><span class="cb-wl-status">pending</span></li>`).join("")}
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

      ${tt()}

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
  `,Lt(),document.getElementById("cb-run").addEventListener("click",Et),document.getElementById("cb-stop").addEventListener("click",()=>{z=!0}),document.getElementById("cb-run-compare")?.addEventListener("click",St),document.getElementById("cb-copy-compare-json")?.addEventListener("click",()=>xe(JSON.stringify(T,null,2))),document.getElementById("cb-download-compare-json")?.addEventListener("click",()=>H("backend-comparison.json",JSON.stringify(T,null,2),"application/json")),document.getElementById("cb-download-compare-csv")?.addEventListener("click",()=>H("backend-comparison.csv",Ct(T),"text/csv")),document.getElementById("cb-copy-json")?.addEventListener("click",_t),document.getElementById("cb-download-json")?.addEventListener("click",()=>H("benchmark.json",JSON.stringify(E,null,2),"application/json")),document.getElementById("cb-download-csv")?.addEventListener("click",()=>H("benchmark.csv",Rt(E),"text/csv")),st(()=>ue),Me().then(({label:n})=>{ue=n;const s=document.getElementById("cb-backend");s&&(s.textContent=n)})}let z=!1;function R(e,t){const n=document.getElementById(`cb-row-${e}`);if(n){const s=n.querySelector(".cb-wl-status");s&&(s.textContent=t)}}async function Et(){z=!1;const e=document.getElementById("cb-run"),t=document.getElementById("cb-stop");e.disabled=!0,t.disabled=!1;const{backend:n,label:s}=await Me(),o=Q(),a=me(),r=[];for(const i of we){if(z){R(i.id,"stopped");continue}const u=[];let d;try{for(let l=0;l<J&&!z;l++){R(i.id,`run ${l+1}/${J}…`);const m=performance.now();await i.run(n);const y=performance.now()-m;u.push(y),R(i.id,`run ${l+1}/${J}: ${y.toFixed(0)} ms`),await new Promise(f=>setTimeout(f,50))}}catch(l){d=l instanceof Error?l.message:String(l),R(i.id,`ERROR: ${d}`)}const c=u.slice(dt),b=c.length>0?K(c):NaN;r.push({id:i.id,label:i.label,warmupMs:u[0]??NaN,runs:c,medianMs:b,...d?{error:d}:{}}),!d&&c.length>0&&R(i.id,`median ${b.toFixed(0)} ms`)}E={device:o,backend:s,initialLoadMs:a,workloads:r,toolUrl:location.href},Ft(E),e.disabled=!1,t.disabled=!0}function Ft(e){const t=document.getElementById("cb-results");t.style.display="";const n=document.getElementById("cb-results-body");n.innerHTML="";for(const s of e.workloads){const o=document.createElement("tr");o.innerHTML=`
      <td>${s.label}${s.error?` <span class="cb-err">(${s.error})</span>`:""}</td>
      <td class="num">${isFinite(s.warmupMs)?s.warmupMs.toFixed(0):"—"}</td>
      <td class="num">${s.runs[0]!=null?s.runs[0].toFixed(0):"—"}</td>
      <td class="num">${s.runs[1]!=null?s.runs[1].toFixed(0):"—"}</td>
      <td class="num">${s.runs[2]!=null?s.runs[2].toFixed(0):"—"}</td>
      <td class="num"><b>${isFinite(s.medianMs)?s.medianMs.toFixed(0):"—"}</b></td>
    `,n.appendChild(o)}document.getElementById("cb-json-pre").textContent=JSON.stringify(e,null,2)}function _t(){E&&xe(JSON.stringify(E,null,2))}function xe(e){navigator.clipboard.writeText(e).then(()=>alert("Copied to clipboard")).catch(()=>alert("Copy failed — open browser console"))}function Ct(e){if(!e||e.length===0)return"";const t=["molecule","basis","nbasis","iters","js_eri_ms","js_scf_ms","js_total_ms","wasm_eri_ms","wasm_scf_ms","wasm_total_ms","speedup_wasm_vs_js","simd_eri_ms","simd_scf_ms","simd_total_ms","speedup_simd_vs_js"],n=a=>isFinite(a)?a.toFixed(2):"",s=(a,r)=>isFinite(a)&&isFinite(r)&&r>0?(a/r).toFixed(2):"",o=e.map(a=>[a.molecule,a.basis,String(a.nbasis),String(a.js.iters),n(a.js.eriMedianMs),n(a.js.scfMedianMs),n(a.js.totalMedianMs),n(a.wasmBase.eriMedianMs),n(a.wasmBase.scfMedianMs),n(a.wasmBase.totalMedianMs),s(a.js.totalMedianMs,a.wasmBase.totalMedianMs),n(a.wasmSimd.eriMedianMs),n(a.wasmSimd.scfMedianMs),n(a.wasmSimd.totalMedianMs),s(a.js.totalMedianMs,a.wasmSimd.totalMedianMs)].map(r=>`"${r}"`).join(","));return t.join(",")+`
`+o.join(`
`)+`
`}function H(e,t,n){const s=new Blob([t],{type:n}),o=URL.createObjectURL(s),a=document.createElement("a");a.href=o,a.download=e,a.click(),setTimeout(()=>URL.revokeObjectURL(o),1e3)}function Rt(e){if(!e)return"";const t=["workload_id","workload_label","backend","warmup_ms","run1_ms","run2_ms","run3_ms","median_ms","error"],n=e.workloads.map(o=>[o.id,o.label,e.backend,isFinite(o.warmupMs)?o.warmupMs.toFixed(2):"",o.runs[0]?.toFixed(2)??"",o.runs[1]?.toFixed(2)??"",o.runs[2]?.toFixed(2)??"",isFinite(o.medianMs)?o.medianMs.toFixed(2):"",o.error??""].map(a=>`"${String(a).replace(/"/g,'""')}"`).join(","));return[`# device.ua,"${e.device.ua.replace(/"/g,'""')}"`,`# device.platform,"${e.device.platform}"`,`# device.deviceMemoryGB,"${e.device.deviceMemoryGB??""}"`,`# device.hardwareConcurrency,"${e.device.hardwareConcurrency}"`,`# device.screen,"${e.device.screen}"`,`# backend,"${e.backend}"`,`# initialLoadMs,"${e.initialLoadMs??""}"`,`# timestamp,"${e.device.timestamp}"`].join(`
`)+`
`+t.join(",")+`
`+n.join(`
`)+`
`}function Lt(){if(document.getElementById("cb-styles"))return;const e=document.createElement("style");e.id="cb-styles",e.textContent=`
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
  `,document.head.appendChild(e)}Bt();
