import"./styles-B8WF9G-F.js";import{p as L}from"./parseXYZ-ByfvMHmk.js";import{g as G,h as W,z as Et,o as _t,f as Ft,n as Ct,C as Rt,F as Lt,M as It}from"./ri-DT3jfd0v.js";import{c as Ot}from"./gradient-Bmg7fDyY.js";import{c as jt}from"./hessian-DpHI2yIv.js";import"./properties-Odo1T7fO.js";const At="modulepreload",Nt=function(t){return"/GANSU-Lite/"+t},ot={},Ht=function(e,s,n){let o=Promise.resolve();if(s&&s.length>0){let u=function(d){return Promise.all(d.map(c=>Promise.resolve(c).then(f=>({status:"fulfilled",value:f}),f=>({status:"rejected",reason:f}))))};document.getElementsByTagName("link");const r=document.querySelector("meta[property=csp-nonce]"),i=r?.nonce||r?.getAttribute("nonce");o=u(s.map(d=>{if(d=Nt(d),d in ot)return;ot[d]=!0;const c=d.endsWith(".css"),f=c?'[rel="stylesheet"]':"";if(document.querySelector(`link[href="${d}"]${f}`))return;const l=document.createElement("link");if(l.rel=c?"stylesheet":At,c||(l.as="script"),l.crossOrigin="",l.href=d,i&&l.setAttribute("nonce",i),document.head.appendChild(l),c)return new Promise((m,h)=>{l.addEventListener("load",m),l.addEventListener("error",()=>h(new Error(`Unable to preload CSS for ${d}`)))})}))}function a(r){const i=new Event("vite:preloadError",{cancelable:!0});if(i.payload=r,window.dispatchEvent(i),!i.defaultPrevented)throw r}return o.then(r=>{for(const i of r||[])i.status==="rejected"&&a(i.reason);return e().catch(a)})};function K(){const t=navigator;return{ua:t.userAgent,platform:t.platform??"unknown",deviceMemoryGB:t.deviceMemory??null,hardwareConcurrency:t.hardwareConcurrency??0,screen:`${screen.width}×${screen.height}`,language:t.language||"unknown",timestamp:new Date().toISOString()}}function mt(){try{const t=performance.getEntriesByType("navigation")[0];return t?Math.round(t.loadEventEnd-t.fetchStart):null}catch{return null}}const j=.9572,A=104.52*Math.PI/180,N=3,Tt=[["O",0,0,0],["H",j*Math.sin(A/2),0,j*Math.cos(A/2)],["H",-j*Math.sin(A/2),0,j*Math.cos(A/2)]];function zt(t){const e=Math.ceil(Math.cbrt(t))+1,s=[];for(let n=-e;n<=e;n++)for(let o=-e;o<=e;o++)for(let a=-e;a<=e;a++)s.push([n,o,a]);return s.sort((n,o)=>{const a=n[0]*n[0]+n[1]*n[1]+n[2]*n[2],r=o[0]*o[0]+o[1]*o[1]+o[2]*o[2];return a!==r?a-r:n[2]!==o[2]?n[2]-o[2]:n[1]!==o[1]?n[1]-o[1]:n[0]-o[0]}),s.slice(0,t)}function Pt(t){const e=Math.PI*(3-Math.sqrt(5));return[t*e,t*1.2345678,t*.7654321]}function Gt(t,[e,s,n]){const o=(r,i)=>[r[0]*Math.cos(i)-r[1]*Math.sin(i),r[0]*Math.sin(i)+r[1]*Math.cos(i),r[2]];return o(((r,i)=>[r[0]*Math.cos(i)+r[2]*Math.sin(i),r[1],-r[0]*Math.sin(i)+r[2]*Math.cos(i)])(o([t[0],t[1],t[2]],n),s),e)}function Wt(t){const e=zt(t),s=[];for(let n=0;n<t;n++){const[o,a,r]=e[n],i=Pt(n),u=o*N,d=a*N,c=r*N;for(const[f,l,m,h]of Tt){const[y,g,b]=Gt([l,m,h],i);s.push(`${f.padEnd(2)} ${(u+y).toFixed(6).padStart(12)} ${(d+g).toFixed(6).padStart(12)} ${(c+b).toFixed(6).padStart(12)}`)}}return`${t*3}
(H2O)${t} — cubic lattice a=${N} A, golden-angle orientations
${s.join(`
`)}`}const at=1.526,$=1.09,Ut=109.47*Math.PI/180;function Dt(t){const e=Ut/2,s=at*Math.sin(e),n=at*Math.cos(e),o=[];for(let c=0;c<t;c++)o.push([c*s,0,c%2===0?0:n]);const a=[];for(const c of o)a.push(_("C",c));for(let c=0;c<t;c++){const[f,l,m]=o[c],h=c%2===0?1:-1,y=m+h*$*Math.cos(e),g=$*Math.sin(e);a.push(_("H",[f,l+g,y])),a.push(_("H",[f,l-g,y]))}const r=o[0],i=o[t-1],u=[-Math.sin(e),0,t>1&&o[1][2]>r[2]?-Math.cos(e):Math.cos(e)],d=[Math.sin(e),0,t>1&&o[t-2][2]>i[2]?-Math.cos(e):Math.cos(e)];return a.push(_("H",[r[0]+$*u[0],r[1],r[2]+$*u[2]])),a.push(_("H",[i[0]+$*d[0],i[1],i[2]+$*d[2]])),`${t+2*t+2}
C${t}H${2*t+2} all-anti
${a.join(`
`)}`}function _(t,[e,s,n]){return`${t.padEnd(2)} ${e.toFixed(6).padStart(12)} ${s.toFixed(6).padStart(12)} ${n.toFixed(6).padStart(12)}`}const k=[{id:"water",label:"Water cluster (H₂O)ₙ",name:t=>`(H2O)${t}`,natoms:t=>3*t,xyz:Wt,ladder:[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32]},{id:"alkane",label:"n-Alkane CₙH₂ₙ₊₂",name:t=>`C${t}H${2*t+2}`,natoms:t=>3*t+2,xyz:Dt,ladder:[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20]}],Q="gansu-stress-v1",tt="gansu-stress-probe-v1",pt=180,V="6-31g(d,p)";function Jt(){try{const t=localStorage.getItem(Q);return t?JSON.parse(t):null}catch{return null}}function X(t){try{localStorage.setItem(Q,JSON.stringify(t))}catch{}}function qt(){try{localStorage.removeItem(Q)}catch{}}function et(t){return t==null?"—":t<1024**2?`${(t/1024).toFixed(0)} KB`:t<1024**3?`${(t/1024**2).toFixed(0)} MB`:`${(t/1024**3).toFixed(2)} GB`}function bt(t){return t==null?"—":(t/1e3).toFixed(2)}const st={ok:"ok",timeout:"timeout (exceeded budget)",error:"error / OOM","worker-died":"worker killed (OOM)","tab-crash":"tab crash"};function rt(t,e,s,n,o,a){return new Promise(r=>{const i=new Worker(new URL("/GANSU-Lite/assets/stressWorker-B8BED22c.js",import.meta.url),{type:"module"});a(i);let u=!1;const d=performance.now();let c=d,f="no progress reported";const l=g=>{u||(u=!0,clearTimeout(h),i.terminate(),a(null),r(g))},m={nbasis:null,totalMs:null,scfMs:null,iterations:null,converged:null,energy:null,eriBytes:null,heapUsedBytes:null,heapLimitBytes:null},h=setTimeout(()=>{const g=(performance.now()-c)/1e3;l({...m,status:"timeout",totalMs:performance.now()-d,detail:`exceeded ${(s/1e3).toFixed(0)} s budget; silent for ${g.toFixed(0)} s after "${f}"`+(g>30?" — a long silence usually means the worker was killed for memory":"")})},s);i.onmessage=g=>{const b=g.data;if(c=performance.now(),b.type==="progress"){f=b.message,o(b.message,b.elapsedMs);return}if(b.type==="done"){l({status:"ok",nbasis:b.nbasis,totalMs:b.totalMs,scfMs:b.scfMs,iterations:b.iterations,converged:b.converged,energy:b.energy,eriBytes:b.eriBytes,heapUsedBytes:b.heapUsedBytes,heapLimitBytes:b.heapLimitBytes,...b.converged?{}:{detail:"SCF did not converge"}});return}b.type==="error"&&l({...m,status:"error",nbasis:b.nbasis,totalMs:b.elapsedMs,detail:`${b.name}: ${b.message} (phase: ${b.phase})`})},i.onerror=g=>{l({...m,status:"worker-died",totalMs:performance.now()-d,detail:g.message||"worker terminated by the browser"})};const y={type:"stress-run",xyzText:t,basisGBS:e,charge:0,eriBackend:"auto",baseUrl:"/GANSU-Lite/",...n?{dftConfig:{functional:"B3LYP",gridLevel:"medium"}}:{}};i.postMessage(y)})}function Yt(){const t=new Set;for(const e of k[0].ladder)t.add(e*25);for(const e of[550,600,650,700])t.add(e);return[...t].sort((e,s)=>e-s)}function Vt(t){return new Promise(e=>{const s=new Worker(new URL("/GANSU-Lite/assets/stressWorker-B8BED22c.js",import.meta.url),{type:"module"}),n=[],o=()=>{s.terminate(),e(n)};s.onmessage=r=>{const i=r.data;if(i.type==="probe-result"){const u={nbasis:i.nbasis,bytes:i.bytes,ok:i.ok,ms:i.ms,gbPerSec:i.gbPerSec,verdict:i.verdict,...i.detail?{detail:i.detail}:{}};n.push(u),t(u);return}i.type==="probe-done"&&o()},s.onerror=()=>o();const a={type:"memory-probe",nbasisLadder:Yt()};s.postMessage(a)})}function ft(t,e){const s=e.verdict==="resident",n={resident:"resident in RAM",swapping:"SWAPPING — past usable memory",failed:"ALLOCATION FAILED",skipped:"skipped (past cap)"}[e.verdict],o=document.createElement("tr");o.innerHTML=`<td class="num">${e.nbasis}</td><td class="num">${et(e.bytes)}</td><td class="num">${e.ms.toFixed(0)}</td><td class="num">${e.gbPerSec>0?e.gbPerSec.toFixed(2):"—"}</td><td${s?"":' class="cb-err"'}>${n}${e.detail?`<br><span class="cb-err">${ht(e.detail)}</span>`:""}</td>`,t.appendChild(o)}let x=[];function Xt(t){try{localStorage.setItem(tt,JSON.stringify(t))}catch{}}function Zt(){try{const t=localStorage.getItem(tt);return t?JSON.parse(t):[]}catch{return[]}}async function Kt(){const t=document.getElementById("cb-stress-probe"),e=document.getElementById("cb-probe-body"),s=document.getElementById("cb-probe-section");t&&(t.disabled=!0),s&&(s.style.display=""),e&&(e.innerHTML=""),x=[],w("probing memory ceiling…"),await Vt(r=>{x.push(r),Xt(x),e&&ft(e,r)});const n=x.filter(r=>r.verdict==="resident"),o=x.find(r=>r.verdict!=="resident"),a=document.getElementById("cb-probe-summary");if(a){const r=n.length?n[n.length-1]:null;if(!r)a.textContent="Could not hold even the smallest probe array — something is wrong.";else{const i=o?{swapping:`at ${o.nbasis} the array no longer fits in real memory`,failed:`allocation failed outright at ${o.nbasis}`,skipped:`the ladder ran past the probe cap at ${o.nbasis}`,resident:""}[o.verdict]:"the ladder was exhausted without hitting a ceiling";a.textContent=`Memory ceiling: ${r.nbasis} basis functions — a ${et(r.bytes)} ERI array still resident in RAM (${r.gbPerSec.toFixed(1)} GB/s); ${i}.`}}w("memory probe done"),t&&(t.disabled=!1)}function Qt(t){const e=t.filter(o=>o.status==="ok"),s=e.length?e[e.length-1]:null,n=t.find(o=>o.status!=="ok")??null;return{largestCompleted:s&&s.nbasis!=null&&s.totalMs!=null?{name:s.name,natoms:s.natoms,nbasis:s.nbasis,seconds:+(s.totalMs/1e3).toFixed(1)}:null,firstFailure:n?{name:n.name,natoms:n.natoms,nbasis:n.nbasis,mode:st[n.status]}:null}}let p={series:k[0],points:[],running:!1,stopRequested:!1,activeWorker:null,output:null};function te(){return`
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
            ${k.map(t=>`<option value="${t.id}">${t.label}</option>`).join("")}
          </select>
        </label>
        <label>Theory
          <select id="cb-stress-theory">
            <option value="rhf">RHF/6-31G(d,p)</option>
            <option value="b3lyp">B3LYP/6-31G(d,p)</option>
          </select>
        </label>
        <label>Budget per point (s)
          <input id="cb-stress-budget" type="number" min="10" max="3600" step="10" value="${pt}">
        </label>
        <label>Ladder (n)
          <input id="cb-stress-ladder" type="text" value="${k[0].ladder.join(", ")}">
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
        The SCF's dominant allocation is the unique-ERI array (N⁴/8 doubles). Probing that
        allocation directly finds the memory ceiling in seconds, instead of waiting out an SCF that
        would take tens of minutes at the same size. Allocation success alone is not the test — a
        desktop OS will back a multi-gigabyte array with the page file and report success while
        thrashing — so the probe writes one value per page and stops when commit throughput
        collapses against the small-array baseline. Only sizes marked <em>resident in RAM</em> count.
      </p>
      <button id="cb-stress-probe">Probe memory ceiling</button>
      <div id="cb-probe-section" style="display:none;margin-top:10px">
        <table class="cb-results-table" style="font-size:0.8rem;max-width:520px">
          <thead><tr><th class="num">Basis fns</th><th class="num">ERI array</th><th class="num">ms</th><th class="num">GB/s</th><th>Result</th></tr></thead>
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
  `}function ee(t){const e=a=>document.getElementById(a),s=e("cb-stress-series"),n=e("cb-stress-ladder");s?.addEventListener("change",()=>{const a=k.find(r=>r.id===s.value);a&&n&&(p.series=a,n.value=a.ladder.join(", "))}),e("cb-stress-run")?.addEventListener("click",()=>{ae(t)}),e("cb-stress-stop")?.addEventListener("click",()=>{p.stopRequested=!0,p.activeWorker?.terminate(),w("stopping…")}),e("cb-stress-reset")?.addEventListener("click",()=>{qt();try{localStorage.removeItem(tt)}catch{}x=[],p.points=[],p.output=null,P();const a=e("cb-stress-recovered");a&&(a.style.display="none");const r=e("cb-stress-summary");r&&(r.style.display="none"),w("saved progress cleared")}),e("cb-stress-probe")?.addEventListener("click",()=>{Kt()}),e("cb-stress-copy-row")?.addEventListener("click",()=>it(yt())),e("cb-stress-copy-json")?.addEventListener("click",()=>it(JSON.stringify(p.output,null,2))),e("cb-stress-download-json")?.addEventListener("click",()=>ct("stress-test.json",JSON.stringify(p.output,null,2),"application/json")),e("cb-stress-download-csv")?.addEventListener("click",()=>ct("stress-test.csv",ce(p.output),"text/csv")),se();const o=Zt();if(o.length){x=o;const a=document.getElementById("cb-probe-section"),r=document.getElementById("cb-probe-body");if(a&&r){a.style.display="";for(const i of o)ft(r,i)}}}function se(){const t=Jt();if(!t)return;p.points=t.points??[];const e=k.find(s=>s.id===t.seriesId);if(e&&(p.series=e),t.inFlight){const s=t.inFlight;p.points.some(o=>o.n===s.n)||p.points.push({n:s.n,name:s.name,natoms:s.natoms,nbasis:null,status:"tab-crash",totalMs:null,scfMs:null,iterations:null,converged:null,energy:null,eriBytes:null,heapUsedBytes:null,heapLimitBytes:null,detail:"tab died while this point was running"});const n=document.getElementById("cb-stress-recovered");n&&(n.style.display="",n.textContent=`Recovered from a previous session: the tab died while running ${s.name}. Results up to that point are restored below and ${s.name} is recorded as a tab crash.`),X({...t,points:p.points,inFlight:null})}p.points.length&&P()}function w(t){const e=document.getElementById("cb-stress-status");e&&(e.textContent=t)}async function ne(){const t=await fetch(`/GANSU-Lite/basis/${V}.gbs`);if(!t.ok)throw new Error(`Cannot load basis ${V}: HTTP ${t.status}`);return t.text()}function oe(t,e){const s=t.split(/[,\s]+/).map(n=>parseInt(n,10)).filter(n=>Number.isFinite(n)&&n>0);return s.length?s:e}async function ae(t){if(p.running)return;const e=document.getElementById("cb-stress-run"),s=document.getElementById("cb-stress-stop"),n=Math.max(10,Number(document.getElementById("cb-stress-budget")?.value)||pt),o=document.getElementById("cb-stress-theory")?.value==="b3lyp",a=oe(document.getElementById("cb-stress-ladder")?.value??"",p.series.ladder);p.running=!0,p.stopRequested=!1,p.points=[],p.output=null,e&&(e.disabled=!0),s&&(s.disabled=!1),P();let r;try{w("loading basis set…"),r=await ne()}catch(c){w(`FAILED: ${c instanceof Error?c.message:String(c)}`),p.running=!1,e&&(e.disabled=!1),s&&(s.disabled=!0);return}const i=p.series,u=o?"B3LYP":"RHF",d=n*1e3;try{w("warm-up…"),await rt(i.xyz(a[0]),r,d,o,()=>{},c=>{p.activeWorker=c})}catch{}for(const c of a){if(p.stopRequested){w("stopped");break}const f=i.name(c),l=i.natoms(c);X({seriesId:i.id,method:u,budgetSeconds:n,points:p.points,inFlight:{n:c,name:f,natoms:l,startedAt:Date.now()}}),w(`${f} — starting…`);const m=await rt(i.xyz(c),r,d,o,(y,g)=>w(`${f} — ${y} (${(g/1e3).toFixed(1)} s)`),y=>{p.activeWorker=y}),h={n:c,name:f,natoms:l,...m};if(p.points.push(h),P(),X({seriesId:i.id,method:u,budgetSeconds:n,points:p.points,inFlight:null}),m.status!=="ok"){w(`${f} failed: ${st[m.status]} — stopping ladder`);break}w(`${f} ok — ${bt(m.totalMs)} s`),await new Promise(y=>setTimeout(y,200))}p.output={device:K(),backend:t(),seriesId:i.id,seriesLabel:i.label,method:u,basis:V,budgetSeconds:n,points:p.points,summary:Qt(p.points),memoryProbe:x,toolUrl:location.href},ie(),p.running=!1,p.activeWorker=null,e&&(e.disabled=!1),s&&(s.disabled=!0),p.stopRequested||w("done")}function P(){const t=document.getElementById("cb-stress-body");if(t){t.innerHTML="";for(const e of p.points){const s=document.createElement("tr"),n=e.status!=="ok";s.innerHTML=`
      <td>${e.name}</td>
      <td class="num">${e.natoms}</td>
      <td class="num">${e.nbasis??"—"}</td>
      <td class="num">${et(e.eriBytes)}</td>
      <td class="num">${e.iterations??"—"}</td>
      <td class="num">${bt(e.totalMs)}</td>
      <td${n?' class="cb-err"':""}>${st[e.status]}${e.detail?`<br><span class="cb-err">${ht(e.detail)}</span>`:""}</td>
    `,t.appendChild(s)}}}function ht(t){return t.replace(/[&<>"]/g,e=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"})[e])}function yt(){const t=p.output;if(!t)return"";const e=t.summary.largestCompleted,s=t.summary.firstFailure;return`| ${[re(t.device),e?e.name:"—",e?String(e.natoms):"—",e?String(e.nbasis):"—",e?e.seconds.toFixed(1):"—",s?s.name:"(no failure in ladder)",s?s.mode:"—"].join(" | ")} |`}function re(t){const e=t.ua;return/iPhone/.test(e)?"iPhone":/iPad/.test(e)?"iPad":/Mac OS X|Macintosh/.test(e)?"Apple-silicon laptop":/Android/.test(e)?"Android":/Windows/.test(e)?"Windows":t.platform}function ie(){const t=document.getElementById("cb-stress-summary"),e=document.getElementById("cb-stress-row");if(!t||!e||!p.output)return;t.style.display="";const s=p.output.summary,n="| Device | Largest completed | Atoms | Basis fns | Time (s) | First failing size | Failure mode |",o="|---|---|---|---|---|---|---|",a=[];s.firstFailure||a.push("# The ladder finished without a failure — extend it to find the limit.");const r=p.points.filter(i=>i.status==="ok"&&i.converged===!1);r.length&&a.push(`# SCF did not converge for: ${r.map(i=>i.name).join(", ")}`),e.textContent=[n,o,yt(),...a].join(`
`)}function it(t){navigator.clipboard.writeText(t).then(()=>alert("Copied to clipboard")).catch(()=>alert("Copy failed — open browser console"))}function ct(t,e,s){const n=new Blob([e],{type:s}),o=URL.createObjectURL(n),a=document.createElement("a");a.href=o,a.download=t,a.click(),setTimeout(()=>URL.revokeObjectURL(o),1e3)}function ce(t){if(!t)return"";const e=["n","system","atoms","nbasis","eri_bytes","iterations","converged","setup_plus_scf_ms","scf_ms","energy_hartree","status","detail"],s=t.points.map(o=>[o.n,o.name,o.natoms,o.nbasis??"",o.eriBytes??"",o.iterations??"",o.converged==null?"":String(o.converged),o.totalMs?.toFixed(2)??"",o.scfMs?.toFixed(2)??"",o.energy?.toFixed(8)??"",o.status,o.detail??""].map(a=>`"${String(a).replace(/"/g,'""')}"`).join(","));return[`# device.ua,"${t.device.ua.replace(/"/g,'""')}"`,`# device.platform,"${t.device.platform}"`,`# device.deviceMemoryGB,"${t.device.deviceMemoryGB??""}"`,`# device.hardwareConcurrency,"${t.device.hardwareConcurrency}"`,`# backend,"${t.backend}"`,`# series,"${t.seriesLabel}"`,`# method,"${t.method}/${t.basis}"`,`# budget_seconds,"${t.budgetSeconds}"`,`# timestamp,"${t.device.timestamp}"`].join(`
`)+`
`+e.join(",")+`
`+s.join(`
`)+`
`}const D=4,le=1,J=1.8897259886,gt=`3
H2O standard
O   0.000000   0.000000   0.117790
H   0.000000   0.755450  -0.471160
H   0.000000  -0.755450  -0.471160`,de=`3
H2O displaced
O   0.000000   0.000000   0.000000
H   0.000000   0.850000   0.620000
H   0.000000  -0.850000   0.620000`,ue=`3
CO2 linear
C   0.000000   0.000000   0.000000
O   0.000000   0.000000   1.160000
O   0.000000   0.000000  -1.160000`,q=new Map;async function I(t){if(q.has(t))return q.get(t);const e=`/GANSU-Lite/basis/${t}.gbs`,s=await fetch(e);if(!s.ok)throw new Error(`Failed to load basis ${t}: ${s.status}`);const n=Ft.fromGBS(await s.text());return q.set(t,n),n}async function me(t){const e=await I("6-31g(d,p)"),s=L(gt),n=new G(s,e,0);await W(n,"RHF").solve({eriBackend:t})}async function pe(t){const e=await I("sto-3g"),s=L(de);let n=new Float64Array(s.length*3);for(let r=0;r<s.length;r++)n[3*r]=s[r].coordinate.x,n[3*r+1]=s[r].coordinate.y,n[3*r+2]=s[r].coordinate.z;const o=.4,a=25;for(let r=0;r<a;r++){const i=s.map((l,m)=>({...l,coordinate:{x:n[3*m],y:n[3*m+1],z:n[3*m+2]}})),u=new G(i,e,0),d=W(u,"RHF");await d.solve({eriBackend:t});const c=Ot(u.primitiveShells,u.atoms,u.cgtoNormalizationFactors,u.numBasis,u.numAlphaSpins,d.density,d.coefficients,d.orbitalEnergies).total;let f=0;for(let l=0;l<c.length;l++)Math.abs(c[l])>f&&(f=Math.abs(c[l]));if(f<5e-4)break;for(let l=0;l<c.length;l++)n[l]-=o*c[l]}}async function be(t){const e=await I("sto-3g"),s=L(ue),n=new Float64Array(s.length*3);for(let o=0;o<s.length;o++)n[3*o]=s[o].coordinate.x*J,n[3*o+1]=s[o].coordinate.y*J,n[3*o+2]=s[o].coordinate.z*J;await jt(s.map(o=>o.atomicNumber),n,e,0,5e-4,void 0,void 0)}async function fe(t){const e=await I("6-31g(d,p)"),s=L(gt),n=new G(s,e,0);await W(n,"RHF",{functional:"B3LYP",gridLevel:"medium"}).solve({eriBackend:t})}const wt=[{id:"h2o_hf_631gdp",label:"H2O HF/6-31G(d,p)",run:me},{id:"h2o_geomopt_sto3g",label:"H2O geom opt HF/STO-3G",run:pe},{id:"co2_hess_sto3g",label:"CO2 Hessian HF/STO-3G",run:be},{id:"h2o_b3lyp_631gdp",label:"H2O B3LYP/6-31G(d,p)",run:fe}];function Z(t){const e=[...t].sort((n,o)=>n-o),s=e.length;return s===0?NaN:s%2===1?e[(s-1)/2]:(e[s/2-1]+e[s/2])/2}async function vt(t){const e=await fetch(t);if(!e.ok)throw new Error(`Cannot fetch ${t}`);const s=await e.arrayBuffer(),n=await WebAssembly.compile(s);let o;const a={"./wasm_eri_bg.js":{__wbindgen_init_externref_table:()=>{const l=o.exports.__wbindgen_externrefs,m=l.grow(4);l.set(0,void 0),l.set(m+0,void 0),l.set(m+1,null),l.set(m+2,!0),l.set(m+3,!1)}}};o=await WebAssembly.instantiate(n,a);const r=o.exports;r.__wbindgen_start&&r.__wbindgen_start();const i=r.memory,u=r.__wbindgen_malloc,d=r.__wbindgen_free,c=l=>{const m=u(l.length*8,8)>>>0;return new Float64Array(i.buffer).set(l,m/8),[m,l.length]},f=l=>{const[m,h]=l,y=new Float64Array(i.buffer).slice(m/8,m/8+h);return d(m,h*8,8),y};return{computeERIs:(l,m,h,y)=>{const[g,b]=c(l),[v,S]=c(m);return f(r.compute_eris_wasm(g,b,v,S,h,y))},computeFockRhf:(l,m,h,y)=>{const[g,b]=c(l),[v,S]=c(m),[E,O]=c(h);return f(r.compute_fock_rhf(g,b,v,S,E,O,y))}}}let F,C;async function he(){if(F!==void 0)return F;try{F=await vt("/GANSU-Lite/wasm/wasm_eri_bg.wasm")}catch{F=null}return F}async function ye(){if(C!==void 0)return C;try{C=await vt("/GANSU-Lite/wasm/wasm_eri_simd_bg.wasm")}catch{C=null}return C}async function Mt(){return await Et("/GANSU-Lite/"),_t()?Ct()==="wasm-simd"?{backend:"wasm",label:"WASM-SIMD"}:{backend:"wasm",label:"WASM"}:{backend:"js",label:"JS"}}function lt(t,e,s,n){let o=t,a=e;if(o<a){const c=o;o=a,a=c}let r=s,i=n;if(r<i){const c=r;r=i,i=c}let u=o*(o+1)/2+a,d=r*(r+1)/2+i;if(u<d){const c=u;u=d,d=c}return u*(u+1)/2+d}function ge(t,e,s,n){const o=new Float64Array(n*n);for(let a=0;a<n;a++)for(let r=a;r<n;r++){let i=0;for(let d=0;d<n;d++)for(let c=0;c<n;c++){const f=e[d*n+c];Math.abs(f)<1e-10||(i+=f*(t[lt(a,r,d,c)]-.5*t[lt(a,d,r,c)]))}const u=s[a*n+r]+i;o[a*n+r]=u,o[r*n+a]=u}return o}async function we(t,e,s,n){const o=t.numBasis,a=new Lt(e,o),r=W(t,"RHF");r.computeNuclearRepulsionEnergy(),r.computeCoreHamiltonianMatrix(),r.computeTransformMatrix(),r.eri=a,r.guessInitialFockMatrix();const i=r.coreHamiltonianMatrix.data;let u=0,d=0,c=0;const f=performance.now();for(let l=0;l<200;l++){r.computeCoefficientMatrix(),r.computeDensityMatrix();const m=r.densityMatrix.data,h=await n(e,m,i,o);if(r.fockMatrix=new It(o,o,h),c=r.computeEnergy(),l>0&&Math.abs(c-u)<1e-8){d=l+1;break}r.updateFockMatrix(),u=c}return{eriMs:s,scfMs:performance.now()-f,iters:d,energy:c}}const Y=[{molecule:"H2O",xyzFile:"H2O.xyz",basis:"6-31g(d,p)"},{molecule:"CO2",xyzFile:"CO2.xyz",basis:"6-31g(d,p)"},{molecule:"Benzene",xyzFile:"Benzene.xyz",basis:"6-31g(d,p)"}],dt=3;async function ve(t){const e=await I(t.basis),s=`/GANSU-Lite/xyz/${t.xyzFile}`,n=await fetch(s);if(!n.ok)throw new Error(`Cannot load ${s}`);const o=await n.text(),a=L(o);return{mol:new G(a,e,0),xyz:o}}async function Me(t,e,s,n){const{mol:o}=await ve(t),a=o.numBasis,r=o.primitiveShells,i=o.cgtoNormalizationFactors,u=(await Ht(async()=>{const{packShells:h}=await import("./ri-DT3jfd0v.js").then(y=>y.W);return{packShells:h}},[])).packShells(r),d=new Float64Array(i);async function c(h,y,g){const b=[],v=[];let S=0,E;try{for(let M=0;M<dt;M++){n(`  ${h} run ${M+1}/${dt}…`);const St=performance.now(),$t=await y(),kt=performance.now()-St,U=await we(o,$t,kt,g);b.push(U.eriMs),v.push(U.scfMs),M===0&&(S=U.iters),await new Promise(Bt=>setTimeout(Bt,30))}}catch(M){E=M instanceof Error?M.message:String(M)}const O=b.length?Z(b):NaN,nt=v.length?Z(v):NaN;return{eriMedianMs:O,scfMedianMs:nt,totalMedianMs:O+nt,iters:S,eriRuns:b,scfRuns:v,...E?{error:E}:{}}}const f=await c("JS",()=>Rt(r,i,a,1e-10,void 0,"js"),(h,y,g,b)=>ge(h,y,g,b)),l=e?await c("WASM baseline",()=>e.computeERIs(u,d,a,1e-10),(h,y,g,b)=>e.computeFockRhf(h,y,g,b)):{eriMedianMs:NaN,scfMedianMs:NaN,totalMedianMs:NaN,iters:0,eriRuns:[],scfRuns:[],error:"WASM baseline not available"},m=s?await c("WASM+SIMD",()=>s.computeERIs(u,d,a,1e-10),(h,y,g,b)=>s.computeFockRhf(h,y,g,b)):{eriMedianMs:NaN,scfMedianMs:NaN,totalMedianMs:NaN,iters:0,eriRuns:[],scfRuns:[],error:"SIMD not available"};return{molecule:t.molecule,basis:t.basis,nbasis:a,js:f,wasmBase:l,wasmSimd:m}}let T=null;async function xe(){const t=document.getElementById("cb-run-compare"),e=document.getElementById("cb-stop-compare");t&&(t.disabled=!0),e&&(e.disabled=!1);const s=r=>{const i=document.getElementById("cb-compare-status");i&&(i.textContent=r)};s("Loading WASM binaries…");const n=await he(),o=await ye(),a=[];for(let r=0;r<Y.length;r++){const i=Y[r];s(`[${r+1}/${Y.length}] ${i.molecule}/${i.basis}…`);try{a.push(await Me(i,n,o,s)),Se(a)}catch(u){const d=u instanceof Error?u.message:String(u);s(`  ERROR: ${d}`)}}T=a,s("Done."),t&&(t.disabled=!1),e&&(e.disabled=!0)}function Se(t){const e=document.getElementById("cb-compare-section");if(!e)return;e.style.display="";const s=document.getElementById("cb-compare-body");if(!s)return;s.innerHTML="";const n=a=>isFinite(a)?a.toFixed(0):"—",o=(a,r)=>isFinite(a)&&isFinite(r)&&r>0?`${(a/r).toFixed(1)}×`:"—";for(const a of t){const r=document.createElement("tr");r.innerHTML=`
      <td>${a.molecule}</td><td>${a.basis}</td><td class="num">${a.nbasis}</td><td class="num">${a.js.iters||"—"}</td>
      <td class="num">${n(a.js.eriMedianMs)}</td>
      <td class="num">${n(a.js.scfMedianMs)}</td>
      <td class="num"><b>${n(a.js.totalMedianMs)}</b></td>
      <td class="num">${n(a.wasmBase.eriMedianMs)}</td>
      <td class="num">${n(a.wasmBase.scfMedianMs)}</td>
      <td class="num"><b>${n(a.wasmBase.totalMedianMs)}</b></td>
      <td class="num">${o(a.js.totalMedianMs,a.wasmBase.totalMedianMs)}</td>
      <td class="num">${n(a.wasmSimd.eriMedianMs)}</td>
      <td class="num">${n(a.wasmSimd.scfMedianMs)}</td>
      <td class="num"><b>${n(a.wasmSimd.totalMedianMs)}</b></td>
      <td class="num">${o(a.js.totalMedianMs,a.wasmSimd.totalMedianMs)}</td>
    `,s.appendChild(r)}}const $e=document.getElementById("app");let B=null,ut="detecting…";function ke(){const t=K(),e=mt();$e.innerHTML=`
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
          <tr><th>Initial load</th><td>${e!=null?`${e} ms`:"unavailable"}</td></tr>
          <tr><th>Backend</th><td id="cb-backend">detecting…</td></tr>
        </table>
      </section>

      <section class="cb-panel">
        <h2>Workloads</h2>
        <p class="cb-note">All fixed: SCF tol = 1e-7, DIIS, default initial guess. Single Web Worker thread (no parallelism — single-core perf + memory bandwidth).</p>
        <ul class="cb-workloads">
          ${wt.map(s=>`<li id="cb-row-${s.id}"><span class="cb-wl-label">${s.label}</span><span class="cb-wl-status">pending</span></li>`).join("")}
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

      ${te()}

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
  `,Re(),document.getElementById("cb-run").addEventListener("click",Be),document.getElementById("cb-stop").addEventListener("click",()=>{z=!0}),document.getElementById("cb-run-compare")?.addEventListener("click",xe),document.getElementById("cb-copy-compare-json")?.addEventListener("click",()=>xt(JSON.stringify(T,null,2))),document.getElementById("cb-download-compare-json")?.addEventListener("click",()=>H("backend-comparison.json",JSON.stringify(T,null,2),"application/json")),document.getElementById("cb-download-compare-csv")?.addEventListener("click",()=>H("backend-comparison.csv",Fe(T),"text/csv")),document.getElementById("cb-copy-json")?.addEventListener("click",_e),document.getElementById("cb-download-json")?.addEventListener("click",()=>H("benchmark.json",JSON.stringify(B,null,2),"application/json")),document.getElementById("cb-download-csv")?.addEventListener("click",()=>H("benchmark.csv",Ce(B),"text/csv")),ee(()=>ut),Mt().then(({label:s})=>{ut=s;const n=document.getElementById("cb-backend");n&&(n.textContent=s)})}let z=!1;function R(t,e){const s=document.getElementById(`cb-row-${t}`);if(s){const n=s.querySelector(".cb-wl-status");n&&(n.textContent=e)}}async function Be(){z=!1;const t=document.getElementById("cb-run"),e=document.getElementById("cb-stop");t.disabled=!0,e.disabled=!1;const{backend:s,label:n}=await Mt(),o=K(),a=mt(),r=[];for(const i of wt){if(z){R(i.id,"stopped");continue}const u=[];let d;try{for(let l=0;l<D&&!z;l++){R(i.id,`run ${l+1}/${D}…`);const m=performance.now();await i.run(s);const h=performance.now()-m;u.push(h),R(i.id,`run ${l+1}/${D}: ${h.toFixed(0)} ms`),await new Promise(y=>setTimeout(y,50))}}catch(l){d=l instanceof Error?l.message:String(l),R(i.id,`ERROR: ${d}`)}const c=u.slice(le),f=c.length>0?Z(c):NaN;r.push({id:i.id,label:i.label,warmupMs:u[0]??NaN,runs:c,medianMs:f,...d?{error:d}:{}}),!d&&c.length>0&&R(i.id,`median ${f.toFixed(0)} ms`)}B={device:o,backend:n,initialLoadMs:a,workloads:r,toolUrl:location.href},Ee(B),t.disabled=!1,e.disabled=!0}function Ee(t){const e=document.getElementById("cb-results");e.style.display="";const s=document.getElementById("cb-results-body");s.innerHTML="";for(const n of t.workloads){const o=document.createElement("tr");o.innerHTML=`
      <td>${n.label}${n.error?` <span class="cb-err">(${n.error})</span>`:""}</td>
      <td class="num">${isFinite(n.warmupMs)?n.warmupMs.toFixed(0):"—"}</td>
      <td class="num">${n.runs[0]!=null?n.runs[0].toFixed(0):"—"}</td>
      <td class="num">${n.runs[1]!=null?n.runs[1].toFixed(0):"—"}</td>
      <td class="num">${n.runs[2]!=null?n.runs[2].toFixed(0):"—"}</td>
      <td class="num"><b>${isFinite(n.medianMs)?n.medianMs.toFixed(0):"—"}</b></td>
    `,s.appendChild(o)}document.getElementById("cb-json-pre").textContent=JSON.stringify(t,null,2)}function _e(){B&&xt(JSON.stringify(B,null,2))}function xt(t){navigator.clipboard.writeText(t).then(()=>alert("Copied to clipboard")).catch(()=>alert("Copy failed — open browser console"))}function Fe(t){if(!t||t.length===0)return"";const e=["molecule","basis","nbasis","iters","js_eri_ms","js_scf_ms","js_total_ms","wasm_eri_ms","wasm_scf_ms","wasm_total_ms","speedup_wasm_vs_js","simd_eri_ms","simd_scf_ms","simd_total_ms","speedup_simd_vs_js"],s=a=>isFinite(a)?a.toFixed(2):"",n=(a,r)=>isFinite(a)&&isFinite(r)&&r>0?(a/r).toFixed(2):"",o=t.map(a=>[a.molecule,a.basis,String(a.nbasis),String(a.js.iters),s(a.js.eriMedianMs),s(a.js.scfMedianMs),s(a.js.totalMedianMs),s(a.wasmBase.eriMedianMs),s(a.wasmBase.scfMedianMs),s(a.wasmBase.totalMedianMs),n(a.js.totalMedianMs,a.wasmBase.totalMedianMs),s(a.wasmSimd.eriMedianMs),s(a.wasmSimd.scfMedianMs),s(a.wasmSimd.totalMedianMs),n(a.js.totalMedianMs,a.wasmSimd.totalMedianMs)].map(r=>`"${r}"`).join(","));return e.join(",")+`
`+o.join(`
`)+`
`}function H(t,e,s){const n=new Blob([e],{type:s}),o=URL.createObjectURL(n),a=document.createElement("a");a.href=o,a.download=t,a.click(),setTimeout(()=>URL.revokeObjectURL(o),1e3)}function Ce(t){if(!t)return"";const e=["workload_id","workload_label","backend","warmup_ms","run1_ms","run2_ms","run3_ms","median_ms","error"],s=t.workloads.map(o=>[o.id,o.label,t.backend,isFinite(o.warmupMs)?o.warmupMs.toFixed(2):"",o.runs[0]?.toFixed(2)??"",o.runs[1]?.toFixed(2)??"",o.runs[2]?.toFixed(2)??"",isFinite(o.medianMs)?o.medianMs.toFixed(2):"",o.error??""].map(a=>`"${String(a).replace(/"/g,'""')}"`).join(","));return[`# device.ua,"${t.device.ua.replace(/"/g,'""')}"`,`# device.platform,"${t.device.platform}"`,`# device.deviceMemoryGB,"${t.device.deviceMemoryGB??""}"`,`# device.hardwareConcurrency,"${t.device.hardwareConcurrency}"`,`# device.screen,"${t.device.screen}"`,`# backend,"${t.backend}"`,`# initialLoadMs,"${t.initialLoadMs??""}"`,`# timestamp,"${t.device.timestamp}"`].join(`
`)+`
`+e.join(",")+`
`+s.join(`
`)+`
`}function Re(){if(document.getElementById("cb-styles"))return;const t=document.createElement("style");t.id="cb-styles",t.textContent=`
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
  `,document.head.appendChild(t)}ke();
