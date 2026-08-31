import"./styles-B8WF9G-F.js";import{p as R}from"./parseXYZ-ByfvMHmk.js";import{g as G,h as U,z as xt,o as St,f as $t,n as kt,C as Bt,F as Et,M as _t}from"./ri-DT3jfd0v.js";import{c as Ft}from"./gradient-Bmg7fDyY.js";import{c as Ct}from"./hessian-DpHI2yIv.js";import"./properties-Odo1T7fO.js";const Rt="modulepreload",Lt=function(t){return"/GANSU-Lite/"+t},et={},jt=function(e,s,n){let o=Promise.resolve();if(s&&s.length>0){let u=function(d){return Promise.all(d.map(c=>Promise.resolve(c).then(b=>({status:"fulfilled",value:b}),b=>({status:"rejected",reason:b}))))};document.getElementsByTagName("link");const r=document.querySelector("meta[property=csp-nonce]"),i=r?.nonce||r?.getAttribute("nonce");o=u(s.map(d=>{if(d=Lt(d),d in et)return;et[d]=!0;const c=d.endsWith(".css"),b=c?'[rel="stylesheet"]':"";if(document.querySelector(`link[href="${d}"]${b}`))return;const l=document.createElement("link");if(l.rel=c?"stylesheet":Rt,c||(l.as="script"),l.crossOrigin="",l.href=d,i&&l.setAttribute("nonce",i),document.head.appendChild(l),c)return new Promise((m,h)=>{l.addEventListener("load",m),l.addEventListener("error",()=>h(new Error(`Unable to preload CSS for ${d}`)))})}))}function a(r){const i=new Event("vite:preloadError",{cancelable:!0});if(i.payload=r,window.dispatchEvent(i),!i.defaultPrevented)throw r}return o.then(r=>{for(const i of r||[])i.status==="rejected"&&a(i.reason);return e().catch(a)})};function Z(){const t=navigator;return{ua:t.userAgent,platform:t.platform??"unknown",deviceMemoryGB:t.deviceMemory??null,hardwareConcurrency:t.hardwareConcurrency??0,screen:`${screen.width}×${screen.height}`,language:t.language||"unknown",timestamp:new Date().toISOString()}}function lt(){try{const t=performance.getEntriesByType("navigation")[0];return t?Math.round(t.loadEventEnd-t.fetchStart):null}catch{return null}}const I=.9572,O=104.52*Math.PI/180,A=3,It=[["O",0,0,0],["H",I*Math.sin(O/2),0,I*Math.cos(O/2)],["H",-I*Math.sin(O/2),0,I*Math.cos(O/2)]];function Ot(t){const e=Math.ceil(Math.cbrt(t))+1,s=[];for(let n=-e;n<=e;n++)for(let o=-e;o<=e;o++)for(let a=-e;a<=e;a++)s.push([n,o,a]);return s.sort((n,o)=>{const a=n[0]*n[0]+n[1]*n[1]+n[2]*n[2],r=o[0]*o[0]+o[1]*o[1]+o[2]*o[2];return a!==r?a-r:n[2]!==o[2]?n[2]-o[2]:n[1]!==o[1]?n[1]-o[1]:n[0]-o[0]}),s.slice(0,t)}function At(t){const e=Math.PI*(3-Math.sqrt(5));return[t*e,t*1.2345678,t*.7654321]}function Nt(t,[e,s,n]){const o=(r,i)=>[r[0]*Math.cos(i)-r[1]*Math.sin(i),r[0]*Math.sin(i)+r[1]*Math.cos(i),r[2]];return o(((r,i)=>[r[0]*Math.cos(i)+r[2]*Math.sin(i),r[1],-r[0]*Math.sin(i)+r[2]*Math.cos(i)])(o([t[0],t[1],t[2]],n),s),e)}function Ht(t){const e=Ot(t),s=[];for(let n=0;n<t;n++){const[o,a,r]=e[n],i=At(n),u=o*A,d=a*A,c=r*A;for(const[b,l,m,h]of It){const[g,y,f]=Nt([l,m,h],i);s.push(`${b.padEnd(2)} ${(u+g).toFixed(6).padStart(12)} ${(d+y).toFixed(6).padStart(12)} ${(c+f).toFixed(6).padStart(12)}`)}}return`${t*3}
(H2O)${t} — cubic lattice a=${A} A, golden-angle orientations
${s.join(`
`)}`}const st=1.526,S=1.09,Tt=109.47*Math.PI/180;function zt(t){const e=Tt/2,s=st*Math.sin(e),n=st*Math.cos(e),o=[];for(let c=0;c<t;c++)o.push([c*s,0,c%2===0?0:n]);const a=[];for(const c of o)a.push(B("C",c));for(let c=0;c<t;c++){const[b,l,m]=o[c],h=c%2===0?1:-1,g=m+h*S*Math.cos(e),y=S*Math.sin(e);a.push(B("H",[b,l+y,g])),a.push(B("H",[b,l-y,g]))}const r=o[0],i=o[t-1],u=[-Math.sin(e),0,t>1&&o[1][2]>r[2]?-Math.cos(e):Math.cos(e)],d=[Math.sin(e),0,t>1&&o[t-2][2]>i[2]?-Math.cos(e):Math.cos(e)];return a.push(B("H",[r[0]+S*u[0],r[1],r[2]+S*u[2]])),a.push(B("H",[i[0]+S*d[0],i[1],i[2]+S*d[2]])),`${t+2*t+2}
C${t}H${2*t+2} all-anti
${a.join(`
`)}`}function B(t,[e,s,n]){return`${t.padEnd(2)} ${e.toFixed(6).padStart(12)} ${s.toFixed(6).padStart(12)} ${n.toFixed(6).padStart(12)}`}const C=[{id:"water",label:"Water cluster (H₂O)ₙ",name:t=>`(H2O)${t}`,natoms:t=>3*t,xyz:Ht,ladder:[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32]},{id:"alkane",label:"n-Alkane CₙH₂ₙ₊₂",name:t=>`C${t}H${2*t+2}`,natoms:t=>3*t+2,xyz:zt,ladder:[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20]}],K="gansu-stress-v1",dt=180,Y="6-31g(d,p)";function Gt(){try{const t=localStorage.getItem(K);return t?JSON.parse(t):null}catch{return null}}function V(t){try{localStorage.setItem(K,JSON.stringify(t))}catch{}}function Ut(){try{localStorage.removeItem(K)}catch{}}function Wt(t){return t==null?"—":t<1024**2?`${(t/1024).toFixed(0)} KB`:t<1024**3?`${(t/1024**2).toFixed(0)} MB`:`${(t/1024**3).toFixed(2)} GB`}function ut(t){return t==null?"—":(t/1e3).toFixed(2)}const Q={ok:"ok",timeout:"timeout (exceeded budget)",error:"error / OOM","worker-died":"worker killed (OOM)","tab-crash":"tab crash"};function nt(t,e,s,n,o,a){return new Promise(r=>{const i=new Worker(new URL("/GANSU-Lite/assets/stressWorker-ChpibzU7.js",import.meta.url),{type:"module"});a(i);let u=!1;const d=performance.now();let c=d,b="no progress reported";const l=y=>{u||(u=!0,clearTimeout(h),i.terminate(),a(null),r(y))},m={nbasis:null,totalMs:null,scfMs:null,iterations:null,converged:null,energy:null,eriBytes:null,heapUsedBytes:null,heapLimitBytes:null},h=setTimeout(()=>{const y=(performance.now()-c)/1e3;l({...m,status:"timeout",totalMs:performance.now()-d,detail:`exceeded ${(s/1e3).toFixed(0)} s budget; silent for ${y.toFixed(0)} s after "${b}"`+(y>30?" — a long silence usually means the worker was killed for memory":"")})},s);i.onmessage=y=>{const f=y.data;if(c=performance.now(),f.type==="progress"){b=f.message,o(f.message,f.elapsedMs);return}if(f.type==="done"){l({status:"ok",nbasis:f.nbasis,totalMs:f.totalMs,scfMs:f.scfMs,iterations:f.iterations,converged:f.converged,energy:f.energy,eriBytes:f.eriBytes,heapUsedBytes:f.heapUsedBytes,heapLimitBytes:f.heapLimitBytes,...f.converged?{}:{detail:"SCF did not converge"}});return}l({...m,status:"error",nbasis:f.nbasis,totalMs:f.elapsedMs,detail:`${f.name}: ${f.message} (phase: ${f.phase})`})},i.onerror=y=>{l({...m,status:"worker-died",totalMs:performance.now()-d,detail:y.message||"worker terminated by the browser"})};const g={type:"stress-run",xyzText:t,basisGBS:e,charge:0,eriBackend:"auto",baseUrl:"/GANSU-Lite/",...n?{dftConfig:{functional:"B3LYP",gridLevel:"medium"}}:{}};i.postMessage(g)})}function Pt(t){const e=t.filter(o=>o.status==="ok"),s=e.length?e[e.length-1]:null,n=t.find(o=>o.status!=="ok")??null;return{largestCompleted:s&&s.nbasis!=null&&s.totalMs!=null?{name:s.name,natoms:s.natoms,nbasis:s.nbasis,seconds:+(s.totalMs/1e3).toFixed(1)}:null,firstFailure:n?{name:n.name,natoms:n.natoms,nbasis:n.nbasis,mode:Q[n.status]}:null}}let p={series:C[0],points:[],running:!1,stopRequested:!1,activeWorker:null,output:null};function Dt(){return`
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
            ${C.map(t=>`<option value="${t.id}">${t.label}</option>`).join("")}
          </select>
        </label>
        <label>Theory
          <select id="cb-stress-theory">
            <option value="rhf">RHF/6-31G(d,p)</option>
            <option value="b3lyp">B3LYP/6-31G(d,p)</option>
          </select>
        </label>
        <label>Budget per point (s)
          <input id="cb-stress-budget" type="number" min="10" max="3600" step="10" value="${dt}">
        </label>
        <label>Ladder (n)
          <input id="cb-stress-ladder" type="text" value="${C[0].ladder.join(", ")}">
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
      <div id="cb-stress-summary" style="display:none">
        <h3 style="font-size:0.95rem;margin:14px 0 6px">Row for the paper table</h3>
        <pre id="cb-stress-row" class="cb-stress-row"></pre>
        <button id="cb-stress-copy-row">Copy row</button>
        <button id="cb-stress-copy-json">Copy JSON</button>
        <button id="cb-stress-download-json">Download JSON</button>
        <button id="cb-stress-download-csv">Download CSV</button>
      </div>
    </section>
  `}function Jt(t){const e=o=>document.getElementById(o),s=e("cb-stress-series"),n=e("cb-stress-ladder");s?.addEventListener("change",()=>{const o=C.find(a=>a.id===s.value);o&&n&&(p.series=o,n.value=o.ladder.join(", "))}),e("cb-stress-run")?.addEventListener("click",()=>{Xt(t)}),e("cb-stress-stop")?.addEventListener("click",()=>{p.stopRequested=!0,p.activeWorker?.terminate(),w("stopping…")}),e("cb-stress-reset")?.addEventListener("click",()=>{Ut(),p.points=[],p.output=null,z();const o=e("cb-stress-recovered");o&&(o.style.display="none");const a=e("cb-stress-summary");a&&(a.style.display="none"),w("saved progress cleared")}),e("cb-stress-copy-row")?.addEventListener("click",()=>ot(mt())),e("cb-stress-copy-json")?.addEventListener("click",()=>ot(JSON.stringify(p.output,null,2))),e("cb-stress-download-json")?.addEventListener("click",()=>at("stress-test.json",JSON.stringify(p.output,null,2),"application/json")),e("cb-stress-download-csv")?.addEventListener("click",()=>at("stress-test.csv",te(p.output),"text/csv")),qt()}function qt(){const t=Gt();if(!t)return;p.points=t.points??[];const e=C.find(s=>s.id===t.seriesId);if(e&&(p.series=e),t.inFlight){const s=t.inFlight;p.points.some(o=>o.n===s.n)||p.points.push({n:s.n,name:s.name,natoms:s.natoms,nbasis:null,status:"tab-crash",totalMs:null,scfMs:null,iterations:null,converged:null,energy:null,eriBytes:null,heapUsedBytes:null,heapLimitBytes:null,detail:"tab died while this point was running"});const n=document.getElementById("cb-stress-recovered");n&&(n.style.display="",n.textContent=`Recovered from a previous session: the tab died while running ${s.name}. Results up to that point are restored below and ${s.name} is recorded as a tab crash.`),V({...t,points:p.points,inFlight:null})}p.points.length&&z()}function w(t){const e=document.getElementById("cb-stress-status");e&&(e.textContent=t)}async function Yt(){const t=await fetch(`/GANSU-Lite/basis/${Y}.gbs`);if(!t.ok)throw new Error(`Cannot load basis ${Y}: HTTP ${t.status}`);return t.text()}function Vt(t,e){const s=t.split(/[,\s]+/).map(n=>parseInt(n,10)).filter(n=>Number.isFinite(n)&&n>0);return s.length?s:e}async function Xt(t){if(p.running)return;const e=document.getElementById("cb-stress-run"),s=document.getElementById("cb-stress-stop"),n=Math.max(10,Number(document.getElementById("cb-stress-budget")?.value)||dt),o=document.getElementById("cb-stress-theory")?.value==="b3lyp",a=Vt(document.getElementById("cb-stress-ladder")?.value??"",p.series.ladder);p.running=!0,p.stopRequested=!1,p.points=[],p.output=null,e&&(e.disabled=!0),s&&(s.disabled=!1),z();let r;try{w("loading basis set…"),r=await Yt()}catch(c){w(`FAILED: ${c instanceof Error?c.message:String(c)}`),p.running=!1,e&&(e.disabled=!1),s&&(s.disabled=!0);return}const i=p.series,u=o?"B3LYP":"RHF",d=n*1e3;try{w("warm-up…"),await nt(i.xyz(a[0]),r,d,o,()=>{},c=>{p.activeWorker=c})}catch{}for(const c of a){if(p.stopRequested){w("stopped");break}const b=i.name(c),l=i.natoms(c);V({seriesId:i.id,method:u,budgetSeconds:n,points:p.points,inFlight:{n:c,name:b,natoms:l,startedAt:Date.now()}}),w(`${b} — starting…`);const m=await nt(i.xyz(c),r,d,o,(g,y)=>w(`${b} — ${g} (${(y/1e3).toFixed(1)} s)`),g=>{p.activeWorker=g}),h={n:c,name:b,natoms:l,...m};if(p.points.push(h),z(),V({seriesId:i.id,method:u,budgetSeconds:n,points:p.points,inFlight:null}),m.status!=="ok"){w(`${b} failed: ${Q[m.status]} — stopping ladder`);break}w(`${b} ok — ${ut(m.totalMs)} s`),await new Promise(g=>setTimeout(g,200))}p.output={device:Z(),backend:t(),seriesId:i.id,seriesLabel:i.label,method:u,basis:Y,budgetSeconds:n,points:p.points,summary:Pt(p.points),toolUrl:location.href},Qt(),p.running=!1,p.activeWorker=null,e&&(e.disabled=!1),s&&(s.disabled=!0),p.stopRequested||w("done")}function z(){const t=document.getElementById("cb-stress-body");if(t){t.innerHTML="";for(const e of p.points){const s=document.createElement("tr"),n=e.status!=="ok";s.innerHTML=`
      <td>${e.name}</td>
      <td class="num">${e.natoms}</td>
      <td class="num">${e.nbasis??"—"}</td>
      <td class="num">${Wt(e.eriBytes)}</td>
      <td class="num">${e.iterations??"—"}</td>
      <td class="num">${ut(e.totalMs)}</td>
      <td${n?' class="cb-err"':""}>${Q[e.status]}${e.detail?`<br><span class="cb-err">${Zt(e.detail)}</span>`:""}</td>
    `,t.appendChild(s)}}}function Zt(t){return t.replace(/[&<>"]/g,e=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"})[e])}function mt(){const t=p.output;if(!t)return"";const e=t.summary.largestCompleted,s=t.summary.firstFailure;return`| ${[Kt(t.device),e?e.name:"—",e?String(e.natoms):"—",e?String(e.nbasis):"—",e?e.seconds.toFixed(1):"—",s?s.name:"(no failure in ladder)",s?s.mode:"—"].join(" | ")} |`}function Kt(t){const e=t.ua;return/iPhone/.test(e)?"iPhone":/iPad/.test(e)?"iPad":/Mac OS X|Macintosh/.test(e)?"Apple-silicon laptop":/Android/.test(e)?"Android":/Windows/.test(e)?"Windows":t.platform}function Qt(){const t=document.getElementById("cb-stress-summary"),e=document.getElementById("cb-stress-row");if(!t||!e||!p.output)return;t.style.display="";const s=p.output.summary,n="| Device | Largest completed | Atoms | Basis fns | Time (s) | First failing size | Failure mode |",o="|---|---|---|---|---|---|---|",a=[];s.firstFailure||a.push("# The ladder finished without a failure — extend it to find the limit.");const r=p.points.filter(i=>i.status==="ok"&&i.converged===!1);r.length&&a.push(`# SCF did not converge for: ${r.map(i=>i.name).join(", ")}`),e.textContent=[n,o,mt(),...a].join(`
`)}function ot(t){navigator.clipboard.writeText(t).then(()=>alert("Copied to clipboard")).catch(()=>alert("Copy failed — open browser console"))}function at(t,e,s){const n=new Blob([e],{type:s}),o=URL.createObjectURL(n),a=document.createElement("a");a.href=o,a.download=t,a.click(),setTimeout(()=>URL.revokeObjectURL(o),1e3)}function te(t){if(!t)return"";const e=["n","system","atoms","nbasis","eri_bytes","iterations","converged","setup_plus_scf_ms","scf_ms","energy_hartree","status","detail"],s=t.points.map(o=>[o.n,o.name,o.natoms,o.nbasis??"",o.eriBytes??"",o.iterations??"",o.converged==null?"":String(o.converged),o.totalMs?.toFixed(2)??"",o.scfMs?.toFixed(2)??"",o.energy?.toFixed(8)??"",o.status,o.detail??""].map(a=>`"${String(a).replace(/"/g,'""')}"`).join(","));return[`# device.ua,"${t.device.ua.replace(/"/g,'""')}"`,`# device.platform,"${t.device.platform}"`,`# device.deviceMemoryGB,"${t.device.deviceMemoryGB??""}"`,`# device.hardwareConcurrency,"${t.device.hardwareConcurrency}"`,`# backend,"${t.backend}"`,`# series,"${t.seriesLabel}"`,`# method,"${t.method}/${t.basis}"`,`# budget_seconds,"${t.budgetSeconds}"`,`# timestamp,"${t.device.timestamp}"`].join(`
`)+`
`+e.join(",")+`
`+s.join(`
`)+`
`}const P=4,ee=1,D=1.8897259886,pt=`3
H2O standard
O   0.000000   0.000000   0.117790
H   0.000000   0.755450  -0.471160
H   0.000000  -0.755450  -0.471160`,se=`3
H2O displaced
O   0.000000   0.000000   0.000000
H   0.000000   0.850000   0.620000
H   0.000000  -0.850000   0.620000`,ne=`3
CO2 linear
C   0.000000   0.000000   0.000000
O   0.000000   0.000000   1.160000
O   0.000000   0.000000  -1.160000`,J=new Map;async function L(t){if(J.has(t))return J.get(t);const e=`/GANSU-Lite/basis/${t}.gbs`,s=await fetch(e);if(!s.ok)throw new Error(`Failed to load basis ${t}: ${s.status}`);const n=$t.fromGBS(await s.text());return J.set(t,n),n}async function oe(t){const e=await L("6-31g(d,p)"),s=R(pt),n=new G(s,e,0);await U(n,"RHF").solve({eriBackend:t})}async function ae(t){const e=await L("sto-3g"),s=R(se);let n=new Float64Array(s.length*3);for(let r=0;r<s.length;r++)n[3*r]=s[r].coordinate.x,n[3*r+1]=s[r].coordinate.y,n[3*r+2]=s[r].coordinate.z;const o=.4,a=25;for(let r=0;r<a;r++){const i=s.map((l,m)=>({...l,coordinate:{x:n[3*m],y:n[3*m+1],z:n[3*m+2]}})),u=new G(i,e,0),d=U(u,"RHF");await d.solve({eriBackend:t});const c=Ft(u.primitiveShells,u.atoms,u.cgtoNormalizationFactors,u.numBasis,u.numAlphaSpins,d.density,d.coefficients,d.orbitalEnergies).total;let b=0;for(let l=0;l<c.length;l++)Math.abs(c[l])>b&&(b=Math.abs(c[l]));if(b<5e-4)break;for(let l=0;l<c.length;l++)n[l]-=o*c[l]}}async function re(t){const e=await L("sto-3g"),s=R(ne),n=new Float64Array(s.length*3);for(let o=0;o<s.length;o++)n[3*o]=s[o].coordinate.x*D,n[3*o+1]=s[o].coordinate.y*D,n[3*o+2]=s[o].coordinate.z*D;await Ct(s.map(o=>o.atomicNumber),n,e,0,5e-4,void 0,void 0)}async function ie(t){const e=await L("6-31g(d,p)"),s=R(pt),n=new G(s,e,0);await U(n,"RHF",{functional:"B3LYP",gridLevel:"medium"}).solve({eriBackend:t})}const bt=[{id:"h2o_hf_631gdp",label:"H2O HF/6-31G(d,p)",run:oe},{id:"h2o_geomopt_sto3g",label:"H2O geom opt HF/STO-3G",run:ae},{id:"co2_hess_sto3g",label:"CO2 Hessian HF/STO-3G",run:re},{id:"h2o_b3lyp_631gdp",label:"H2O B3LYP/6-31G(d,p)",run:ie}];function X(t){const e=[...t].sort((n,o)=>n-o),s=e.length;return s===0?NaN:s%2===1?e[(s-1)/2]:(e[s/2-1]+e[s/2])/2}async function ft(t){const e=await fetch(t);if(!e.ok)throw new Error(`Cannot fetch ${t}`);const s=await e.arrayBuffer(),n=await WebAssembly.compile(s);let o;const a={"./wasm_eri_bg.js":{__wbindgen_init_externref_table:()=>{const l=o.exports.__wbindgen_externrefs,m=l.grow(4);l.set(0,void 0),l.set(m+0,void 0),l.set(m+1,null),l.set(m+2,!0),l.set(m+3,!1)}}};o=await WebAssembly.instantiate(n,a);const r=o.exports;r.__wbindgen_start&&r.__wbindgen_start();const i=r.memory,u=r.__wbindgen_malloc,d=r.__wbindgen_free,c=l=>{const m=u(l.length*8,8)>>>0;return new Float64Array(i.buffer).set(l,m/8),[m,l.length]},b=l=>{const[m,h]=l,g=new Float64Array(i.buffer).slice(m/8,m/8+h);return d(m,h*8,8),g};return{computeERIs:(l,m,h,g)=>{const[y,f]=c(l),[v,x]=c(m);return b(r.compute_eris_wasm(y,f,v,x,h,g))},computeFockRhf:(l,m,h,g)=>{const[y,f]=c(l),[v,x]=c(m),[k,j]=c(h);return b(r.compute_fock_rhf(y,f,v,x,k,j,g))}}}let E,_;async function ce(){if(E!==void 0)return E;try{E=await ft("/GANSU-Lite/wasm/wasm_eri_bg.wasm")}catch{E=null}return E}async function le(){if(_!==void 0)return _;try{_=await ft("/GANSU-Lite/wasm/wasm_eri_simd_bg.wasm")}catch{_=null}return _}async function ht(){return await xt("/GANSU-Lite/"),St()?kt()==="wasm-simd"?{backend:"wasm",label:"WASM-SIMD"}:{backend:"wasm",label:"WASM"}:{backend:"js",label:"JS"}}function rt(t,e,s,n){let o=t,a=e;if(o<a){const c=o;o=a,a=c}let r=s,i=n;if(r<i){const c=r;r=i,i=c}let u=o*(o+1)/2+a,d=r*(r+1)/2+i;if(u<d){const c=u;u=d,d=c}return u*(u+1)/2+d}function de(t,e,s,n){const o=new Float64Array(n*n);for(let a=0;a<n;a++)for(let r=a;r<n;r++){let i=0;for(let d=0;d<n;d++)for(let c=0;c<n;c++){const b=e[d*n+c];Math.abs(b)<1e-10||(i+=b*(t[rt(a,r,d,c)]-.5*t[rt(a,d,r,c)]))}const u=s[a*n+r]+i;o[a*n+r]=u,o[r*n+a]=u}return o}async function ue(t,e,s,n){const o=t.numBasis,a=new Et(e,o),r=U(t,"RHF");r.computeNuclearRepulsionEnergy(),r.computeCoreHamiltonianMatrix(),r.computeTransformMatrix(),r.eri=a,r.guessInitialFockMatrix();const i=r.coreHamiltonianMatrix.data;let u=0,d=0,c=0;const b=performance.now();for(let l=0;l<200;l++){r.computeCoefficientMatrix(),r.computeDensityMatrix();const m=r.densityMatrix.data,h=await n(e,m,i,o);if(r.fockMatrix=new _t(o,o,h),c=r.computeEnergy(),l>0&&Math.abs(c-u)<1e-8){d=l+1;break}r.updateFockMatrix(),u=c}return{eriMs:s,scfMs:performance.now()-b,iters:d,energy:c}}const q=[{molecule:"H2O",xyzFile:"H2O.xyz",basis:"6-31g(d,p)"},{molecule:"CO2",xyzFile:"CO2.xyz",basis:"6-31g(d,p)"},{molecule:"Benzene",xyzFile:"Benzene.xyz",basis:"6-31g(d,p)"}],it=3;async function me(t){const e=await L(t.basis),s=`/GANSU-Lite/xyz/${t.xyzFile}`,n=await fetch(s);if(!n.ok)throw new Error(`Cannot load ${s}`);const o=await n.text(),a=R(o);return{mol:new G(a,e,0),xyz:o}}async function pe(t,e,s,n){const{mol:o}=await me(t),a=o.numBasis,r=o.primitiveShells,i=o.cgtoNormalizationFactors,u=(await jt(async()=>{const{packShells:h}=await import("./ri-DT3jfd0v.js").then(g=>g.W);return{packShells:h}},[])).packShells(r),d=new Float64Array(i);async function c(h,g,y){const f=[],v=[];let x=0,k;try{for(let M=0;M<it;M++){n(`  ${h} run ${M+1}/${it}…`);const yt=performance.now(),wt=await g(),vt=performance.now()-yt,W=await ue(o,wt,vt,y);f.push(W.eriMs),v.push(W.scfMs),M===0&&(x=W.iters),await new Promise(Mt=>setTimeout(Mt,30))}}catch(M){k=M instanceof Error?M.message:String(M)}const j=f.length?X(f):NaN,tt=v.length?X(v):NaN;return{eriMedianMs:j,scfMedianMs:tt,totalMedianMs:j+tt,iters:x,eriRuns:f,scfRuns:v,...k?{error:k}:{}}}const b=await c("JS",()=>Bt(r,i,a,1e-10,void 0,"js"),(h,g,y,f)=>de(h,g,y,f)),l=e?await c("WASM baseline",()=>e.computeERIs(u,d,a,1e-10),(h,g,y,f)=>e.computeFockRhf(h,g,y,f)):{eriMedianMs:NaN,scfMedianMs:NaN,totalMedianMs:NaN,iters:0,eriRuns:[],scfRuns:[],error:"WASM baseline not available"},m=s?await c("WASM+SIMD",()=>s.computeERIs(u,d,a,1e-10),(h,g,y,f)=>s.computeFockRhf(h,g,y,f)):{eriMedianMs:NaN,scfMedianMs:NaN,totalMedianMs:NaN,iters:0,eriRuns:[],scfRuns:[],error:"SIMD not available"};return{molecule:t.molecule,basis:t.basis,nbasis:a,js:b,wasmBase:l,wasmSimd:m}}let H=null;async function be(){const t=document.getElementById("cb-run-compare"),e=document.getElementById("cb-stop-compare");t&&(t.disabled=!0),e&&(e.disabled=!1);const s=r=>{const i=document.getElementById("cb-compare-status");i&&(i.textContent=r)};s("Loading WASM binaries…");const n=await ce(),o=await le(),a=[];for(let r=0;r<q.length;r++){const i=q[r];s(`[${r+1}/${q.length}] ${i.molecule}/${i.basis}…`);try{a.push(await pe(i,n,o,s)),fe(a)}catch(u){const d=u instanceof Error?u.message:String(u);s(`  ERROR: ${d}`)}}H=a,s("Done."),t&&(t.disabled=!1),e&&(e.disabled=!0)}function fe(t){const e=document.getElementById("cb-compare-section");if(!e)return;e.style.display="";const s=document.getElementById("cb-compare-body");if(!s)return;s.innerHTML="";const n=a=>isFinite(a)?a.toFixed(0):"—",o=(a,r)=>isFinite(a)&&isFinite(r)&&r>0?`${(a/r).toFixed(1)}×`:"—";for(const a of t){const r=document.createElement("tr");r.innerHTML=`
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
    `,s.appendChild(r)}}const he=document.getElementById("app");let $=null,ct="detecting…";function ge(){const t=Z(),e=lt();he.innerHTML=`
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
          ${bt.map(s=>`<li id="cb-row-${s.id}"><span class="cb-wl-label">${s.label}</span><span class="cb-wl-status">pending</span></li>`).join("")}
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

      ${Dt()}

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
  `,Se(),document.getElementById("cb-run").addEventListener("click",ye),document.getElementById("cb-stop").addEventListener("click",()=>{T=!0}),document.getElementById("cb-run-compare")?.addEventListener("click",be),document.getElementById("cb-copy-compare-json")?.addEventListener("click",()=>gt(JSON.stringify(H,null,2))),document.getElementById("cb-download-compare-json")?.addEventListener("click",()=>N("backend-comparison.json",JSON.stringify(H,null,2),"application/json")),document.getElementById("cb-download-compare-csv")?.addEventListener("click",()=>N("backend-comparison.csv",Me(H),"text/csv")),document.getElementById("cb-copy-json")?.addEventListener("click",ve),document.getElementById("cb-download-json")?.addEventListener("click",()=>N("benchmark.json",JSON.stringify($,null,2),"application/json")),document.getElementById("cb-download-csv")?.addEventListener("click",()=>N("benchmark.csv",xe($),"text/csv")),Jt(()=>ct),ht().then(({label:s})=>{ct=s;const n=document.getElementById("cb-backend");n&&(n.textContent=s)})}let T=!1;function F(t,e){const s=document.getElementById(`cb-row-${t}`);if(s){const n=s.querySelector(".cb-wl-status");n&&(n.textContent=e)}}async function ye(){T=!1;const t=document.getElementById("cb-run"),e=document.getElementById("cb-stop");t.disabled=!0,e.disabled=!1;const{backend:s,label:n}=await ht(),o=Z(),a=lt(),r=[];for(const i of bt){if(T){F(i.id,"stopped");continue}const u=[];let d;try{for(let l=0;l<P&&!T;l++){F(i.id,`run ${l+1}/${P}…`);const m=performance.now();await i.run(s);const h=performance.now()-m;u.push(h),F(i.id,`run ${l+1}/${P}: ${h.toFixed(0)} ms`),await new Promise(g=>setTimeout(g,50))}}catch(l){d=l instanceof Error?l.message:String(l),F(i.id,`ERROR: ${d}`)}const c=u.slice(ee),b=c.length>0?X(c):NaN;r.push({id:i.id,label:i.label,warmupMs:u[0]??NaN,runs:c,medianMs:b,...d?{error:d}:{}}),!d&&c.length>0&&F(i.id,`median ${b.toFixed(0)} ms`)}$={device:o,backend:n,initialLoadMs:a,workloads:r,toolUrl:location.href},we($),t.disabled=!1,e.disabled=!0}function we(t){const e=document.getElementById("cb-results");e.style.display="";const s=document.getElementById("cb-results-body");s.innerHTML="";for(const n of t.workloads){const o=document.createElement("tr");o.innerHTML=`
      <td>${n.label}${n.error?` <span class="cb-err">(${n.error})</span>`:""}</td>
      <td class="num">${isFinite(n.warmupMs)?n.warmupMs.toFixed(0):"—"}</td>
      <td class="num">${n.runs[0]!=null?n.runs[0].toFixed(0):"—"}</td>
      <td class="num">${n.runs[1]!=null?n.runs[1].toFixed(0):"—"}</td>
      <td class="num">${n.runs[2]!=null?n.runs[2].toFixed(0):"—"}</td>
      <td class="num"><b>${isFinite(n.medianMs)?n.medianMs.toFixed(0):"—"}</b></td>
    `,s.appendChild(o)}document.getElementById("cb-json-pre").textContent=JSON.stringify(t,null,2)}function ve(){$&&gt(JSON.stringify($,null,2))}function gt(t){navigator.clipboard.writeText(t).then(()=>alert("Copied to clipboard")).catch(()=>alert("Copy failed — open browser console"))}function Me(t){if(!t||t.length===0)return"";const e=["molecule","basis","nbasis","iters","js_eri_ms","js_scf_ms","js_total_ms","wasm_eri_ms","wasm_scf_ms","wasm_total_ms","speedup_wasm_vs_js","simd_eri_ms","simd_scf_ms","simd_total_ms","speedup_simd_vs_js"],s=a=>isFinite(a)?a.toFixed(2):"",n=(a,r)=>isFinite(a)&&isFinite(r)&&r>0?(a/r).toFixed(2):"",o=t.map(a=>[a.molecule,a.basis,String(a.nbasis),String(a.js.iters),s(a.js.eriMedianMs),s(a.js.scfMedianMs),s(a.js.totalMedianMs),s(a.wasmBase.eriMedianMs),s(a.wasmBase.scfMedianMs),s(a.wasmBase.totalMedianMs),n(a.js.totalMedianMs,a.wasmBase.totalMedianMs),s(a.wasmSimd.eriMedianMs),s(a.wasmSimd.scfMedianMs),s(a.wasmSimd.totalMedianMs),n(a.js.totalMedianMs,a.wasmSimd.totalMedianMs)].map(r=>`"${r}"`).join(","));return e.join(",")+`
`+o.join(`
`)+`
`}function N(t,e,s){const n=new Blob([e],{type:s}),o=URL.createObjectURL(n),a=document.createElement("a");a.href=o,a.download=t,a.click(),setTimeout(()=>URL.revokeObjectURL(o),1e3)}function xe(t){if(!t)return"";const e=["workload_id","workload_label","backend","warmup_ms","run1_ms","run2_ms","run3_ms","median_ms","error"],s=t.workloads.map(o=>[o.id,o.label,t.backend,isFinite(o.warmupMs)?o.warmupMs.toFixed(2):"",o.runs[0]?.toFixed(2)??"",o.runs[1]?.toFixed(2)??"",o.runs[2]?.toFixed(2)??"",isFinite(o.medianMs)?o.medianMs.toFixed(2):"",o.error??""].map(a=>`"${String(a).replace(/"/g,'""')}"`).join(","));return[`# device.ua,"${t.device.ua.replace(/"/g,'""')}"`,`# device.platform,"${t.device.platform}"`,`# device.deviceMemoryGB,"${t.device.deviceMemoryGB??""}"`,`# device.hardwareConcurrency,"${t.device.hardwareConcurrency}"`,`# device.screen,"${t.device.screen}"`,`# backend,"${t.backend}"`,`# initialLoadMs,"${t.initialLoadMs??""}"`,`# timestamp,"${t.device.timestamp}"`].join(`
`)+`
`+e.join(",")+`
`+s.join(`
`)+`
`}function Se(){if(document.getElementById("cb-styles"))return;const t=document.createElement("style");t.id="cb-styles",t.textContent=`
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
  `,document.head.appendChild(t)}ge();
