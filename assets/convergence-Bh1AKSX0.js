import"./styles-B8WF9G-F.js";import{p as V}from"./parseXYZ-BE3rbSP_.js";import{g as D,f as A}from"./ri-Cuj4t-H2.js";import{a as Z,H as j,b as _}from"./theoryControls-AWp8GxPQ.js";import{i as U,c as X,r as Y,t as d,a as J,b as Q,g as W,d as ee}from"./nav-BKPoqOVV.js";const O=[{id:"h2",labelKey:"conv.scenH2",descKey:"conv.descH2",charge:0,mult:1,xyz:`2
H2
H  0.0  0.0  0.0
H  0.0  0.0  0.740000`},{id:"hf",labelKey:"conv.scenHF",descKey:"conv.descHF",charge:0,mult:1,xyz:`2
HF
H  0.0  0.0  0.0
F  0.0  0.0  0.917000`},{id:"h2o",labelKey:"conv.scenH2O",descKey:"conv.descH2O",charge:0,mult:1,xyz:(()=>{const t=104*Math.PI/180,o=t/2,r=.96*Math.sin(o),e=.96*Math.cos(o);return`3
H2O
O  0.0  0.0  0.0
H  ${r.toFixed(6)}  0.0  ${e.toFixed(6)}
H  ${(-r).toFixed(6)}  0.0  ${e.toFixed(6)}`})()},{id:"lih",labelKey:"conv.scenLiH",descKey:"conv.descLiH",charge:0,mult:1,xyz:`2
LiH
Li  0.0  0.0  0.0
H   0.0  0.0  1.596000`},{id:"n2",labelKey:"conv.scenN2",descKey:"conv.descN2",charge:0,mult:1,xyz:`2
N2
N  0.0  0.0  0.0
N  0.0  0.0  1.098000`},{id:"nh3",labelKey:"conv.scenNH3",descKey:"conv.descNH3",charge:0,mult:1,xyz:(()=>{const o=Math.sqrt(.878983),r=Math.sqrt(3)/2;return["4","NH3",`N   0.000000  0.000000  ${.381.toFixed(6)}`,`H   ${o.toFixed(6)}  0.000000  0.000000`,`H   ${(-o/2).toFixed(6)}  ${(o*r).toFixed(6)}  0.000000`,`H   ${(-o/2).toFixed(6)}  ${(-o*r).toFixed(6)}  0.000000`].join(`
`)})()}],b=["STO-3G","STO-6G","3-21G","6-31G","6-31G(d,p)","cc-pVDZ","def2-SVP","def2-TZVP"],te=new Set(["cc-pVTZ","def2-TZVP"]);let k=O[0],g=[],u=!1,L=!1,$=!1,R=0,S="HF";const C=new Map;async function oe(a){const t=C.get(a);if(t)return t;const o=`/GANSU-Lite/basis/${a.toLowerCase()}.gbs`,r=await fetch(o);if(!r.ok)throw new Error(`Failed to load basis set: ${a}`);const e=await r.text(),n=A.fromGBS(e);return C.set(a,n),n}const m=document.getElementById("app");function v(){m.innerHTML=`
    <div class="opt-page">
      ${Y("convergence")}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${d("conv.molecule")}</h2>
          <div class="opt-scenario-grid" id="mol-grid"></div>

          <div class="theory-row" style="margin:8px 0;display:flex;align-items:center;gap:8px;font-size:0.9rem;">
            <span>Theory:</span>${Z("theory-sel",S,"",j)}
          </div>

          <button id="run-btn" class="opt-run-btn" ${u?"disabled":""}>
            ${u?d("conv.running"):d("conv.run")}
          </button>
          ${u?`<button id="stop-btn" class="opt-stop-btn">${d("conv.stop")}</button>`:""}

          <div id="progress-area"></div>
          <div id="result-table"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="graph-container">
            ${!$&&!u?`<p class="opt-hint">${d("conv.waiting")}</p>`:""}
          </div>
        </div>
      </div>
    </div>`,se();const a=m.querySelector("#mol-grid"),t=document.createElement("div");t.className="opt-category-row";for(const r of O){const e=document.createElement("div");e.className="opt-scenario-card"+(r.id===k.id?" selected":""),e.innerHTML=`<strong>${d(r.labelKey)}</strong><span class="conv-desc">${d(r.descKey)}</span>`,e.addEventListener("click",()=>{u||(k=r,g=[],$=!1,v())}),t.appendChild(e)}a.appendChild(t),m.querySelector("#nav-theme").addEventListener("click",()=>{J(),v()}),m.querySelector("#nav-lang").addEventListener("click",()=>{Q(),v()});const o=m.querySelector("#theory-sel");o&&o.addEventListener("change",()=>{S=o.value}),m.querySelector("#run-btn").addEventListener("click",()=>{u||re()}),m.querySelector("#stop-btn")?.addEventListener("click",()=>{L=!0}),g.length>0&&(P(),G())}async function re(){u=!0,L=!1,g=[],$=!1,v();const a=Math.floor((k.mult-1)/2),t=performance.now(),o=V(k.xyz),r=S!=="HF";for(let e=0;e<b.length&&!L;e++){const n=b[e];if(r&&te.has(n)){g.push({basis:n,nbasis:0,energy:NaN,timeMs:0,error:"skipped (DFT × TZ too slow)"});continue}ne(n,e,b.length);try{const i=performance.now(),f=await oe(n),h=new D(o,f,k.charge,a),x=await(await _(h,f,S)).solve({eriBackend:"js"}),H=performance.now()-i;g.push({basis:n,nbasis:h.numBasis,energy:x,timeMs:H}),P(),G()}catch(i){const f=i instanceof Error?i.message:String(i);console.error(`Error computing ${n}:`,i),g.push({basis:n,nbasis:0,energy:NaN,timeMs:0,error:f})}await new Promise(i=>setTimeout(i,0))}R=performance.now()-t,u=!1,$=!0,v()}function ne(a,t,o){const r=m.querySelector("#progress-area");if(!r)return;const e=(t/o*100).toFixed(0);r.innerHTML=`
    <div class="opt-progress"><div class="opt-progress-bar" style="width:${e}%"></div></div>
    <div class="opt-progress-text">${a}... (${t}/${o})</div>`}function G(){const a=m.querySelector("#result-table");if(!a||g.length===0)return;const t=g.filter(n=>isFinite(n.energy)),o=t.length>0?t[t.length-1]:null,r=$?` (${(R/1e3).toFixed(1)}s)`:"";let e=`<div class="opt-summary">
    <h3>${$?d("conv.done")+r:d("conv.running")}</h3>
    <table>
      <tr>
        <th>${d("conv.colBasis")}</th>
        <th style="text-align:right">M</th>
        <th style="text-align:right">${d("conv.colEnergy")}</th>
        <th style="text-align:right">${d("conv.colTime")}</th>
      </tr>`;for(const n of g){const i=o&&n.basis===o.basis,f=n.timeMs<1e3?`${n.timeMs.toFixed(0)} ms`:`${(n.timeMs/1e3).toFixed(1)} s`;n.error?e+=`<tr style="color:var(--color-error)">
        <td>${n.basis}</td>
        <td></td>
        <td colspan="2" style="font-size:0.7rem">${n.error}</td>
      </tr>`:e+=`<tr${i?' style="font-weight:600"':""}>
        <td>${n.basis}</td>
        <td style="text-align:right">${n.nbasis}</td>
        <td style="text-align:right;font-family:monospace">${n.energy.toFixed(6)}</td>
        <td style="text-align:right">${f}</td>
      </tr>`}e+="</table></div>",a.innerHTML=e}function P(){const a=m.querySelector("#graph-container");if(!a||g.length===0)return;const t=W(),o=520,r=360,e=80,n=24,i=36,f=44,h=o-e-n,p=r-i-f,x=g.filter(s=>isFinite(s.energy));if(x.length===0)return;const H=x.map(s=>s.energy),E=Math.min(...H),N=Math.max(...H),K=Math.max((N-E)*.15,.005),F=E-K,B=N+K,z=s=>e+s/(b.length-1)*h,M=s=>i+p-(s-F)/(B-F)*p;let l=`<svg width="${o}" height="${r}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${o} ${r}" style="max-width:100%;">`;l+=`<rect x="${e}" y="${i}" width="${h}" height="${p}" fill="${t.surface}" rx="2"/>`;for(let s=0;s<=5;s++){const c=F+(B-F)*s/5,y=M(c);l+=`<line x1="${e}" y1="${y}" x2="${e+h}" y2="${y}" stroke="${t.grid}" stroke-width="0.5"/>`,l+=`<text x="${e-6}" y="${y+3}" text-anchor="end" font-size="9" fill="${t.dim}">${c.toFixed(4)}</text>`}l+=`<line x1="${e}" y1="${i}" x2="${e}" y2="${i+p}" stroke="${t.axis}" stroke-width="1"/>`,l+=`<line x1="${e}" y1="${i+p}" x2="${e+h}" y2="${i+p}" stroke="${t.axis}" stroke-width="1"/>`;for(let s=0;s<b.length;s++){const c=z(s);l+=`<line x1="${c}" y1="${i+p}" x2="${c}" y2="${i+p+4}" stroke="${t.axis}" stroke-width="1"/>`,l+=`<text x="${c}" y="${i+p+16}" text-anchor="middle" font-size="8.5" fill="${t.dim}">${b[s]}</text>`}const T=ee()?"#00d4ff":"#0077cc";if(x.length>=2){let s="";for(let c=0;c<x.length;c++){const y=b.indexOf(x[c].basis),w=z(y),q=M(x[c].energy);s+=c===0?`M${w.toFixed(1)},${q.toFixed(1)}`:` L${w.toFixed(1)},${q.toFixed(1)}`}l+=`<path d="${s}" fill="none" stroke="${T}" stroke-width="2"/>`}for(const s of x){const c=b.indexOf(s.basis),y=z(c),w=M(s.energy);l+=`<circle cx="${y}" cy="${w}" r="4" fill="${T}"/>`,l+=`<text x="${y}" y="${w-8}" text-anchor="middle" font-size="8" fill="${t.dim}">${s.nbasis}</text>`}if(x.length>=2){const s=x[x.length-1].energy,c=M(s);l+=`<line x1="${e}" y1="${c}" x2="${e+h}" y2="${c}" stroke="${T}" stroke-width="0.8" stroke-dasharray="4,3" opacity="0.5"/>`}l+=`<text x="${e+h/2}" y="20" text-anchor="middle" font-size="12" font-weight="600" fill="${t.titleSvg}">${d("conv.graphTitle")}</text>`,l+=`<text x="${e+h/2}" y="${r-4}" text-anchor="middle" font-size="10" fill="${t.dim}">${d("conv.xBasis")}</text>`,l+=`<text x="14" y="${i+p/2}" text-anchor="middle" font-size="10" fill="${t.dim}" transform="rotate(-90,14,${i+p/2})">${d("conv.yEnergy")}</text>`,l+="</svg>",a.innerHTML=l}let I=!1;function se(){if(I)return;I=!0;const a=document.createElement("style");a.textContent=`
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      background: var(--color-bg);
      color: var(--color-text);
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
      min-height: 100vh;
    }
    .opt-page { max-width: 960px; margin: 0 auto; padding: 16px 20px; }

    .opt-header {
      display: flex; justify-content: space-between; align-items: center;
      margin-bottom: 20px; padding-bottom: 12px;
      border-bottom: 1px solid var(--color-border);
    }
    .opt-header-left { display: flex; flex-direction: column; gap: 2px; }
    .opt-back-link { font-size: 0.75rem; color: var(--color-link); text-decoration: none; }
    .opt-back-link:hover { text-decoration: underline; }
    .opt-title { font-size: 1rem; font-weight: 600; }
    .opt-header-right { display: flex; gap: 6px; }
    .opt-header-right button {
      background: none; border: 1px solid var(--color-border); border-radius: 6px;
      padding: 4px 8px; cursor: pointer; color: var(--color-text); font-size: 0.8rem;
    }
    .opt-header-right button:hover { background: var(--color-surface-alt); }

    .opt-content { display: flex; gap: 20px; }
    .opt-panel {
      background: var(--color-surface); border: 1px solid var(--color-border);
      border-radius: 10px; padding: 16px 18px;
    }
    .opt-controls { flex: 0 0 300px; }
    .opt-graph-panel { flex: 1; min-width: 0; display: flex; flex-direction: column; align-items: center; padding-top: 8px; }
    .opt-controls h2 {
      font-size: 0.78rem; font-weight: 600; color: var(--color-text-secondary);
      text-transform: uppercase; letter-spacing: 0.04em;
      margin: 14px 0 8px;
    }
    .opt-controls h2:first-child { margin-top: 0; }

    .opt-scenario-grid { display: flex; flex-direction: column; gap: 4px; }
    .opt-category-row { display: flex; flex-wrap: wrap; gap: 5px; }
    .opt-scenario-card {
      padding: 5px 10px; border: 1px solid var(--color-border); border-radius: 6px;
      cursor: pointer; transition: all 0.15s; flex: 0 0 auto;
    }
    .opt-scenario-card:hover { background: var(--color-surface-alt); }
    .opt-scenario-card.selected {
      border-color: var(--color-accent); background: var(--color-surface-alt);
      box-shadow: 0 0 0 1px var(--color-accent);
    }
    .opt-scenario-card strong { display: block; font-size: 0.78rem; white-space: nowrap; }
    .conv-desc { display: block !important; font-size: 0.68rem; color: var(--color-text-dim); }

    .opt-run-btn {
      width: 100%; margin-top: 14px; padding: 10px; border: none; border-radius: 8px;
      font-size: 0.85rem; font-weight: 600; cursor: pointer;
      background: var(--color-accent); color: var(--color-accent-on);
      transition: background 0.15s;
    }
    .opt-run-btn:hover:not([disabled]) { background: var(--color-accent-hover); }
    .opt-run-btn[disabled] { opacity: 0.6; cursor: not-allowed; }
    .opt-stop-btn {
      width: 100%; margin-top: 6px; padding: 8px; border: 1px solid var(--color-error); border-radius: 8px;
      font-size: 0.8rem; font-weight: 600; cursor: pointer;
      background: none; color: var(--color-error);
    }

    .opt-progress {
      height: 6px; background: var(--color-progress-bg, #e0e4ea);
      border-radius: 3px; margin-top: 12px; overflow: hidden;
    }
    .opt-progress-bar {
      height: 100%; background: var(--color-accent); border-radius: 3px;
      transition: width 0.2s ease;
    }
    .opt-progress-text {
      font-size: 0.72rem; color: var(--color-text-dim); text-align: center; margin-top: 4px;
    }

    .opt-summary {
      margin-top: 14px; padding: 12px; background: var(--color-surface-alt);
      border-radius: 8px; border: 1px solid var(--color-border);
    }
    .opt-summary h3 {
      font-size: 0.82rem; color: var(--color-converged); margin-bottom: 8px;
    }
    .opt-summary table { width: 100%; font-size: 0.78rem; border-collapse: collapse; }
    .opt-summary th {
      font-size: 0.68rem; color: var(--color-text-dim); text-align: left;
      padding: 2px 4px; font-weight: 500;
    }
    .opt-summary td { padding: 3px 4px; }

    #graph-container { width: 100%; text-align: center; }
    .opt-hint { color: var(--color-text-dim); font-size: 0.85rem; padding: 60px 20px; }

    @media (max-width: 700px) {
      .opt-content { flex-direction: column; }
      .opt-controls { flex: none; }
    }
  `,document.head.appendChild(a)}U();X();v();
