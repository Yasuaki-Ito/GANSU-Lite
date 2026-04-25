import{i as j,c as P,r as A,t as d,a as U,b as D,g as V,d as X}from"./nav-BaYBNpRD.js";import{p as Y}from"./parseXYZ-DTmoiCYR.js";import{t as Z,u as _,d as J}from"./builder-B-MUn5dT.js";const O=[{id:"h2",labelKey:"conv.scenH2",descKey:"conv.descH2",charge:0,mult:1,xyz:`2
H2
H  0.0  0.0  0.0
H  0.0  0.0  0.740000`},{id:"hf",labelKey:"conv.scenHF",descKey:"conv.descHF",charge:0,mult:1,xyz:`2
HF
H  0.0  0.0  0.0
F  0.0  0.0  0.917000`},{id:"h2o",labelKey:"conv.scenH2O",descKey:"conv.descH2O",charge:0,mult:1,xyz:(()=>{const e=104*Math.PI/180,r=e/2,o=.96*Math.sin(r),t=.96*Math.cos(r);return`3
H2O
O  0.0  0.0  0.0
H  ${o.toFixed(6)}  0.0  ${t.toFixed(6)}
H  ${(-o).toFixed(6)}  0.0  ${t.toFixed(6)}`})()},{id:"lih",labelKey:"conv.scenLiH",descKey:"conv.descLiH",charge:0,mult:1,xyz:`2
LiH
Li  0.0  0.0  0.0
H   0.0  0.0  1.596000`},{id:"n2",labelKey:"conv.scenN2",descKey:"conv.descN2",charge:0,mult:1,xyz:`2
N2
N  0.0  0.0  0.0
N  0.0  0.0  1.098000`},{id:"nh3",labelKey:"conv.scenNH3",descKey:"conv.descNH3",charge:0,mult:1,xyz:(()=>{const r=Math.sqrt(.878983),o=Math.sqrt(3)/2;return["4","NH3",`N   0.000000  0.000000  ${.381.toFixed(6)}`,`H   ${r.toFixed(6)}  0.000000  0.000000`,`H   ${(-r/2).toFixed(6)}  ${(r*o).toFixed(6)}  0.000000`,`H   ${(-r/2).toFixed(6)}  ${(-r*o).toFixed(6)}  0.000000`].join(`
`)})()}],y=["STO-3G","STO-6G","3-21G","6-31G","6-31G(d,p)","cc-pVDZ"];let k=O[0],g=[],u=!1,L=!1,$=!1,R=0;const q=new Map;async function Q(i){const e=q.get(i);if(e)return e;const r=`/GANSU-Lite/basis/${i.toLowerCase()}.gbs`,o=await fetch(r);if(!o.ok)throw new Error(`Failed to load basis set: ${i}`);const t=await o.text(),n=J.fromGBS(t);return q.set(i,n),n}const h=document.getElementById("app");function v(){h.innerHTML=`
    <div class="opt-page">
      ${A("convergence")}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${d("conv.molecule")}</h2>
          <div class="opt-scenario-grid" id="mol-grid"></div>

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
    </div>`,et();const i=h.querySelector("#mol-grid"),e=document.createElement("div");e.className="opt-category-row";for(const r of O){const o=document.createElement("div");o.className="opt-scenario-card"+(r.id===k.id?" selected":""),o.innerHTML=`<strong>${d(r.labelKey)}</strong><span class="conv-desc">${d(r.descKey)}</span>`,o.addEventListener("click",()=>{u||(k=r,g=[],$=!1,v())}),e.appendChild(o)}i.appendChild(e),h.querySelector("#nav-theme").addEventListener("click",()=>{U(),v()}),h.querySelector("#nav-lang").addEventListener("click",()=>{D(),v()}),h.querySelector("#run-btn").addEventListener("click",()=>{u||W()}),h.querySelector("#stop-btn")?.addEventListener("click",()=>{L=!0}),g.length>0&&(I(),G())}async function W(){u=!0,L=!1,g=[],$=!1,v();const i=Math.floor((k.mult-1)/2),e=performance.now(),r=Y(k.xyz);for(let o=0;o<y.length&&!L;o++){const t=y[o];tt(t,o,y.length);try{const n=performance.now(),a=await Q(t),b=new Z(r,a,k.charge,i),p=await _(b,"RHF").solve({eriBackend:"js"}),x=performance.now()-n;g.push({basis:t,nbasis:b.numBasis,energy:p,timeMs:x}),I(),G()}catch(n){const a=n instanceof Error?n.message:String(n);console.error(`Error computing ${t}:`,n),g.push({basis:t,nbasis:0,energy:NaN,timeMs:0,error:a})}await new Promise(n=>setTimeout(n,0))}R=performance.now()-e,u=!1,$=!0,v()}function tt(i,e,r){const o=h.querySelector("#progress-area");if(!o)return;const t=(e/r*100).toFixed(0);o.innerHTML=`
    <div class="opt-progress"><div class="opt-progress-bar" style="width:${t}%"></div></div>
    <div class="opt-progress-text">${i}... (${e}/${r})</div>`}function G(){const i=h.querySelector("#result-table");if(!i||g.length===0)return;const e=g.filter(n=>isFinite(n.energy)),r=e.length>0?e[e.length-1]:null,o=$?` (${(R/1e3).toFixed(1)}s)`:"";let t=`<div class="opt-summary">
    <h3>${$?d("conv.done")+o:d("conv.running")}</h3>
    <table>
      <tr>
        <th>${d("conv.colBasis")}</th>
        <th style="text-align:right">M</th>
        <th style="text-align:right">${d("conv.colEnergy")}</th>
        <th style="text-align:right">${d("conv.colTime")}</th>
      </tr>`;for(const n of g){const a=r&&n.basis===r.basis,b=n.timeMs<1e3?`${n.timeMs.toFixed(0)} ms`:`${(n.timeMs/1e3).toFixed(1)} s`;n.error?t+=`<tr style="color:var(--color-error)">
        <td>${n.basis}</td>
        <td></td>
        <td colspan="2" style="font-size:0.7rem">${n.error}</td>
      </tr>`:t+=`<tr${a?' style="font-weight:600"':""}>
        <td>${n.basis}</td>
        <td style="text-align:right">${n.nbasis}</td>
        <td style="text-align:right;font-family:monospace">${n.energy.toFixed(6)}</td>
        <td style="text-align:right">${b}</td>
      </tr>`}t+="</table></div>",i.innerHTML=t}function I(){const i=h.querySelector("#graph-container");if(!i||g.length===0)return;const e=V(),r=520,o=360,t=80,n=24,a=36,b=44,m=r-t-n,p=o-a-b,x=g.filter(s=>isFinite(s.energy));if(x.length===0)return;const S=x.map(s=>s.energy),N=Math.min(...S),T=Math.max(...S),E=Math.max((T-N)*.15,.005),H=N-E,K=T+E,z=s=>t+s/(y.length-1)*m,F=s=>a+p-(s-H)/(K-H)*p;let l=`<svg width="${r}" height="${o}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${r} ${o}" style="max-width:100%;">`;l+=`<rect x="${t}" y="${a}" width="${m}" height="${p}" fill="${e.surface}" rx="2"/>`;for(let s=0;s<=5;s++){const c=H+(K-H)*s/5,f=F(c);l+=`<line x1="${t}" y1="${f}" x2="${t+m}" y2="${f}" stroke="${e.grid}" stroke-width="0.5"/>`,l+=`<text x="${t-6}" y="${f+3}" text-anchor="end" font-size="9" fill="${e.dim}">${c.toFixed(4)}</text>`}l+=`<line x1="${t}" y1="${a}" x2="${t}" y2="${a+p}" stroke="${e.axis}" stroke-width="1"/>`,l+=`<line x1="${t}" y1="${a+p}" x2="${t+m}" y2="${a+p}" stroke="${e.axis}" stroke-width="1"/>`;for(let s=0;s<y.length;s++){const c=z(s);l+=`<line x1="${c}" y1="${a+p}" x2="${c}" y2="${a+p+4}" stroke="${e.axis}" stroke-width="1"/>`,l+=`<text x="${c}" y="${a+p+16}" text-anchor="middle" font-size="8.5" fill="${e.dim}">${y[s]}</text>`}const M=X()?"#00d4ff":"#0077cc";if(x.length>=2){let s="";for(let c=0;c<x.length;c++){const f=y.indexOf(x[c].basis),w=z(f),B=F(x[c].energy);s+=c===0?`M${w.toFixed(1)},${B.toFixed(1)}`:` L${w.toFixed(1)},${B.toFixed(1)}`}l+=`<path d="${s}" fill="none" stroke="${M}" stroke-width="2"/>`}for(const s of x){const c=y.indexOf(s.basis),f=z(c),w=F(s.energy);l+=`<circle cx="${f}" cy="${w}" r="4" fill="${M}"/>`,l+=`<text x="${f}" y="${w-8}" text-anchor="middle" font-size="8" fill="${e.dim}">${s.nbasis}</text>`}if(x.length>=2){const s=x[x.length-1].energy,c=F(s);l+=`<line x1="${t}" y1="${c}" x2="${t+m}" y2="${c}" stroke="${M}" stroke-width="0.8" stroke-dasharray="4,3" opacity="0.5"/>`}l+=`<text x="${t+m/2}" y="20" text-anchor="middle" font-size="12" font-weight="600" fill="${e.titleSvg}">${d("conv.graphTitle")}</text>`,l+=`<text x="${t+m/2}" y="${o-4}" text-anchor="middle" font-size="10" fill="${e.dim}">${d("conv.xBasis")}</text>`,l+=`<text x="14" y="${a+p/2}" text-anchor="middle" font-size="10" fill="${e.dim}" transform="rotate(-90,14,${a+p/2})">${d("conv.yEnergy")}</text>`,l+="</svg>",i.innerHTML=l}let C=!1;function et(){if(C)return;C=!0;const i=document.createElement("style");i.textContent=`
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
  `,document.head.appendChild(i)}j();P();v();
