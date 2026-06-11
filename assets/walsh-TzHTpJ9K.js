import"./styles-B8WF9G-F.js";import{p as $t}from"./parseXYZ-ByfvMHmk.js";import{g as ft,f as mt}from"./ri-DT3jfd0v.js";import{a as yt,H as ut,b as bt}from"./theoryControls-CuXHSm5X.js";import{i as wt,c as vt,r as kt,t as m,a as Mt,b as Et,g as dt,d as j}from"./nav-BKPoqOVV.js";function nt(n,o,c){const l=c*Math.PI/180/2,y=o*Math.sin(l),$=o*Math.cos(l);return`3
${n}H2
${n}  0.0  0.0  0.0
H  ${y.toFixed(6)}  0.0  ${$.toFixed(6)}
H  ${(-y).toFixed(6)}  0.0  ${$.toFixed(6)}`}const xt=[{id:"h2o",labelKey:"walsh.scenH2O",descKey:"walsh.descH2O",charge:0,mult:1,angleMin:80,angleMax:180,steps:21,nCore:1,buildXYZ:n=>nt("O",.96,n)},{id:"beh2",labelKey:"walsh.scenBeH2",descKey:"walsh.descBeH2",charge:0,mult:1,angleMin:80,angleMax:180,steps:21,nCore:1,buildXYZ:n=>nt("Be",1.334,n)},{id:"ch2",labelKey:"walsh.scenCH2",descKey:"walsh.descCH2",charge:0,mult:1,angleMin:80,angleMax:180,steps:21,nCore:1,buildXYZ:n=>nt("C",1.11,n)}];let Y=xt[0],r=[],et=0,G=0,T=!1,st=!1,C=!1,pt=0,it="HF",z=0;const lt=new Map;async function zt(n){const o=lt.get(n);if(o)return o;const c=`/GANSU-Lite/basis/${n.toLowerCase()}.gbs`,x=await fetch(c);if(!x.ok)throw new Error(`Failed to load basis set: ${n}`);const l=await x.text(),y=mt.fromGBS(l);return lt.set(n,y),y}const E=document.getElementById("app");function K(){E.innerHTML=`
    <div class="opt-page">
      ${kt("walsh")}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${m("walsh.molecule")}</h2>
          <div class="opt-scenario-grid" id="scen-grid"></div>

          <div class="theory-row" style="margin:12px 0 8px;display:flex;align-items:center;gap:8px;font-size:0.9rem;">
            <span>Theory:</span>${yt("theory-sel",it,"",ut)}
          </div>

          <button id="run-btn" class="opt-run-btn" ${T?"disabled":""}>
            ${T?m("walsh.running"):m("walsh.run")}
          </button>
          ${T?`<button id="stop-btn" class="opt-stop-btn">${m("walsh.stop")}</button>`:""}

          <div id="progress-area"></div>
          <div id="summary-area"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="graph-container">
            ${!C&&!T?`<p class="opt-hint">${m("walsh.waiting")}</p>`:""}
          </div>
          <div id="slider-area"></div>
          <div id="mol-vis"></div>
        </div>
      </div>
    </div>`,qt();const n=E.querySelector("#scen-grid"),o=document.createElement("div");o.className="opt-category-row";for(const x of xt){const l=document.createElement("div");l.className="opt-scenario-card"+(x.id===Y.id?" selected":""),l.innerHTML=`<strong>${m(x.labelKey)}</strong><span class="conv-desc">${m(x.descKey)}</span>`,l.addEventListener("click",()=>{T||(Y=x,r=[],C=!1,K())}),o.appendChild(l)}n.appendChild(o),E.querySelector("#nav-theme").addEventListener("click",()=>{Mt(),K()}),E.querySelector("#nav-lang").addEventListener("click",()=>{Et(),K()});const c=E.querySelector("#theory-sel");c&&c.addEventListener("change",()=>{it=c.value}),E.querySelector("#run-btn").addEventListener("click",()=>{T||Ft()}),E.querySelector("#stop-btn")?.addEventListener("click",()=>{st=!0}),r.length>0&&(at(),Ct(),C&&(Tt(),ht()))}async function Ft(){T=!0,st=!1,r=[],C=!1,et=0,G=0,K();const n=Y,o=performance.now(),c=await zt("STO-3G"),x=(n.angleMax-n.angleMin)/(n.steps-1);for(let l=0;l<n.steps&&!st;l++){const y=n.angleMin+l*x;St(y,l,n.steps);try{const $=n.buildXYZ(y),O=$t($),i=Math.floor((n.mult-1)/2),h=new ft(O,c,n.charge,i);l===0&&(G=h.numBasis,et=h.numAlphaSpins);const g=await bt(h,c,it),L=await g.solve({eriBackend:"js"}),b=Array.from(g.orbitalEnergies);r.push({angle:y,orbitalEnergies:b,totalEnergy:L}),at()}catch($){console.error(`Error at angle ${y}°:`,$)}await new Promise($=>setTimeout($,0))}pt=performance.now()-o,T=!1,C=!0,K()}function St(n,o,c){const x=E.querySelector("#progress-area");if(!x)return;const l=(o/c*100).toFixed(0);x.innerHTML=`
    <div class="opt-progress"><div class="opt-progress-bar" style="width:${l}%"></div></div>
    <div class="opt-progress-text">${n.toFixed(0)}° (${o+1}/${c})</div>`}function Ct(){const n=E.querySelector("#summary-area");if(!n||r.length===0)return;let o=0;for(let $=1;$<r.length;$++)r[$].totalEnergy<r[o].totalEnergy&&(o=$);const c=G-Y.nCore,x=et-Y.nCore,l=(pt/1e3).toFixed(1);let y=`<div class="opt-summary">
    <h3>${C?m("walsh.done"):m("walsh.running")}</h3>
    <table>
      <tr><td style="color:var(--color-text-dim)">${m("walsh.basis")}</td><td>STO-3G (M=${G})</td></tr>
      <tr><td style="color:var(--color-text-dim)">${m("walsh.nOcc")}</td><td>${x} / ${c}</td></tr>
      <tr><td style="color:var(--color-text-dim)">${m("walsh.minAngle")}</td><td>${r[o].angle.toFixed(1)}°</td></tr>
      <tr><td style="color:var(--color-text-dim)">${m("walsh.minEnergy")}</td><td style="font-family:monospace">${r[o].totalEnergy.toFixed(6)} Eh</td></tr>
      ${C?`<tr><td style="color:var(--color-text-dim)">${m("walsh.time")}</td><td>${l}s</td></tr>`:""}
    </table>
  </div>`;n.innerHTML=y}const Ht=["#d62728","#1f77b4","#2ca02c","#9467bd","#ff7f0e","#8c564b","#e377c2","#17becf","#bcbd22","#7f7f7f"],Lt=["#ff6b6b","#4dabf7","#51cf66","#cc5de8","#ffa94d","#d4a373","#f783ac","#38d9f5","#d9e363","#adb5bd"];function at(){const n=E.querySelector("#graph-container");if(!n||r.length===0)return;const o=dt(),c=j()?Lt:Ht,x=Y,l=x.nCore,y=G,$=560,O=380,i=72,h=52,g=36,L=44,b=$-i-h,w=O-g-L,X=[];for(const t of r)for(let s=l;s<t.orbitalEnergies.length;s++)X.push(t.orbitalEnergies[s]);if(X.length===0)return;const D=Math.min(...X),P=Math.max(...X),N=Math.max((P-D)*.08,.05),R=D-N,H=P+N,M=t=>i+(t-x.angleMin)/(x.angleMax-x.angleMin)*b,F=t=>g+w-(t-R)/(H-R)*w;let e=`<svg width="${$}" height="${O}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${$} ${O}" style="max-width:100%;">`;e+=`<rect x="${i}" y="${g}" width="${b}" height="${w}" fill="${o.surface}" rx="2"/>`;for(let t=0;t<=5;t++){const s=R+(H-R)*t/5,a=F(s);e+=`<line x1="${i}" y1="${a}" x2="${i+b}" y2="${a}" stroke="${o.grid}" stroke-width="0.5"/>`,e+=`<text x="${i-6}" y="${a+3}" text-anchor="end" font-size="9" fill="${o.dim}">${s.toFixed(2)}</text>`}for(let t=x.angleMin;t<=x.angleMax+.1;t+=20){const s=Math.round(t),a=M(s);e+=`<line x1="${a}" y1="${g}" x2="${a}" y2="${g+w}" stroke="${o.grid}" stroke-width="0.3"/>`,e+=`<line x1="${a}" y1="${g+w}" x2="${a}" y2="${g+w+4}" stroke="${o.axis}" stroke-width="1"/>`,e+=`<text x="${a}" y="${g+w+16}" text-anchor="middle" font-size="9" fill="${o.dim}">${s}°</text>`}e+=`<line x1="${i}" y1="${g}" x2="${i}" y2="${g+w}" stroke="${o.axis}" stroke-width="1"/>`,e+=`<line x1="${i}" y1="${g+w}" x2="${i+b}" y2="${g+w}" stroke="${o.axis}" stroke-width="1"/>`;for(let t=l;t<y;t++){const s=t-l,a=t<et,u=c[s%c.length],ot=a?2:1.2,rt=a?"":' stroke-dasharray="6,3"';let _="";for(let S=0;S<r.length;S++){if(t>=r[S].orbitalEnergies.length)continue;const Q=M(r[S].angle),tt=F(r[S].orbitalEnergies[t]);_+=S===0?`M${Q.toFixed(1)},${tt.toFixed(1)}`:` L${Q.toFixed(1)},${tt.toFixed(1)}`}if(_&&(e+=`<path d="${_}" fill="none" stroke="${u}" stroke-width="${ot}"${rt}/>`),r.length>0){const S=r[r.length-1];if(t<S.orbitalEnergies.length){const Q=M(S.angle)+4,tt=F(S.orbitalEnergies[t])+3,gt=a?`${s+1}`:`${s+1}*`;e+=`<text x="${Q}" y="${tt}" font-size="8" fill="${u}" font-weight="${a?"600":"400"}">${gt}</text>`}}}const d=g+10;if(e+=`<line x1="${i+8}" y1="${d}" x2="${i+28}" y2="${d}" stroke="${o.dim}" stroke-width="2"/>`,e+=`<text x="${i+32}" y="${d+3}" font-size="7.5" fill="${o.dim}">${m("walsh.occupied")}</text>`,e+=`<line x1="${i+8}" y1="${d+14}" x2="${i+28}" y2="${d+14}" stroke="${o.dim}" stroke-width="1.2" stroke-dasharray="6,3"/>`,e+=`<text x="${i+32}" y="${d+17}" font-size="7.5" fill="${o.dim}">${m("walsh.virtual")}</text>`,C&&z>=0&&z<r.length){const t=r[z],s=M(t.angle),a=j()?"#ffd700":"#cc8800";e+=`<line x1="${s}" y1="${g}" x2="${s}" y2="${g+w}" stroke="${a}" stroke-width="1" stroke-dasharray="3,2" opacity="0.6"/>`;for(let u=l;u<y;u++){if(u>=t.orbitalEnergies.length)continue;const ot=u-l,rt=c[ot%c.length],_=F(t.orbitalEnergies[u]);e+=`<circle cx="${s}" cy="${_}" r="4.5" fill="${rt}" stroke="#fff" stroke-width="1"/>`}}e+=`<text x="${i+b/2}" y="20" text-anchor="middle" font-size="12" font-weight="600" fill="${o.titleSvg}">${m("walsh.graphTitle")}</text>`,e+=`<text x="${i+b/2}" y="${O-4}" text-anchor="middle" font-size="10" fill="${o.dim}">${m("walsh.xAngle")}</text>`,e+=`<text x="14" y="${g+w/2}" text-anchor="middle" font-size="10" fill="${o.dim}" transform="rotate(-90,14,${g+w/2})">${m("walsh.yOrbEnergy")}</text>`,e+="</svg>";const f=140,p=16,k=f-p-36,B=r.map(t=>t.totalEnergy),U=Math.min(...B),W=Math.max(...B),V=Math.max((W-U)*.15,.002),I=U-V,J=W+V,A=t=>p+k-(t-I)/(J-I)*k;let v=`<svg width="${$}" height="${f}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${$} ${f}" style="max-width:100%;margin-top:8px;">`;v+=`<rect x="${i}" y="${p}" width="${b}" height="${k}" fill="${o.surface}" rx="2"/>`;for(let t=0;t<=3;t++){const s=I+(J-I)*t/3,a=A(s);v+=`<line x1="${i}" y1="${a}" x2="${i+b}" y2="${a}" stroke="${o.grid}" stroke-width="0.5"/>`,v+=`<text x="${i-6}" y="${a+3}" text-anchor="end" font-size="8" fill="${o.dim}">${s.toFixed(4)}</text>`}v+=`<line x1="${i}" y1="${p}" x2="${i}" y2="${p+k}" stroke="${o.axis}" stroke-width="1"/>`,v+=`<line x1="${i}" y1="${p+k}" x2="${i+b}" y2="${p+k}" stroke="${o.axis}" stroke-width="1"/>`;for(let t=x.angleMin;t<=x.angleMax+.1;t+=20){const s=Math.round(t),a=M(s);v+=`<line x1="${a}" y1="${p+k}" x2="${a}" y2="${p+k+4}" stroke="${o.axis}" stroke-width="1"/>`,v+=`<text x="${a}" y="${p+k+14}" text-anchor="middle" font-size="8" fill="${o.dim}">${s}°</text>`}const Z=j()?"#00d4ff":"#0077cc";if(r.length>=2){let t="";for(let s=0;s<r.length;s++){const a=M(r[s].angle),u=A(r[s].totalEnergy);t+=s===0?`M${a.toFixed(1)},${u.toFixed(1)}`:` L${a.toFixed(1)},${u.toFixed(1)}`}v+=`<path d="${t}" fill="none" stroke="${Z}" stroke-width="2"/>`}for(const t of r)v+=`<circle cx="${M(t.angle).toFixed(1)}" cy="${A(t.totalEnergy).toFixed(1)}" r="2.5" fill="${Z}"/>`;if(r.length>0){let t=0;for(let u=1;u<r.length;u++)r[u].totalEnergy<r[t].totalEnergy&&(t=u);const s=M(r[t].angle),a=A(r[t].totalEnergy);v+=`<circle cx="${s}" cy="${a}" r="5" fill="none" stroke="${Z}" stroke-width="1.5"/>`,v+=`<text x="${s}" y="${a-8}" text-anchor="middle" font-size="8" fill="${Z}">${r[t].angle.toFixed(0)}°</text>`}if(C&&z>=0&&z<r.length){const t=r[z],s=M(t.angle),a=A(t.totalEnergy),u=j()?"#ffd700":"#cc8800";v+=`<line x1="${s}" y1="${p}" x2="${s}" y2="${p+k}" stroke="${u}" stroke-width="1" stroke-dasharray="3,2" opacity="0.6"/>`,v+=`<circle cx="${s}" cy="${a}" r="5" fill="${Z}" stroke="${u}" stroke-width="2"/>`}v+=`<text x="${i+b/2}" y="10" text-anchor="middle" font-size="10" font-weight="600" fill="${o.titleSvg}">${m("walsh.totalEnergy")}</text>`,v+=`<text x="${i+b/2}" y="${f-2}" text-anchor="middle" font-size="9" fill="${o.dim}">${m("walsh.xAngle")}</text>`,v+=`<text x="14" y="${p+k/2}" text-anchor="middle" font-size="9" fill="${o.dim}" transform="rotate(-90,14,${p+k/2})">${m("walsh.yTotalEnergy")}</text>`,v+="</svg>",n.innerHTML=e+v}function Tt(){const n=E.querySelector("#slider-area");!n||r.length===0||(z>=r.length&&(z=0),n.innerHTML=`
    <div class="walsh-slider">
      <input type="range" id="angle-slider" min="0" max="${r.length-1}" value="${z}" />
      <div class="walsh-slider-label">${r[z].angle.toFixed(0)}°</div>
    </div>`,E.querySelector("#angle-slider")?.addEventListener("input",o=>{z=parseInt(o.target.value);const c=E.querySelector(".walsh-slider-label");c&&(c.textContent=`${r[z].angle.toFixed(0)}°`),at(),ht()}))}const Ot={H:"#999",Be:"#6C0",C:"#555",O:"#F22"};function ht(){const n=E.querySelector("#mol-vis");if(!n||!C||r.length===0)return;const o=dt(),c=Y,x=r[z].angle,l=220,y=180,$=30,i=c.buildXYZ(x).split(`
`),h=[];for(let e=2;e<i.length;e++){const d=i[e].trim().split(/\s+/);d.length>=4&&h.push({sym:d[0],x:parseFloat(d[1]),z:parseFloat(d[3])})}let g=0,L=0,b=0,w=0;for(const e of[c.angleMin,c.angleMax]){const f=c.buildXYZ(e).split(`
`);for(let p=2;p<f.length;p++){const q=f[p].trim().split(/\s+/),k=parseFloat(q[1]),B=parseFloat(q[3]);g=Math.min(g,k),L=Math.max(L,k),b=Math.min(b,B),w=Math.max(w,B)}}const X=L-g||.01,D=w-b||.01,P=Math.min((l-2*$)/X,(y-2*$)/D),N=(g+L)/2,R=(b+w)/2,H=e=>l/2+(e-N)*P,M=e=>y/2-(e-R)*P;let F=`<svg width="${l}" height="${y}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${l} ${y}" style="max-width:100%;">`;for(let e=0;e<h.length;e++)for(let d=e+1;d<h.length;d++){const f=h[d].x-h[e].x,p=h[d].z-h[e].z;Math.sqrt(f*f+p*p)<2.5&&(F+=`<line x1="${H(h[e].x)}" y1="${M(h[e].z)}" x2="${H(h[d].x)}" y2="${M(h[d].z)}" stroke="${o.grid}" stroke-width="4" stroke-linecap="round"/>`)}for(const e of h){const d=H(e.x),f=M(e.z),p=e.sym==="H"?14:18;F+=`<circle cx="${d}" cy="${f}" r="${p}" fill="${Ot[e.sym]||"#888"}" stroke="${o.axis}" stroke-width="1.5"/>`,F+=`<text x="${d}" y="${f+4}" text-anchor="middle" font-size="11" font-weight="bold" fill="#fff">${e.sym}</text>`}if(h.length>=3){const e=H(h[0].x),d=M(h[0].z),f=28,p=Math.atan2(M(h[1].z)-d,H(h[1].x)-e),q=Math.atan2(M(h[2].z)-d,H(h[2].x)-e),k=e+f*Math.cos(p),B=d+f*Math.sin(p),U=e+f*Math.cos(q),W=d+f*Math.sin(q),V=j()?"#ffd700":"#cc8800";F+=`<path d="M ${k.toFixed(1)} ${B.toFixed(1)} A ${f} ${f} 0 0 0 ${U.toFixed(1)} ${W.toFixed(1)}" fill="none" stroke="${V}" stroke-width="1.5"/>`;const I=(p+q)/2,J=e+(f+14)*Math.cos(I),A=d+(f+14)*Math.sin(I);F+=`<text x="${J.toFixed(1)}" y="${A.toFixed(1)}" text-anchor="middle" font-size="11" font-weight="600" fill="${V}">${x.toFixed(0)}°</text>`}F+="</svg>",n.innerHTML=F}let ct=!1;function qt(){if(ct)return;ct=!0;const n=document.createElement("style");n.textContent=`
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
    #mol-vis { width: 100%; text-align: center; }
    .opt-hint { color: var(--color-text-dim); font-size: 0.85rem; padding: 60px 20px; }

    .walsh-slider {
      width: 100%; max-width: 560px; margin: 8px auto 0;
      display: flex; align-items: center; gap: 10px;
    }
    .walsh-slider input[type=range] { flex: 1; cursor: pointer; }
    .walsh-slider-label {
      font-size: 0.85rem; font-weight: 600; min-width: 40px; text-align: center;
      color: var(--color-text);
    }

    @media (max-width: 700px) {
      .opt-content { flex-direction: column; }
      .opt-controls { flex: none; }
    }
  `,document.head.appendChild(n)}wt();vt();K();
