import{i as gt,c as $t,r as ft,t as m,a as mt,b as yt,g as ct,d as V}from"./nav-Lc1uTV-6.js";import{p as ut}from"./parseXYZ-DTmoiCYR.js";import{t as bt,u as wt,d as vt}from"./builder-B-MUn5dT.js";function rt(r,o,l){const p=l*Math.PI/180/2,y=o*Math.sin(p),$=o*Math.cos(p);return`3
${r}H2
${r}  0.0  0.0  0.0
H  ${y.toFixed(6)}  0.0  ${$.toFixed(6)}
H  ${(-y).toFixed(6)}  0.0  ${$.toFixed(6)}`}const dt=[{id:"h2o",labelKey:"walsh.scenH2O",descKey:"walsh.descH2O",charge:0,mult:1,angleMin:80,angleMax:180,steps:21,nCore:1,buildXYZ:r=>rt("O",.96,r)},{id:"beh2",labelKey:"walsh.scenBeH2",descKey:"walsh.descBeH2",charge:0,mult:1,angleMin:80,angleMax:180,steps:21,nCore:1,buildXYZ:r=>rt("Be",1.334,r)},{id:"ch2",labelKey:"walsh.scenCH2",descKey:"walsh.descCH2",charge:0,mult:1,angleMin:80,angleMax:180,steps:21,nCore:1,buildXYZ:r=>rt("C",1.11,r)}];let Y=dt[0],n=[],et=0,_=0,O=!1,st=!1,C=!1,xt=0,z=0;const at=new Map;async function kt(r){const o=at.get(r);if(o)return o;const l=`/GANSU-Lite/basis/${r.toLowerCase()}.gbs`,d=await fetch(l);if(!d.ok)throw new Error(`Failed to load basis set: ${r}`);const p=await d.text(),y=vt.fromGBS(p);return at.set(r,y),y}const E=document.getElementById("app");function K(){E.innerHTML=`
    <div class="opt-page">
      ${ft("walsh")}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${m("walsh.molecule")}</h2>
          <div class="opt-scenario-grid" id="scen-grid"></div>

          <button id="run-btn" class="opt-run-btn" ${O?"disabled":""}>
            ${O?m("walsh.running"):m("walsh.run")}
          </button>
          ${O?`<button id="stop-btn" class="opt-stop-btn">${m("walsh.stop")}</button>`:""}

          <div id="progress-area"></div>
          <div id="summary-area"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="graph-container">
            ${!C&&!O?`<p class="opt-hint">${m("walsh.waiting")}</p>`:""}
          </div>
          <div id="slider-area"></div>
          <div id="mol-vis"></div>
        </div>
      </div>
    </div>`,Lt();const r=E.querySelector("#scen-grid"),o=document.createElement("div");o.className="opt-category-row";for(const l of dt){const d=document.createElement("div");d.className="opt-scenario-card"+(l.id===Y.id?" selected":""),d.innerHTML=`<strong>${m(l.labelKey)}</strong><span class="conv-desc">${m(l.descKey)}</span>`,d.addEventListener("click",()=>{O||(Y=l,n=[],C=!1,K())}),o.appendChild(d)}r.appendChild(o),E.querySelector("#nav-theme").addEventListener("click",()=>{mt(),K()}),E.querySelector("#nav-lang").addEventListener("click",()=>{yt(),K()}),E.querySelector("#run-btn").addEventListener("click",()=>{O||Mt()}),E.querySelector("#stop-btn")?.addEventListener("click",()=>{st=!0}),n.length>0&&(it(),Et(),C&&(Ct(),pt()))}async function Mt(){O=!0,st=!1,n=[],C=!1,et=0,_=0,K();const r=Y,o=performance.now(),l=await kt("STO-3G"),d=(r.angleMax-r.angleMin)/(r.steps-1);for(let p=0;p<r.steps&&!st;p++){const y=r.angleMin+p*d;zt(y,p,r.steps);try{const $=r.buildXYZ(y),T=ut($),i=Math.floor((r.mult-1)/2),h=new bt(T,l,r.charge,i);p===0&&(_=h.numBasis,et=h.numAlphaSpins);const g=wt(h,"RHF"),L=await g.solve({eriBackend:"js"}),b=Array.from(g.orbitalEnergies);n.push({angle:y,orbitalEnergies:b,totalEnergy:L}),it()}catch($){console.error(`Error at angle ${y}°:`,$)}await new Promise($=>setTimeout($,0))}xt=performance.now()-o,O=!1,C=!0,K()}function zt(r,o,l){const d=E.querySelector("#progress-area");if(!d)return;const p=(o/l*100).toFixed(0);d.innerHTML=`
    <div class="opt-progress"><div class="opt-progress-bar" style="width:${p}%"></div></div>
    <div class="opt-progress-text">${r.toFixed(0)}° (${o+1}/${l})</div>`}function Et(){const r=E.querySelector("#summary-area");if(!r||n.length===0)return;let o=0;for(let $=1;$<n.length;$++)n[$].totalEnergy<n[o].totalEnergy&&(o=$);const l=_-Y.nCore,d=et-Y.nCore,p=(xt/1e3).toFixed(1);let y=`<div class="opt-summary">
    <h3>${C?m("walsh.done"):m("walsh.running")}</h3>
    <table>
      <tr><td style="color:var(--color-text-dim)">${m("walsh.basis")}</td><td>STO-3G (M=${_})</td></tr>
      <tr><td style="color:var(--color-text-dim)">${m("walsh.nOcc")}</td><td>${d} / ${l}</td></tr>
      <tr><td style="color:var(--color-text-dim)">${m("walsh.minAngle")}</td><td>${n[o].angle.toFixed(1)}°</td></tr>
      <tr><td style="color:var(--color-text-dim)">${m("walsh.minEnergy")}</td><td style="font-family:monospace">${n[o].totalEnergy.toFixed(6)} Eh</td></tr>
      ${C?`<tr><td style="color:var(--color-text-dim)">${m("walsh.time")}</td><td>${p}s</td></tr>`:""}
    </table>
  </div>`;r.innerHTML=y}const Ft=["#d62728","#1f77b4","#2ca02c","#9467bd","#ff7f0e","#8c564b","#e377c2","#17becf","#bcbd22","#7f7f7f"],St=["#ff6b6b","#4dabf7","#51cf66","#cc5de8","#ffa94d","#d4a373","#f783ac","#38d9f5","#d9e363","#adb5bd"];function it(){const r=E.querySelector("#graph-container");if(!r||n.length===0)return;const o=ct(),l=V()?St:Ft,d=Y,p=d.nCore,y=_,$=560,T=380,i=72,h=52,g=36,L=44,b=$-i-h,w=T-g-L,R=[];for(const t of n)for(let s=p;s<t.orbitalEnergies.length;s++)R.push(t.orbitalEnergies[s]);if(R.length===0)return;const N=Math.min(...R),P=Math.max(...R),W=Math.max((P-N)*.08,.05),X=N-W,H=P+W,M=t=>i+(t-d.angleMin)/(d.angleMax-d.angleMin)*b,F=t=>g+w-(t-X)/(H-X)*w;let e=`<svg width="${$}" height="${T}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${$} ${T}" style="max-width:100%;">`;e+=`<rect x="${i}" y="${g}" width="${b}" height="${w}" fill="${o.surface}" rx="2"/>`;for(let t=0;t<=5;t++){const s=X+(H-X)*t/5,a=F(s);e+=`<line x1="${i}" y1="${a}" x2="${i+b}" y2="${a}" stroke="${o.grid}" stroke-width="0.5"/>`,e+=`<text x="${i-6}" y="${a+3}" text-anchor="end" font-size="9" fill="${o.dim}">${s.toFixed(2)}</text>`}for(let t=d.angleMin;t<=d.angleMax+.1;t+=20){const s=Math.round(t),a=M(s);e+=`<line x1="${a}" y1="${g}" x2="${a}" y2="${g+w}" stroke="${o.grid}" stroke-width="0.3"/>`,e+=`<line x1="${a}" y1="${g+w}" x2="${a}" y2="${g+w+4}" stroke="${o.axis}" stroke-width="1"/>`,e+=`<text x="${a}" y="${g+w+16}" text-anchor="middle" font-size="9" fill="${o.dim}">${s}°</text>`}e+=`<line x1="${i}" y1="${g}" x2="${i}" y2="${g+w}" stroke="${o.axis}" stroke-width="1"/>`,e+=`<line x1="${i}" y1="${g+w}" x2="${i+b}" y2="${g+w}" stroke="${o.axis}" stroke-width="1"/>`;for(let t=p;t<y;t++){const s=t-p,a=t<et,u=l[s%l.length],ot=a?2:1.2,nt=a?"":' stroke-dasharray="6,3"';let G="";for(let S=0;S<n.length;S++){if(t>=n[S].orbitalEnergies.length)continue;const Q=M(n[S].angle),tt=F(n[S].orbitalEnergies[t]);G+=S===0?`M${Q.toFixed(1)},${tt.toFixed(1)}`:` L${Q.toFixed(1)},${tt.toFixed(1)}`}if(G&&(e+=`<path d="${G}" fill="none" stroke="${u}" stroke-width="${ot}"${nt}/>`),n.length>0){const S=n[n.length-1];if(t<S.orbitalEnergies.length){const Q=M(S.angle)+4,tt=F(S.orbitalEnergies[t])+3,ht=a?`${s+1}`:`${s+1}*`;e+=`<text x="${Q}" y="${tt}" font-size="8" fill="${u}" font-weight="${a?"600":"400"}">${ht}</text>`}}}const c=g+10;if(e+=`<line x1="${i+8}" y1="${c}" x2="${i+28}" y2="${c}" stroke="${o.dim}" stroke-width="2"/>`,e+=`<text x="${i+32}" y="${c+3}" font-size="7.5" fill="${o.dim}">${m("walsh.occupied")}</text>`,e+=`<line x1="${i+8}" y1="${c+14}" x2="${i+28}" y2="${c+14}" stroke="${o.dim}" stroke-width="1.2" stroke-dasharray="6,3"/>`,e+=`<text x="${i+32}" y="${c+17}" font-size="7.5" fill="${o.dim}">${m("walsh.virtual")}</text>`,C&&z>=0&&z<n.length){const t=n[z],s=M(t.angle),a=V()?"#ffd700":"#cc8800";e+=`<line x1="${s}" y1="${g}" x2="${s}" y2="${g+w}" stroke="${a}" stroke-width="1" stroke-dasharray="3,2" opacity="0.6"/>`;for(let u=p;u<y;u++){if(u>=t.orbitalEnergies.length)continue;const ot=u-p,nt=l[ot%l.length],G=F(t.orbitalEnergies[u]);e+=`<circle cx="${s}" cy="${G}" r="4.5" fill="${nt}" stroke="#fff" stroke-width="1"/>`}}e+=`<text x="${i+b/2}" y="20" text-anchor="middle" font-size="12" font-weight="600" fill="${o.titleSvg}">${m("walsh.graphTitle")}</text>`,e+=`<text x="${i+b/2}" y="${T-4}" text-anchor="middle" font-size="10" fill="${o.dim}">${m("walsh.xAngle")}</text>`,e+=`<text x="14" y="${g+w/2}" text-anchor="middle" font-size="10" fill="${o.dim}" transform="rotate(-90,14,${g+w/2})">${m("walsh.yOrbEnergy")}</text>`,e+="</svg>";const f=140,x=16,k=f-x-36,q=n.map(t=>t.totalEnergy),D=Math.min(...q),U=Math.max(...q),Z=Math.max((U-D)*.15,.002),I=D-Z,J=U+Z,A=t=>x+k-(t-I)/(J-I)*k;let v=`<svg width="${$}" height="${f}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${$} ${f}" style="max-width:100%;margin-top:8px;">`;v+=`<rect x="${i}" y="${x}" width="${b}" height="${k}" fill="${o.surface}" rx="2"/>`;for(let t=0;t<=3;t++){const s=I+(J-I)*t/3,a=A(s);v+=`<line x1="${i}" y1="${a}" x2="${i+b}" y2="${a}" stroke="${o.grid}" stroke-width="0.5"/>`,v+=`<text x="${i-6}" y="${a+3}" text-anchor="end" font-size="8" fill="${o.dim}">${s.toFixed(4)}</text>`}v+=`<line x1="${i}" y1="${x}" x2="${i}" y2="${x+k}" stroke="${o.axis}" stroke-width="1"/>`,v+=`<line x1="${i}" y1="${x+k}" x2="${i+b}" y2="${x+k}" stroke="${o.axis}" stroke-width="1"/>`;for(let t=d.angleMin;t<=d.angleMax+.1;t+=20){const s=Math.round(t),a=M(s);v+=`<line x1="${a}" y1="${x+k}" x2="${a}" y2="${x+k+4}" stroke="${o.axis}" stroke-width="1"/>`,v+=`<text x="${a}" y="${x+k+14}" text-anchor="middle" font-size="8" fill="${o.dim}">${s}°</text>`}const j=V()?"#00d4ff":"#0077cc";if(n.length>=2){let t="";for(let s=0;s<n.length;s++){const a=M(n[s].angle),u=A(n[s].totalEnergy);t+=s===0?`M${a.toFixed(1)},${u.toFixed(1)}`:` L${a.toFixed(1)},${u.toFixed(1)}`}v+=`<path d="${t}" fill="none" stroke="${j}" stroke-width="2"/>`}for(const t of n)v+=`<circle cx="${M(t.angle).toFixed(1)}" cy="${A(t.totalEnergy).toFixed(1)}" r="2.5" fill="${j}"/>`;if(n.length>0){let t=0;for(let u=1;u<n.length;u++)n[u].totalEnergy<n[t].totalEnergy&&(t=u);const s=M(n[t].angle),a=A(n[t].totalEnergy);v+=`<circle cx="${s}" cy="${a}" r="5" fill="none" stroke="${j}" stroke-width="1.5"/>`,v+=`<text x="${s}" y="${a-8}" text-anchor="middle" font-size="8" fill="${j}">${n[t].angle.toFixed(0)}°</text>`}if(C&&z>=0&&z<n.length){const t=n[z],s=M(t.angle),a=A(t.totalEnergy),u=V()?"#ffd700":"#cc8800";v+=`<line x1="${s}" y1="${x}" x2="${s}" y2="${x+k}" stroke="${u}" stroke-width="1" stroke-dasharray="3,2" opacity="0.6"/>`,v+=`<circle cx="${s}" cy="${a}" r="5" fill="${j}" stroke="${u}" stroke-width="2"/>`}v+=`<text x="${i+b/2}" y="10" text-anchor="middle" font-size="10" font-weight="600" fill="${o.titleSvg}">${m("walsh.totalEnergy")}</text>`,v+=`<text x="${i+b/2}" y="${f-2}" text-anchor="middle" font-size="9" fill="${o.dim}">${m("walsh.xAngle")}</text>`,v+=`<text x="14" y="${x+k/2}" text-anchor="middle" font-size="9" fill="${o.dim}" transform="rotate(-90,14,${x+k/2})">${m("walsh.yTotalEnergy")}</text>`,v+="</svg>",r.innerHTML=e+v}function Ct(){const r=E.querySelector("#slider-area");!r||n.length===0||(z>=n.length&&(z=0),r.innerHTML=`
    <div class="walsh-slider">
      <input type="range" id="angle-slider" min="0" max="${n.length-1}" value="${z}" />
      <div class="walsh-slider-label">${n[z].angle.toFixed(0)}°</div>
    </div>`,E.querySelector("#angle-slider")?.addEventListener("input",o=>{z=parseInt(o.target.value);const l=E.querySelector(".walsh-slider-label");l&&(l.textContent=`${n[z].angle.toFixed(0)}°`),it(),pt()}))}const Ht={H:"#999",Be:"#6C0",C:"#555",O:"#F22"};function pt(){const r=E.querySelector("#mol-vis");if(!r||!C||n.length===0)return;const o=ct(),l=Y,d=n[z].angle,p=220,y=180,$=30,i=l.buildXYZ(d).split(`
`),h=[];for(let e=2;e<i.length;e++){const c=i[e].trim().split(/\s+/);c.length>=4&&h.push({sym:c[0],x:parseFloat(c[1]),z:parseFloat(c[3])})}let g=0,L=0,b=0,w=0;for(const e of[l.angleMin,l.angleMax]){const f=l.buildXYZ(e).split(`
`);for(let x=2;x<f.length;x++){const B=f[x].trim().split(/\s+/),k=parseFloat(B[1]),q=parseFloat(B[3]);g=Math.min(g,k),L=Math.max(L,k),b=Math.min(b,q),w=Math.max(w,q)}}const R=L-g||.01,N=w-b||.01,P=Math.min((p-2*$)/R,(y-2*$)/N),W=(g+L)/2,X=(b+w)/2,H=e=>p/2+(e-W)*P,M=e=>y/2-(e-X)*P;let F=`<svg width="${p}" height="${y}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${p} ${y}" style="max-width:100%;">`;for(let e=0;e<h.length;e++)for(let c=e+1;c<h.length;c++){const f=h[c].x-h[e].x,x=h[c].z-h[e].z;Math.sqrt(f*f+x*x)<2.5&&(F+=`<line x1="${H(h[e].x)}" y1="${M(h[e].z)}" x2="${H(h[c].x)}" y2="${M(h[c].z)}" stroke="${o.grid}" stroke-width="4" stroke-linecap="round"/>`)}for(const e of h){const c=H(e.x),f=M(e.z),x=e.sym==="H"?14:18;F+=`<circle cx="${c}" cy="${f}" r="${x}" fill="${Ht[e.sym]||"#888"}" stroke="${o.axis}" stroke-width="1.5"/>`,F+=`<text x="${c}" y="${f+4}" text-anchor="middle" font-size="11" font-weight="bold" fill="#fff">${e.sym}</text>`}if(h.length>=3){const e=H(h[0].x),c=M(h[0].z),f=28,x=Math.atan2(M(h[1].z)-c,H(h[1].x)-e),B=Math.atan2(M(h[2].z)-c,H(h[2].x)-e),k=e+f*Math.cos(x),q=c+f*Math.sin(x),D=e+f*Math.cos(B),U=c+f*Math.sin(B),Z=V()?"#ffd700":"#cc8800";F+=`<path d="M ${k.toFixed(1)} ${q.toFixed(1)} A ${f} ${f} 0 0 0 ${D.toFixed(1)} ${U.toFixed(1)}" fill="none" stroke="${Z}" stroke-width="1.5"/>`;const I=(x+B)/2,J=e+(f+14)*Math.cos(I),A=c+(f+14)*Math.sin(I);F+=`<text x="${J.toFixed(1)}" y="${A.toFixed(1)}" text-anchor="middle" font-size="11" font-weight="600" fill="${Z}">${d.toFixed(0)}°</text>`}F+="</svg>",r.innerHTML=F}let lt=!1;function Lt(){if(lt)return;lt=!0;const r=document.createElement("style");r.textContent=`
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
  `,document.head.appendChild(r)}gt();$t();K();
