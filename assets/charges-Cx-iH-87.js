import{i as O,c as A,r as j,t as s,a as I,b as G,g as P}from"./nav-Lc1uTV-6.js";import{p as U}from"./parseXYZ-DTmoiCYR.js";import{t as Y,u as _,b as D,B as K,d as V}from"./builder-B-MUn5dT.js";import{a as Q,b as W}from"./properties-CNs3u-bv.js";import{a as X}from"./moleculeViewer3D-B1aGFjc3.js";const z=[{id:"h2",labelKey:"chg.scenH2",descKey:"chg.descH2",category:"homo",charge:0,mult:1,xyz:`2
H2
H  0.0  0.0  0.0
H  0.0  0.0  0.740000`},{id:"n2",labelKey:"chg.scenN2",descKey:"chg.descN2",category:"homo",charge:0,mult:1,xyz:`2
N2
N  0.0  0.0  0.0
N  0.0  0.0  1.098000`},{id:"f2",labelKey:"chg.scenF2",descKey:"chg.descF2",category:"homo",charge:0,mult:1,xyz:`2
F2
F  0.0  0.0  0.0
F  0.0  0.0  1.412000`},{id:"lih",labelKey:"chg.scenLiH",descKey:"chg.descLiH",category:"hetero",charge:0,mult:1,xyz:`2
LiH
Li  0.0  0.0  0.0
H   0.0  0.0  1.596000`},{id:"lif",labelKey:"chg.scenLiF",descKey:"chg.descLiF",category:"hetero",charge:0,mult:1,xyz:`2
LiF
Li  0.0  0.0  0.0
F   0.0  0.0  1.564000`},{id:"hf",labelKey:"chg.scenHF",descKey:"chg.descHF",category:"hetero",charge:0,mult:1,xyz:`2
HF
H  0.0  0.0  0.0
F  0.0  0.0  0.917000`},{id:"co",labelKey:"chg.scenCO",descKey:"chg.descCO",category:"hetero",charge:0,mult:1,xyz:`2
CO
C  0.0  0.0  0.0
O  0.0  0.0  1.128000`},{id:"h2o",labelKey:"chg.scenH2O",descKey:"chg.descH2O",category:"poly",charge:0,mult:1,xyz:(()=>{const e=104*Math.PI/180,o=e/2,r=.96*Math.sin(o),n=.96*Math.cos(o);return`3
H2O
O  0.0  0.0  0.0
H  ${r.toFixed(6)}  0.0  ${n.toFixed(6)}
H  ${(-r).toFixed(6)}  0.0  ${n.toFixed(6)}`})()},{id:"nh3",labelKey:"chg.scenNH3",descKey:"chg.descNH3",category:"poly",charge:0,mult:1,xyz:(()=>{const o=Math.sqrt(.878983),r=Math.sqrt(3)/2;return["4","NH3",`N   0.000000  0.000000  ${.381.toFixed(6)}`,`H   ${o.toFixed(6)}  0.000000  0.000000`,`H   ${(-o/2).toFixed(6)}  ${(o*r).toFixed(6)}  0.000000`,`H   ${(-o/2).toFixed(6)}  ${(-o*r).toFixed(6)}  0.000000`].join(`
`)})()},{id:"ch4",labelKey:"chg.scenCH4",descKey:"chg.descCH4",category:"poly",charge:0,mult:1,xyz:(()=>{const e=1.089/Math.sqrt(3);return["5","CH4","C   0.000000  0.000000  0.000000",`H   ${e.toFixed(6)}  ${e.toFixed(6)}  ${e.toFixed(6)}`,`H   ${e.toFixed(6)}  ${(-e).toFixed(6)}  ${(-e).toFixed(6)}`,`H   ${(-e).toFixed(6)}  ${e.toFixed(6)}  ${(-e).toFixed(6)}`,`H   ${(-e).toFixed(6)}  ${(-e).toFixed(6)}  ${e.toFixed(6)}`].join(`
`)})()},{id:"bh3",labelKey:"chg.scenBH3",descKey:"chg.descBH3",category:"poly",charge:0,mult:1,xyz:`4
BH3
B  0.0  0.0  0.0
${[0,120,240].map(r=>r*Math.PI/180).map(r=>`H  ${(1.19*Math.cos(r)).toFixed(6)}  ${(1.19*Math.sin(r)).toFixed(6)}  0.0`).join(`
`)}`},{id:"heh",labelKey:"chg.scenHeH",descKey:"chg.descHeH",category:"ion",charge:1,mult:1,xyz:`2
HeH+
He  0.0  0.0  0.0
H   0.0  0.0  0.774000`},{id:"bh4",labelKey:"chg.scenBH4",descKey:"chg.descBH4",category:"ion",charge:-1,mult:1,xyz:(()=>{const e=1.255/Math.sqrt(3);return["5","BH4-","B   0.000000  0.000000  0.000000",`H   ${e.toFixed(6)}  ${e.toFixed(6)}  ${e.toFixed(6)}`,`H   ${e.toFixed(6)}  ${(-e).toFixed(6)}  ${(-e).toFixed(6)}`,`H   ${(-e).toFixed(6)}  ${e.toFixed(6)}  ${(-e).toFixed(6)}`,`H   ${(-e).toFixed(6)}  ${(-e).toFixed(6)}  ${e.toFixed(6)}`].join(`
`)})()}],Z={homo:"chg.catHomo",hetero:"chg.catHetero",poly:"chg.catPoly",ion:"chg.catIon"};let F=z[0],y=new Map,h=!1,C=!1;const N=new Map;async function J(t){const e=N.get(t);if(e)return e;const o=`/GANSU-Lite/basis/${t.toLowerCase()}.gbs`,r=await fetch(o);if(!r.ok)throw new Error(`Failed to load basis set: ${t}`);const n=await r.text(),a=V.fromGBS(n);return N.set(t,a),a}const g=document.getElementById("app");function v(){const t=y.has(F.id);g.innerHTML=`
    <div class="opt-page">
      ${j("charges")}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${s("chg.molecule")}</h2>
          <div class="opt-scenario-grid" id="scen-grid"></div>

          <div class="ladder-basis-info">
            <span style="color:var(--color-text-dim);font-size:0.72rem">${s("chg.basis")}: STO-3G</span>
          </div>

          <button id="run-btn" class="opt-run-btn" ${h?"disabled":""}>
            ${h?s("chg.running"):t?s("chg.rerun"):s("chg.run")}
          </button>
          <button id="run-all-btn" class="opt-run-all-btn" ${h?"disabled":""}>
            ${s("chg.runAll")}
          </button>
          ${h?`<button id="stop-btn" class="opt-stop-btn">${s("chg.stop")}</button>`:""}

          <div id="progress-area"></div>
          <div id="summary-area"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="mol-vis">
            ${!t&&!h?`<p class="opt-hint">${s("chg.waiting")}</p>`:""}
          </div>
          <div id="compare-chart"></div>
        </div>
      </div>
    </div>`,ce();const e=g.querySelector("#scen-grid"),o=["homo","hetero","poly","ion"];for(const r of o){const n=z.filter(c=>c.category===r);if(n.length===0)continue;const a=document.createElement("div");a.className="chg-cat-label",a.textContent=s(Z[r]),e.appendChild(a);const m=document.createElement("div");m.className="opt-category-row";for(const c of n){const f=document.createElement("div"),u=y.has(c.id);f.className="opt-scenario-card"+(c.id===F.id?" selected":"")+(u?" computed":""),f.innerHTML=`<strong>${s(c.labelKey)}</strong><span class="conv-desc">${s(c.descKey)}</span>`,f.addEventListener("click",()=>{h||(F=c,v())}),m.appendChild(f)}e.appendChild(m)}g.querySelector("#nav-theme").addEventListener("click",()=>{I(),v()}),g.querySelector("#nav-lang").addEventListener("click",()=>{G(),v()}),g.querySelector("#run-btn").addEventListener("click",()=>{h||ee(F)}),g.querySelector("#run-all-btn").addEventListener("click",()=>{h||te()}),g.querySelector("#stop-btn")?.addEventListener("click",()=>{C=!0}),t&&(ne(y.get(F.id)),re(y.get(F.id))),y.size>1&&ae()}async function E(t){const e=await J("STO-3G"),o=U(t.xyz),r=Math.floor((t.mult-1)/2),n=new Y(o,e,t.charge,r),a=performance.now(),m=_(n,"RHF");await m.solve({eriBackend:"js"});const c=performance.now()-a,f=Q(m.density,m.overlap,n.atoms,n.atomToBasisRange),u=W(m.density,n.atoms,n.primitiveShells,n.cgtoNormalizationFactors,n.numBasis),p=[],b=[];for(const x of n.atoms)p.push(D(x.atomicNumber)),b.push({x:x.coordinate.x*K,y:x.coordinate.y*K,z:x.coordinate.z*K});return{molId:t.id,charges:Array.from(f),dipole:u,atomSymbols:p,atomCoords:b,timeMs:c}}async function ee(t){h=!0,v();try{const e=await E(t);y.set(t.id,e)}catch(e){console.error("Charge computation error:",e)}h=!1,v()}async function te(){h=!0,C=!1,v();for(let t=0;t<z.length&&!C;t++){const e=z[t];oe(t,z.length,s(e.labelKey));try{const o=await E(e);y.set(e.id,o)}catch(o){console.error(`Error computing ${e.id}:`,o)}await new Promise(o=>setTimeout(o,0))}h=!1,v()}function oe(t,e,o){const r=g.querySelector("#progress-area");if(!r)return;const n=(t/e*100).toFixed(0);r.innerHTML=`
    <div class="opt-progress"><div class="opt-progress-bar" style="width:${n}%"></div></div>
    <div class="opt-progress-text">${o} (${t+1}/${e})</div>`}function re(t){const e=g.querySelector("#summary-area");if(!e)return;let o=`<div class="opt-summary">
    <h3>${s("chg.done")}</h3>
    <table>
      <tr>
        <th>${s("chg.colAtom")}</th>
        <th>${s("chg.colCharge")}</th>
      </tr>`;for(let r=0;r<t.atomSymbols.length;r++){const n=t.charges[r],a=L(n);o+=`<tr>
      <td>${t.atomSymbols[r]}${r+1}</td>
      <td style="font-family:monospace;color:${a};font-weight:600">${n>=0?"+":""}${n.toFixed(4)}</td>
    </tr>`}o+=`</table>
    <div style="margin-top:8px;font-size:0.75rem;">
      <span style="color:var(--color-text-dim)">${s("chg.dipole")}:</span>
      <strong>${t.dipole.debye.toFixed(3)} D</strong>
    </div>
    <div style="font-size:0.68rem;color:var(--color-text-dim);margin-top:2px">
      ${s("chg.time")}: ${t.timeMs<1e3?t.timeMs.toFixed(0)+"ms":(t.timeMs/1e3).toFixed(1)+"s"}
    </div>
  </div>`,e.innerHTML=o}function L(t){const o=Math.max(-.6,Math.min(.6,t))/.6;if(o>0){const r=Math.round(255*(1-o)),n=Math.round(255*(1-o));return`rgb(${r},${n},255)`}else{const r=-o,n=Math.round(255*(1-r)),a=Math.round(255*(1-r));return`rgb(255,${n},${a})`}}function ne(t){const e=g.querySelector("#mol-vis");e&&X(e,t)}function ae(){const t=g.querySelector("#compare-chart");if(!t||y.size<2)return;const e=P(),o=[];for(const l of z){const i=y.get(l.id);if(!i)continue;let $=0;for(let w=1;w<i.atomSymbols.length;w++)i.atomSymbols[w]!=="H"&&i.atomSymbols[$]==="H"&&($=w);const M=i.atomSymbols[$];o.push({label:`${s(l.labelKey)} (${M})`,charge:i.charges[$]})}if(o.length<2)return;const r=560,n=200,a=82,m=20,c=28,f=60,u=r-a-m,p=n-c-f,b=Math.max(.1,...o.map(l=>Math.abs(l.charge)))*1.2,x=u/o.length,S=l=>c+p/2-l/b*(p/2);let d=`<svg width="${r}" height="${n}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${r} ${n}" style="max-width:100%;margin-top:12px;">`;d+=`<rect x="${a}" y="${c}" width="${u}" height="${p}" fill="${e.surface}" rx="2"/>`;const H=S(0);d+=`<line x1="${a}" y1="${H}" x2="${a+u}" y2="${H}" stroke="${e.axis}" stroke-width="1"/>`;const T=[b*.5,-b*.5,b,-b];for(const l of T){const i=S(l);i>=c&&i<=c+p&&(d+=`<line x1="${a}" y1="${i}" x2="${a+u}" y2="${i}" stroke="${e.grid}" stroke-width="0.5"/>`,d+=`<text x="${a-6}" y="${i+3}" text-anchor="end" font-size="8" fill="${e.dim}">${l>=0?"+":""}${l.toFixed(2)}</text>`)}d+=`<text x="${a-6}" y="${H+3}" text-anchor="end" font-size="8" fill="${e.dim}">0</text>`,d+=`<line x1="${a}" y1="${c}" x2="${a}" y2="${c+p}" stroke="${e.axis}" stroke-width="1"/>`;for(let l=0;l<o.length;l++){const i=o[l],$=a+x*l+x/2,M=a+x*l+x*.2,w=x*.6,k=S(i.charge),R=L(i.charge);i.charge>=0?d+=`<rect x="${M}" y="${k}" width="${w}" height="${H-k}" fill="${R}" opacity="0.85" rx="2"/>`:d+=`<rect x="${M}" y="${H}" width="${w}" height="${k-H}" fill="${R}" opacity="0.85" rx="2"/>`;const q=i.charge>=0?k-5:k+12;d+=`<text x="${$}" y="${q}" text-anchor="middle" font-size="8" fill="${L(i.charge)}" font-weight="700">${i.charge>=0?"+":""}${i.charge.toFixed(3)}</text>`,d+=`<text x="${$}" y="${c+p+12}" text-anchor="middle" font-size="7.5" fill="${e.dim}" transform="rotate(-30,${$},${c+p+12})">${i.label}</text>`}d+=`<text x="${a+u/2}" y="16" text-anchor="middle" font-size="11" font-weight="600" fill="${e.titleSvg}">${s("chg.compareTitle")}</text>`,d+=`<text x="14" y="${c+p/2}" text-anchor="middle" font-size="9" fill="${e.dim}" transform="rotate(-90,14,${c+p/2})">${s("chg.yCharge")}</text>`,d+="</svg>",t.innerHTML=d}let B=!1;function ce(){if(B)return;B=!0;const t=document.createElement("style");t.textContent=`
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      background: var(--color-bg);
      color: var(--color-text);
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
      min-height: 100vh;
    }
    .opt-page { max-width: 1000px; margin: 0 auto; padding: 16px 20px; }

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
    .opt-category-row { display: flex; flex-wrap: wrap; gap: 5px; margin-bottom: 6px; }
    .opt-scenario-card {
      padding: 5px 10px; border: 1px solid var(--color-border); border-radius: 6px;
      cursor: pointer; transition: all 0.15s; flex: 0 0 auto;
    }
    .opt-scenario-card:hover { background: var(--color-surface-alt); }
    .opt-scenario-card.selected {
      border-color: var(--color-accent); background: var(--color-surface-alt);
      box-shadow: 0 0 0 1px var(--color-accent);
    }
    .opt-scenario-card.computed { border-left: 3px solid var(--color-converged); }
    .opt-scenario-card strong { display: block; font-size: 0.78rem; white-space: nowrap; }
    .conv-desc { display: block !important; font-size: 0.68rem; color: var(--color-text-dim); }

    .chg-cat-label {
      font-size: 0.68rem; font-weight: 600; color: var(--color-text-dim);
      text-transform: uppercase; letter-spacing: 0.04em;
      margin-top: 4px;
    }

    .ladder-basis-info { margin-top: 10px; }

    .opt-run-btn {
      width: 100%; margin-top: 14px; padding: 10px; border: none; border-radius: 8px;
      font-size: 0.85rem; font-weight: 600; cursor: pointer;
      background: var(--color-accent); color: var(--color-accent-on);
      transition: background 0.15s;
    }
    .opt-run-btn:hover:not([disabled]) { background: var(--color-accent-hover); }
    .opt-run-btn[disabled] { opacity: 0.6; cursor: not-allowed; }
    .opt-run-all-btn {
      width: 100%; margin-top: 6px; padding: 8px; border: 1px solid var(--color-accent); border-radius: 8px;
      font-size: 0.8rem; font-weight: 600; cursor: pointer;
      background: none; color: var(--color-accent);
    }
    .opt-run-all-btn:hover:not([disabled]) { background: var(--color-surface-alt); }
    .opt-run-all-btn[disabled] { opacity: 0.6; cursor: not-allowed; }
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

    #mol-vis { width: 100%; text-align: center; }
    #compare-chart { width: 100%; text-align: center; }
    .opt-hint { color: var(--color-text-dim); font-size: 0.85rem; padding: 60px 20px; }

    @media (max-width: 700px) {
      .opt-content { flex-direction: column; }
      .opt-controls { flex: none; }
    }
  `,document.head.appendChild(t)}O();A();v();
