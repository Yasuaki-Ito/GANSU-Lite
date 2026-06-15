import"./styles-B8WF9G-F.js";import{p as A}from"./parseXYZ-ByfvMHmk.js";import{g as j,b as I,B as K,f as G}from"./ri-DT3jfd0v.js";import{a as P,b as D}from"./theoryControls-CuXHSm5X.js";import{a as U,b as Y}from"./properties-Odo1T7fO.js";import{i as _,c as V,r as Q,t as i,a as W,b as X,g as Z}from"./nav-oJm8MoOS.js";import{a as J}from"./moleculeViewer3D-D3ce842e.js";const z=[{id:"h2",labelKey:"chg.scenH2",descKey:"chg.descH2",category:"homo",charge:0,mult:1,xyz:`2
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
`)})()}],ee={homo:"chg.catHomo",hetero:"chg.catHetero",poly:"chg.catPoly",ion:"chg.catIon"};let F=z[0],f=new Map,g=!1,L=!1,C="HF";const T=new Map;async function te(t){const e=T.get(t);if(e)return e;const o=`/GANSU-Lite/basis/${t.toLowerCase()}.gbs`,r=await fetch(o);if(!r.ok)throw new Error(`Failed to load basis set: ${t}`);const n=await r.text(),a=G.fromGBS(n);return T.set(t,a),a}const m=document.getElementById("app");function v(){const t=f.has(F.id);m.innerHTML=`
    <div class="opt-page">
      ${Q("charges")}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${i("chg.molecule")}</h2>
          <div class="opt-scenario-grid" id="scen-grid"></div>

          <div class="ladder-basis-info">
            <span style="color:var(--color-text-dim);font-size:0.72rem">${i("chg.basis")}: STO-3G</span>
          </div>

          <div class="theory-row" style="margin:8px 0;display:flex;align-items:center;gap:8px;font-size:0.9rem;">
            <span>Theory:</span>${P("theory-sel",C)}
          </div>

          <button id="run-btn" class="opt-run-btn" ${g?"disabled":""}>
            ${g?i("chg.running"):t?i("chg.rerun"):i("chg.run")}
          </button>
          <button id="run-all-btn" class="opt-run-all-btn" ${g?"disabled":""}>
            ${i("chg.runAll")}
          </button>
          ${g?`<button id="stop-btn" class="opt-stop-btn">${i("chg.stop")}</button>`:""}

          <div id="progress-area"></div>
          <div id="summary-area"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="mol-vis">
            ${!t&&!g?`<p class="opt-hint">${i("chg.waiting")}</p>`:""}
          </div>
          <div id="compare-chart"></div>
        </div>
      </div>
    </div>`,se();const e=m.querySelector("#scen-grid"),o=["homo","hetero","poly","ion"];for(const n of o){const a=z.filter(x=>x.category===n);if(a.length===0)continue;const u=document.createElement("div");u.className="chg-cat-label",u.textContent=i(ee[n]),e.appendChild(u);const s=document.createElement("div");s.className="opt-category-row";for(const x of a){const h=document.createElement("div"),d=f.has(x.id);h.className="opt-scenario-card"+(x.id===F.id?" selected":"")+(d?" computed":""),h.innerHTML=`<strong>${i(x.labelKey)}</strong><span class="conv-desc">${i(x.descKey)}</span>`,h.addEventListener("click",()=>{g||(F=x,v())}),s.appendChild(h)}e.appendChild(s)}m.querySelector("#nav-theme").addEventListener("click",()=>{W(),v()}),m.querySelector("#nav-lang").addEventListener("click",()=>{X(),v()});const r=m.querySelector("#theory-sel");r&&r.addEventListener("change",()=>{C=r.value}),m.querySelector("#run-btn").addEventListener("click",()=>{g||oe(F)}),m.querySelector("#run-all-btn").addEventListener("click",()=>{g||re()}),m.querySelector("#stop-btn")?.addEventListener("click",()=>{L=!0}),t&&(ce(f.get(F.id)),ae(f.get(F.id))),f.size>1&&ie()}async function E(t){const e=await te("STO-3G"),o=A(t.xyz),r=Math.floor((t.mult-1)/2),n=new j(o,e,t.charge,r),a=performance.now(),u=await D(n,e,C);await u.solve({eriBackend:"js"});const s=performance.now()-a,x=U(u.density,u.overlap,n.atoms,n.atomToBasisRange),h=Y(u.density,n.atoms,n.primitiveShells,n.cgtoNormalizationFactors,n.numBasis),d=[],b=[];for(const y of n.atoms)d.push(I(y.atomicNumber)),b.push({x:y.coordinate.x*K,y:y.coordinate.y*K,z:y.coordinate.z*K});return{molId:t.id,charges:Array.from(x),dipole:h,atomSymbols:d,atomCoords:b,timeMs:s}}async function oe(t){g=!0,v();try{const e=await E(t);f.set(t.id,e)}catch(e){console.error("Charge computation error:",e)}g=!1,v()}async function re(){g=!0,L=!1,v();for(let t=0;t<z.length&&!L;t++){const e=z[t];ne(t,z.length,i(e.labelKey));try{const o=await E(e);f.set(e.id,o)}catch(o){console.error(`Error computing ${e.id}:`,o)}await new Promise(o=>setTimeout(o,0))}g=!1,v()}function ne(t,e,o){const r=m.querySelector("#progress-area");if(!r)return;const n=(t/e*100).toFixed(0);r.innerHTML=`
    <div class="opt-progress"><div class="opt-progress-bar" style="width:${n}%"></div></div>
    <div class="opt-progress-text">${o} (${t+1}/${e})</div>`}function ae(t){const e=m.querySelector("#summary-area");if(!e)return;let o=`<div class="opt-summary">
    <h3>${i("chg.done")}</h3>
    <table>
      <tr>
        <th>${i("chg.colAtom")}</th>
        <th>${i("chg.colCharge")}</th>
      </tr>`;for(let r=0;r<t.atomSymbols.length;r++){const n=t.charges[r],a=N(n);o+=`<tr>
      <td>${t.atomSymbols[r]}${r+1}</td>
      <td style="font-family:monospace;color:${a};font-weight:600">${n>=0?"+":""}${n.toFixed(4)}</td>
    </tr>`}o+=`</table>
    <div style="margin-top:8px;font-size:0.75rem;">
      <span style="color:var(--color-text-dim)">${i("chg.dipole")}:</span>
      <strong>${t.dipole.debye.toFixed(3)} D</strong>
    </div>
    <div style="font-size:0.68rem;color:var(--color-text-dim);margin-top:2px">
      ${i("chg.time")}: ${t.timeMs<1e3?t.timeMs.toFixed(0)+"ms":(t.timeMs/1e3).toFixed(1)+"s"}
    </div>
  </div>`,e.innerHTML=o}function N(t){const o=Math.max(-.6,Math.min(.6,t))/.6;if(o>0){const r=Math.round(255*(1-o)),n=Math.round(255*(1-o));return`rgb(${r},${n},255)`}else{const r=-o,n=Math.round(255*(1-r)),a=Math.round(255*(1-r));return`rgb(255,${n},${a})`}}function ce(t){const e=m.querySelector("#mol-vis");e&&J(e,t)}function ie(){const t=m.querySelector("#compare-chart");if(!t||f.size<2)return;const e=Z(),o=[];for(const l of z){const c=f.get(l.id);if(!c)continue;let $=0;for(let w=1;w<c.atomSymbols.length;w++)c.atomSymbols[w]!=="H"&&c.atomSymbols[$]==="H"&&($=w);const k=c.atomSymbols[$];o.push({label:`${i(l.labelKey)} (${k})`,charge:c.charges[$]})}if(o.length<2)return;const r=560,n=200,a=82,u=20,s=28,x=60,h=r-a-u,d=n-s-x,b=Math.max(.1,...o.map(l=>Math.abs(l.charge)))*1.2,y=h/o.length,S=l=>s+d/2-l/b*(d/2);let p=`<svg width="${r}" height="${n}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${r} ${n}" style="max-width:100%;margin-top:12px;">`;p+=`<rect x="${a}" y="${s}" width="${h}" height="${d}" fill="${e.surface}" rx="2"/>`;const H=S(0);p+=`<line x1="${a}" y1="${H}" x2="${a+h}" y2="${H}" stroke="${e.axis}" stroke-width="1"/>`;const q=[b*.5,-b*.5,b,-b];for(const l of q){const c=S(l);c>=s&&c<=s+d&&(p+=`<line x1="${a}" y1="${c}" x2="${a+h}" y2="${c}" stroke="${e.grid}" stroke-width="0.5"/>`,p+=`<text x="${a-6}" y="${c+3}" text-anchor="end" font-size="8" fill="${e.dim}">${l>=0?"+":""}${l.toFixed(2)}</text>`)}p+=`<text x="${a-6}" y="${H+3}" text-anchor="end" font-size="8" fill="${e.dim}">0</text>`,p+=`<line x1="${a}" y1="${s}" x2="${a}" y2="${s+d}" stroke="${e.axis}" stroke-width="1"/>`;for(let l=0;l<o.length;l++){const c=o[l],$=a+y*l+y/2,k=a+y*l+y*.2,w=y*.6,M=S(c.charge),R=N(c.charge);c.charge>=0?p+=`<rect x="${k}" y="${M}" width="${w}" height="${H-M}" fill="${R}" opacity="0.85" rx="2"/>`:p+=`<rect x="${k}" y="${H}" width="${w}" height="${M-H}" fill="${R}" opacity="0.85" rx="2"/>`;const O=c.charge>=0?M-5:M+12;p+=`<text x="${$}" y="${O}" text-anchor="middle" font-size="8" fill="${N(c.charge)}" font-weight="700">${c.charge>=0?"+":""}${c.charge.toFixed(3)}</text>`,p+=`<text x="${$}" y="${s+d+12}" text-anchor="middle" font-size="7.5" fill="${e.dim}" transform="rotate(-30,${$},${s+d+12})">${c.label}</text>`}p+=`<text x="${a+h/2}" y="16" text-anchor="middle" font-size="11" font-weight="600" fill="${e.titleSvg}">${i("chg.compareTitle")}</text>`,p+=`<text x="14" y="${s+d/2}" text-anchor="middle" font-size="9" fill="${e.dim}" transform="rotate(-90,14,${s+d/2})">${i("chg.yCharge")}</text>`,p+="</svg>",t.innerHTML=p}let B=!1;function se(){if(B)return;B=!0;const t=document.createElement("style");t.textContent=`
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
  `,document.head.appendChild(t)}_();V();v();
