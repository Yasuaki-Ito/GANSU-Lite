import{i as le,c as pe,r as ge,t as x,a as he,b as xe,d as me,g as ee}from"./nav-BaYBNpRD.js";const te=[{id:"h2",labelKey:"ladder.scenH2",descKey:"ladder.descH2",charge:0,mult:1,xyz:`2
H2
H  0.0  0.0  0.0
H  0.0  0.0  0.740000`},{id:"beh2",labelKey:"ladder.scenBeH2",descKey:"ladder.descBeH2",charge:0,mult:1,xyz:`3
BeH2
Be  0.0  0.0  0.0
H  0.0  0.0  1.334000
H  0.0  0.0  -1.334000`},{id:"lih",labelKey:"ladder.scenLiH",descKey:"ladder.descLiH",charge:0,mult:1,xyz:`2
LiH
Li  0.0  0.0  0.0
H   0.0  0.0  1.596000`},{id:"bh3",labelKey:"ladder.scenBH3",descKey:"ladder.descBH3",charge:0,mult:1,xyz:`4
BH3
B  0.0  0.0  0.0
${[0,120,240].map(r=>r*Math.PI/180).map(r=>`H  ${(1.19*Math.cos(r)).toFixed(6)}  ${(1.19*Math.sin(r)).toFixed(6)}  0.0`).join(`
`)}`},{id:"hf",labelKey:"ladder.scenHF",descKey:"ladder.descHF",charge:0,mult:1,xyz:`2
HF
H  0.0  0.0  0.0
F  0.0  0.0  0.917000`},{id:"h2o",labelKey:"ladder.scenH2O",descKey:"ladder.descH2O",charge:0,mult:1,xyz:(()=>{const t=104*Math.PI/180,o=t/2,r=.96*Math.sin(o),n=.96*Math.cos(o);return`3
H2O
O  0.0  0.0  0.0
H  ${r.toFixed(6)}  0.0  ${n.toFixed(6)}
H  ${(-r).toFixed(6)}  0.0  ${n.toFixed(6)}`})()},{id:"nh3",labelKey:"ladder.scenNH3",descKey:"ladder.descNH3",charge:0,mult:1,xyz:(()=>{const o=Math.sqrt(.878983),r=Math.sqrt(3)/2;return["4","NH3",`N   0.000000  0.000000  ${.381.toFixed(6)}`,`H   ${o.toFixed(6)}  0.000000  0.000000`,`H   ${(-o/2).toFixed(6)}  ${(o*r).toFixed(6)}  0.000000`,`H   ${(-o/2).toFixed(6)}  ${(-o*r).toFixed(6)}  0.000000`].join(`
`)})()},{id:"ch4",labelKey:"ladder.scenCH4",descKey:"ladder.descCH4",charge:0,mult:1,xyz:(()=>{const t=1.089/Math.sqrt(3);return["5","CH4","C   0.000000  0.000000  0.000000",`H   ${t.toFixed(6)}  ${t.toFixed(6)}  ${t.toFixed(6)}`,`H   ${t.toFixed(6)}  ${(-t).toFixed(6)}  ${(-t).toFixed(6)}`,`H   ${(-t).toFixed(6)}  ${t.toFixed(6)}  ${(-t).toFixed(6)}`,`H   ${(-t).toFixed(6)}  ${(-t).toFixed(6)}  ${t.toFixed(6)}`].join(`
`)})()},{id:"co",labelKey:"ladder.scenCO",descKey:"ladder.descCO",charge:0,mult:1,xyz:`2
CO
C  0.0  0.0  0.0
O  0.0  0.0  1.128000`},{id:"n2",labelKey:"ladder.scenN2",descKey:"ladder.descN2",charge:0,mult:1,xyz:`2
N2
N  0.0  0.0  0.0
N  0.0  0.0  1.098000`},{id:"f2",labelKey:"ladder.scenF2",descKey:"ladder.descF2",charge:0,mult:1,xyz:`2
F2
F  0.0  0.0  0.0
F  0.0  0.0  1.412000`}],fe=["STO-3G","STO-6G","3-21G","6-31G","6-31G(d,p)","cc-pVDZ","cc-pVTZ"],H=[{id:"HF",label:"HF",category:"Post-HF",type:"hf"},{id:"MP2",label:"MP2",category:"Post-HF",type:"postHF"},{id:"MP3",label:"MP3",category:"Post-HF",type:"postHF"},{id:"CCSD",label:"CCSD",category:"Post-HF",type:"postHF"},{id:"SVWN",label:"SVWN",category:"LDA",type:"dft",functional:"SVWN"},{id:"BLYP",label:"BLYP",category:"GGA",type:"dft",functional:"BLYP"},{id:"PBE",label:"PBE",category:"GGA",type:"dft",functional:"PBE"},{id:"B3LYP",label:"B3LYP",category:"Hybrid",type:"dft",functional:"B3LYP"}];let W=te[0],I="STO-3G",T=new Set(H.map(e=>e.id)),$=[],C=!1,K=!1,G=[],O=[],V=[],Y=0,q=null,_=0,D="none";const J=new Map;async function ue(e){const t=J.get(e);if(t)return t;const o=`/GANSU-Lite/basis/${e.toLowerCase()}.gbs`,r=await fetch(o);if(!r.ok)throw new Error(`Failed to load basis set: ${e}`);const n=await r.text();return J.set(e,n),n}const ye={HF:"#1f77b4",MP2:"#2ca02c",MP3:"#ff7f0e",CCSD:"#d62728",SVWN:"#9467bd",BLYP:"#8c564b",PBE:"#e377c2",B3LYP:"#17becf"},$e={HF:"#4dabf7",MP2:"#51cf66",MP3:"#ffa94d",CCSD:"#ff6b6b",SVWN:"#b197fc",BLYP:"#c4956a",PBE:"#f59fda",B3LYP:"#66d9e8"};function R(e){return me()?$e[e]??"#ccc":ye[e]??"#666"}function be(){re(),G=H.map(o=>T.has(o.id)?"pending":"skipped"),O=H.map(()=>0),V=H.map(o=>T.has(o.id)?"":"skipped"),_=performance.now(),Y=0,D="none";const e=document.createElement("div");e.id="acc-progress-overlay",e.className="acc-prog-overlay";const t=document.createElement("div");t.className="acc-prog-modal",t.id="acc-progress-modal",e.appendChild(t),document.body.appendChild(e),q=setInterval(()=>{Y=performance.now()-_;const o=document.getElementById("acc-prog-elapsed");o&&(o.textContent=(Y/1e3).toFixed(1)+"s")},200),j()}function j(){const e=document.getElementById("acc-progress-modal");if(!e)return;const o=(G.filter(s=>s==="done"||s==="skipped").length/H.length*100).toFixed(0),r=(Y/1e3).toFixed(1),n=D==="wasm-simd"?"active":"inactive",d=D==="wasm"?"active":"inactive",c=D==="js"?"active":"inactive";let m=`
    <div class="acc-prog-header">
      <span class="acc-prog-title">${x("acc.running")}</span>
      <div class="pt-badges">
        <span class="pt-badge simd ${n}">SIMD</span>
        <span class="pt-badge ${d}">WASM</span>
        <span class="pt-badge ${c}">JS</span>
      </div>
      <span class="acc-prog-elapsed" id="acc-prog-elapsed">${r}s</span>
    </div>
    <div class="acc-prog-bar-wrap">
      <div class="acc-prog-bar" style="width:${o}%"></div>
    </div>
    <div class="acc-prog-list">`;for(let s=0;s<H.length;s++){const a=H[s],i=G[s],f=R(a.id);let y,v,w;if(i==="done"){const F=O[s]<1e3?`${O[s].toFixed(0)}ms`:`${(O[s]/1e3).toFixed(1)}s`;y='<span class="acc-prog-check">&#10003;</span>',v=F,w="acc-prog-done"}else i==="running"?(y='<span class="acc-prog-spinner"></span>',v=V[s]||"...",w="acc-prog-active"):i==="skipped"?(y='<span class="acc-prog-skip">&#8212;</span>',v="stopped",w="acc-prog-skipped"):(y='<span class="acc-prog-pending">&#9675;</span>',v="",w="acc-prog-pending-row");s===4&&(m+='<div class="acc-prog-sep"></div>'),m+=`
      <div class="acc-prog-row ${w}">
        <span class="acc-prog-icon">${y}</span>
        <span class="acc-prog-label" style="color:${f}">${a.label}</span>
        <span class="acc-prog-cat">${a.category}</span>
        <span class="acc-prog-status">${v}</span>
      </div>`}m+=`</div>
    <button id="acc-prog-cancel" class="acc-prog-cancel-btn">${x("acc.stop")}</button>`,e.innerHTML=m,document.getElementById("acc-prog-cancel")?.addEventListener("click",()=>{N&&(N.terminate(),N=null),He(0),ae()})}function oe(e,t,o){G[e]=t,o!==void 0&&(V[e]=o),j()}function ve(e,t){O[e]=t}function we(e,t){let o=t;t.includes("Converged")?o=t:t.includes("ERIs")?o="ERI "+(t.match(/(\d+)%/)?.[0]??""):t.includes("nuclear")?o="V_nn":t.includes("one-electron")?o="1e integrals":t.includes("two-electron")?o="ERI":t.includes("transformation")?o="S^{-1/2}":t.includes("SAD")?o="SAD guess":t.includes("GWH")?o="GWH guess":t.includes("initial Fock")?o="H_core guess":t.includes("grid")?o="grid switch":t.includes("SCF did not")?o="not converged":t.includes("SCF iteration")&&(o=t.replace(/^.*?(SCF iteration)/,"$1")),V[e]=o,j()}function re(){q&&(clearInterval(q),q=null),document.getElementById("acc-progress-overlay")?.remove()}const E=document.getElementById("app");function B(){const e=W.id==="h2";!e&&I==="cc-pVTZ"&&(I="cc-pVDZ");const t=fe.map(a=>{const i=a==="cc-pVTZ"&&!e?"disabled":"";return`<option value="${a}" ${a===I?"selected":""} ${i}>${a}${i?" (H₂ only)":""}</option>`}).join("");E.innerHTML=`
    <div class="opt-page">
      ${ge("accuracy")}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${x("acc.molecule")}</h2>
          <div class="opt-scenario-grid" id="scen-grid"></div>

          <h2>${x("acc.basis")}</h2>
          <select id="basis-select" class="acc-basis-select" ${C?"disabled":""}>
            ${t}
          </select>

          <h2>${x("acc.colMethod")}</h2>
          <div class="acc-method-toggles" id="method-toggles"></div>

          <button id="run-btn" class="opt-run-btn" ${C?"disabled":""}>
            ${C?x("acc.running"):x("acc.run")}
          </button>

          <div id="summary-area"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="graph-container">
            ${!K&&!C?`<p class="opt-hint">${x("acc.waiting")}</p>`:""}
          </div>
        </div>
      </div>
    </div>`,Ce();const o=E.querySelector("#scen-grid"),r=document.createElement("div");r.className="opt-category-row";for(const a of te){const i=document.createElement("div");i.className="opt-scenario-card"+(a.id===W.id?" selected":""),i.innerHTML=`<strong>${x(a.labelKey)}</strong><span class="conv-desc">${x(a.descKey)}</span>`,i.addEventListener("click",()=>{C||(W=a,$=[],K=!1,B())}),r.appendChild(i)}o.appendChild(r),E.querySelector("#basis-select").addEventListener("change",a=>{I=a.target.value,B()});const n=E.querySelector("#method-toggles"),d=document.createElement("div");d.className="acc-method-group",d.innerHTML='<span class="acc-method-group-label">Post-HF</span>';const c=document.createElement("div");c.className="acc-method-group-btns",d.appendChild(c);const m=document.createElement("div");m.className="acc-method-group",m.innerHTML='<span class="acc-method-group-label">DFT</span>';const s=document.createElement("div");s.className="acc-method-group-btns",m.appendChild(s);for(const a of H){const i=document.createElement("button"),f=a.type==="hf",y=T.has(a.id);i.className="acc-method-btn"+(y?" active":"")+(f?" always-on":""),i.textContent=a.label,i.disabled=C||f,i.style.color=y?R(a.id):"",i.title=f?"HF is always required":"",i.addEventListener("click",()=>{C||f||(T.has(a.id)?T.delete(a.id):T.add(a.id),B())}),a.type==="hf"||a.type==="postHF"?c.appendChild(i):s.appendChild(i)}n.appendChild(d),n.appendChild(m),E.querySelector("#nav-theme").addEventListener("click",()=>{he(),B()}),E.querySelector("#nav-lang").addEventListener("click",()=>{xe(),B()}),E.querySelector("#run-btn").addEventListener("click",()=>{C||ke()}),$.length>0&&(ne(),ce())}let N=null;async function ke(){C=!0,$=[],K=!1,B(),be();const e=W;try{const t=await ue(I),o=new Worker(new URL("/GANSU-Lite/assets/accuracyWorker-DECf1emp.js",import.meta.url),{type:"module"});N=o,o.postMessage({type:"run",xyzText:e.xyz,basisGBS:t,charge:e.charge,multiplicity:e.mult,selectedMethods:Array.from(T),baseUrl:"/GANSU-Lite/"}),await new Promise((r,n)=>{o.onmessage=d=>{const c=d.data;switch(c.type){case"method-state":oe(c.index,c.state,c.info);break;case"method-info":we(c.index,c.msg);break;case"method-time":ve(c.index,c.ms);break;case"backend":D=c.backend,j();break;case"method-result":$.push({method:H[c.methodIndex],totalEnergy:c.totalEnergy,diffFromHF:c.diffFromHF,timeMs:c.timeMs}),ne(),ce();break;case"done":o.terminate(),N=null,r();break;case"error":o.terminate(),N=null,n(new Error(c.message));break}},o.onerror=d=>{o.terminate(),N=null,n(new Error(d.message||"Worker failed"))}})}catch(t){t instanceof Error&&t.message==="Cancelled"||console.error("Accuracy computation error:",t)}ae()}function He(e){for(let t=e;t<H.length;t++)(G[t]==="pending"||G[t]==="running")&&oe(t,"skipped")}function ae(){C=!1,K=!0,setTimeout(()=>{re(),B()},600)}function ce(){const e=E.querySelector("#summary-area");if(!e||$.length===0)return;const t=$.reduce((r,n)=>r+n.timeMs,0);let o=`<div class="opt-summary">
    <h3>${K?x("acc.done"):x("acc.running")}</h3>
    <table>
      <tr>
        <th>${x("acc.colMethod")}</th>
        <th>${x("acc.colCategory")}</th>
        <th>${x("acc.colEnergy")}</th>
        <th>${x("acc.colDiff")}</th>
        <th>${x("acc.colTime")}</th>
      </tr>`;for(const r of $){const n=r.totalEnergy.toFixed(6),d=r.method.type==="hf"?"—":`${r.diffFromHF.toFixed(2)} mEh`,c=r.timeMs<1e3?`${r.timeMs.toFixed(0)}ms`:`${(r.timeMs/1e3).toFixed(1)}s`,m=R(r.method.id);o+=`<tr>
      <td><span style="color:${m};font-weight:600">${r.method.label}</span></td>
      <td style="font-size:0.65rem;color:var(--color-text-dim)">${r.method.category}</td>
      <td style="font-family:monospace">${n}</td>
      <td style="font-family:monospace">${d}</td>
      <td style="text-align:right">${c}</td>
    </tr>`}o+=`</table>
    ${K?`<div style="margin-top:6px;font-size:0.72rem;color:var(--color-text-dim)">${x("acc.totalTime")}: ${(t/1e3).toFixed(1)}s</div>`:""}
  </div>`,e.innerHTML=o}function ne(){const e=E.querySelector("#graph-container");!e||$.length===0||(e.innerHTML=Fe()+Me())}function Fe(){const e=ee(),t=680,o=300,r=82,n=20,d=36,c=52,m=t-r-n,s=o-d-c,a=$.map(g=>g.totalEnergy),i=Math.min(...a),f=Math.max(...a),y=Math.max((f-i)*.15,.002),v=i-y,w=f+y,F=g=>d+s-(g-v)/(w-v)*s,u=$.length,p=m/H.length;let h=`<svg width="${t}" height="${o}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${t} ${o}" style="max-width:100%;">`;h+=`<rect x="${r}" y="${d}" width="${m}" height="${s}" fill="${e.surface}" rx="2"/>`;for(let g=0;g<=5;g++){const b=v+(w-v)*g/5,S=F(b);h+=`<line x1="${r}" y1="${S}" x2="${r+m}" y2="${S}" stroke="${e.grid}" stroke-width="0.5"/>`,h+=`<text x="${r-6}" y="${S+3}" text-anchor="end" font-size="9" fill="${e.dim}">${b.toFixed(4)}</text>`}const P=r+p*4;h+=`<line x1="${P}" y1="${d}" x2="${P}" y2="${d+s}" stroke="${e.grid}" stroke-width="1" stroke-dasharray="4,3"/>`;const l=r+p*2,k=r+p*6;h+=`<text x="${l}" y="${o-4}" text-anchor="middle" font-size="9" fill="${e.dim}" font-weight="600">Post-HF</text>`,h+=`<text x="${k}" y="${o-4}" text-anchor="middle" font-size="9" fill="${e.dim}" font-weight="600">DFT</text>`,h+=`<line x1="${r}" y1="${d}" x2="${r}" y2="${d+s}" stroke="${e.axis}" stroke-width="1"/>`,h+=`<line x1="${r}" y1="${d+s}" x2="${r+m}" y2="${d+s}" stroke="${e.axis}" stroke-width="1"/>`;const M=$[0].totalEnergy;for(let g=0;g<u;g++){const b=$[g],S=R(b.method.id),z=r+p*g+p/2,L=F(b.totalEnergy);if(b.method.type!=="hf"){const A=F(M),Z=r+p*g+p*.15,se=p*.7,ie=Math.min(L,A),de=Math.abs(A-L);h+=`<rect x="${Z}" y="${ie}" width="${se}" height="${de}" fill="${S}" opacity="0.25" rx="2"/>`}const U=r+p*g+p*.1,X=r+p*g+p*.9;if(h+=`<line x1="${U}" y1="${L}" x2="${X}" y2="${L}" stroke="${S}" stroke-width="3" stroke-linecap="round"/>`,h+=`<text x="${z}" y="${L-8}" text-anchor="middle" font-size="8" fill="${S}" font-weight="600">${b.totalEnergy.toFixed(4)}</text>`,h+=`<text x="${z}" y="${d+s+14}" text-anchor="middle" font-size="9" fill="${e.dim}" font-weight="600">${b.method.label}</text>`,b.method.type==="dft"&&(h+=`<text x="${z}" y="${d+s+25}" text-anchor="middle" font-size="7" fill="${e.dim}" opacity="0.7">${b.method.category}</text>`),g>0&&g!==4){const A=F($[g-1].totalEnergy),Z=r+p*(g-1)+p/2;h+=`<line x1="${Z}" y1="${A}" x2="${z}" y2="${L}" stroke="${e.grid}" stroke-width="1" stroke-dasharray="4,3"/>`}}return h+=`<text x="${r+m/2}" y="20" text-anchor="middle" font-size="12" font-weight="600" fill="${e.titleSvg}">${x("acc.graphTitle")}</text>`,h+=`<text x="14" y="${d+s/2}" text-anchor="middle" font-size="10" fill="${e.dim}" transform="rotate(-90,14,${d+s/2})">${x("acc.yEnergy")}</text>`,h+="</svg>",h}function Me(){const e=ee(),t=$.filter(l=>l.method.type!=="hf");if(t.length<1)return"";const o=680,r=200,n=82,d=20,c=28,m=48,s=o-n-d,a=r-c-m,i=t.length,f=s/i,y=t.map(l=>l.diffFromHF),w=Math.max(Math.abs(Math.min(...y)),Math.abs(Math.max(...y)),1)*1.2,F=l=>c+a/2-l/w*(a/2);let u=`<svg width="${o}" height="${r}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${o} ${r}" style="max-width:100%;margin-top:8px;">`;u+=`<rect x="${n}" y="${c}" width="${s}" height="${a}" fill="${e.surface}" rx="2"/>`;const p=t.filter(l=>l.method.type==="postHF").length;if(p>0&&p<i){const l=n+f*p;u+=`<line x1="${l}" y1="${c}" x2="${l}" y2="${c+a}" stroke="${e.grid}" stroke-width="1" stroke-dasharray="4,3"/>`}const h=4;for(let l=-h;l<=h;l++){const k=w*l/h,M=F(k);M<c||M>c+a||(u+=`<line x1="${n}" y1="${M}" x2="${n+s}" y2="${M}" stroke="${e.grid}" stroke-width="${l===0?1:.5}"/>`,u+=`<text x="${n-6}" y="${M+3}" text-anchor="end" font-size="8" fill="${e.dim}">${k.toFixed(0)}</text>`)}const P=F(0);u+=`<line x1="${n}" y1="${P}" x2="${n+s}" y2="${P}" stroke="${e.axis}" stroke-width="1"/>`,u+=`<line x1="${n}" y1="${c}" x2="${n}" y2="${c+a}" stroke="${e.axis}" stroke-width="1"/>`;for(let l=0;l<i;l++){const k=t[l],M=R(k.method.id),g=n+f*l+f/2,b=n+f*l+f*.2,S=f*.6,z=F(y[l]),L=Math.abs(z-P),U=Math.min(z,P);u+=`<rect x="${b}" y="${U}" width="${S}" height="${L}" fill="${M}" opacity="0.75" rx="2"/>`;const X=y[l]<0?z+12:z-5;u+=`<text x="${g}" y="${X}" text-anchor="middle" font-size="9" fill="${M}" font-weight="700">${y[l].toFixed(1)}</text>`,u+=`<text x="${g}" y="${c+a+14}" text-anchor="middle" font-size="9" fill="${e.dim}" font-weight="600">${k.method.label}</text>`,k.method.type==="dft"&&(u+=`<text x="${g}" y="${c+a+25}" text-anchor="middle" font-size="7" fill="${e.dim}" opacity="0.7">${k.method.category}</text>`)}if(p>0&&p<i){const l=n+f*p/2,k=n+f*(p+(i-p)/2);u+=`<text x="${l}" y="${r-4}" text-anchor="middle" font-size="9" fill="${e.dim}" font-weight="600">Post-HF</text>`,u+=`<text x="${k}" y="${r-4}" text-anchor="middle" font-size="9" fill="${e.dim}" font-weight="600">DFT</text>`}return u+=`<text x="${n+s/2}" y="16" text-anchor="middle" font-size="11" font-weight="600" fill="${e.titleSvg}">${x("acc.diffTitle")}</text>`,u+=`<text x="14" y="${c+a/2}" text-anchor="middle" font-size="9" fill="${e.dim}" transform="rotate(-90,14,${c+a/2})">${x("acc.yDiff")}</text>`,u+="</svg>",u}let Q=!1;function Ce(){if(Q)return;Q=!0;const e=document.createElement("style");e.textContent=`
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      background: var(--color-bg);
      color: var(--color-text);
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
      min-height: 100vh;
    }
    .opt-page { max-width: 1040px; margin: 0 auto; padding: 16px 20px; }

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

    .acc-basis-select {
      width: 100%; padding: 6px 8px; border: 1px solid var(--color-border);
      border-radius: 6px; background: var(--color-surface);
      color: var(--color-text); font-size: 0.82rem;
    }

    .opt-run-btn {
      width: 100%; margin-top: 14px; padding: 10px; border: none; border-radius: 8px;
      font-size: 0.85rem; font-weight: 600; cursor: pointer;
      background: var(--color-accent); color: var(--color-accent-on);
      transition: background 0.15s;
    }
    .opt-run-btn:hover:not([disabled]) { background: var(--color-accent-hover); }
    .opt-run-btn[disabled] { opacity: 0.6; cursor: not-allowed; }

    .opt-summary {
      margin-top: 14px; padding: 12px; background: var(--color-surface-alt);
      border-radius: 8px; border: 1px solid var(--color-border);
    }
    .opt-summary h3 {
      font-size: 0.82rem; color: var(--color-converged); margin-bottom: 8px;
    }
    .opt-summary table { width: 100%; font-size: 0.72rem; border-collapse: collapse; }
    .opt-summary th {
      font-size: 0.65rem; color: var(--color-text-dim); text-align: left;
      padding: 2px 4px; font-weight: 500;
    }
    .opt-summary td { padding: 3px 4px; }

    #graph-container { width: 100%; text-align: center; }
    .opt-hint { color: var(--color-text-dim); font-size: 0.85rem; padding: 60px 20px; }

    /* ── Progress modal ── */
    .acc-prog-overlay {
      position: fixed; inset: 0; z-index: 9000;
      background: rgba(0,0,0,0.45);
      display: flex; align-items: center; justify-content: center;
    }
    .acc-prog-modal {
      background: var(--color-surface, #fff); border: 1px solid var(--color-border, #ddd);
      border-radius: 12px; padding: 20px 24px; width: 380px; max-width: 92vw;
      box-shadow: 0 8px 32px rgba(0,0,0,0.25);
    }
    .acc-prog-header {
      display: flex; justify-content: space-between; align-items: center;
      margin-bottom: 10px;
    }
    .acc-prog-title { font-size: 0.9rem; font-weight: 600; color: var(--color-text); }
    .acc-prog-elapsed { font-size: 0.78rem; color: var(--color-text-dim); font-family: monospace; }

    .acc-prog-bar-wrap {
      height: 6px; background: var(--color-progress-bg, #e0e4ea);
      border-radius: 3px; margin-bottom: 14px; overflow: hidden;
    }
    .acc-prog-bar {
      height: 100%; background: var(--color-accent, #4a90d9);
      border-radius: 3px; transition: width 0.3s ease;
    }

    .acc-prog-list { display: flex; flex-direction: column; gap: 3px; }

    .acc-prog-sep {
      height: 1px; background: var(--color-border, #ddd);
      margin: 4px 0;
    }

    .acc-prog-row {
      display: flex; align-items: center; gap: 8px;
      padding: 4px 6px; border-radius: 5px; font-size: 0.78rem;
    }
    .acc-prog-active {
      background: var(--color-surface-alt, #f5f5f5);
    }
    .acc-prog-icon { width: 18px; text-align: center; flex-shrink: 0; }
    .acc-prog-label { font-weight: 600; width: 50px; }
    .acc-prog-cat { font-size: 0.65rem; color: var(--color-text-dim); width: 52px; }
    .acc-prog-status {
      flex: 1; text-align: right; font-size: 0.7rem;
      color: var(--color-text-dim); font-family: monospace;
      overflow: hidden; text-overflow: ellipsis; white-space: nowrap;
    }

    .acc-prog-check { color: var(--color-converged, #2ecc71); font-weight: 700; }
    .acc-prog-skip { color: var(--color-text-dim, #999); }
    .acc-prog-pending { color: var(--color-text-dim, #ccc); font-size: 0.6rem; }
    .acc-prog-pending-row { opacity: 0.5; }
    .acc-prog-skipped { opacity: 0.4; }

    .acc-prog-spinner {
      display: inline-block; width: 14px; height: 14px;
      border: 2px solid var(--color-border, #ddd);
      border-top-color: var(--color-accent, #4a90d9);
      border-radius: 50%;
      animation: acc-spin 0.8s linear infinite;
    }
    @keyframes acc-spin { to { transform: rotate(360deg); } }

    .acc-prog-cancel-btn {
      display: block; width: 100%; margin-top: 14px;
      padding: 8px; border: 1px solid var(--color-error, #e74c3c);
      border-radius: 8px; font-size: 0.8rem; font-weight: 600;
      cursor: pointer; background: none;
      color: var(--color-error, #e74c3c);
    }
    .acc-prog-cancel-btn:hover { background: rgba(231,76,60,0.08); }

    /* ── Method toggle buttons ── */
    .acc-method-toggles { display: flex; flex-direction: column; gap: 5px; }
    .acc-method-group { display: flex; align-items: center; gap: 4px; }
    .acc-method-group-label {
      font-size: 0.6rem; color: var(--color-text-dim); font-weight: 500;
      letter-spacing: 0.03em; width: 40px; flex-shrink: 0;
    }
    .acc-method-group-btns { display: flex; gap: 3px; flex: 1; }
    .acc-method-btn {
      padding: 3px 0; border: 1px solid var(--color-border);
      border-radius: 4px; font-size: 0.7rem; font-weight: 600;
      cursor: pointer; background: var(--color-surface);
      color: var(--color-text-dim); transition: all 0.15s;
      flex: 1; text-align: center; min-width: 0;
    }
    .acc-method-btn:hover:not([disabled]) { background: var(--color-surface-alt); }
    .acc-method-btn.active {
      border-color: currentColor; background: var(--color-surface-alt);
      box-shadow: 0 0 0 1px currentColor;
    }
    .acc-method-btn.always-on {
      opacity: 0.7; cursor: default;
    }
    .acc-method-btn[disabled]:not(.always-on) { opacity: 0.5; cursor: not-allowed; }

    @media (max-width: 700px) {
      .opt-content { flex-direction: column; }
      .opt-controls { flex: none; }
    }
  `,document.head.appendChild(e)}le();pe();B();
