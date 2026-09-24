import"./styles-T9xqG_Ue.js";import{p as nt}from"./parseXYZ-BV3t15l8.js";import{g as it,f as qt,M as Mt}from"./ri-Dh08ihv_.js";import{t as _,a as Zt,b as lt,H as Gt}from"./theoryControls-BlyfxdhA.js";import{t as F,i as jt,c as Dt,r as Ut,a as Vt,b as _t,d as ct,g as Et}from"./nav-BxAkgz0C.js";import{e as zt}from"./radioGroup-C7LT8iiA.js";const Ft={1:.31,2:.28,3:1.28,4:.96,5:.84,6:.76,7:.71,8:.66,9:.57,10:.58,11:1.66,12:1.41,13:1.21,14:1.11,15:1.07,16:1.05,17:1.02,18:1.06,19:2.03,20:1.76,21:1.7,22:1.6,23:1.53,24:1.39,25:1.39,26:1.32,27:1.26,28:1.24,29:1.32,30:1.22,31:1.22,32:1.2,33:1.19,34:1.2,35:1.2,36:1.16},Wt=[[0,0,1,"H"],[0,17,2,"He"],[1,0,3,"Li"],[1,1,4,"Be"],[1,12,5,"B"],[1,13,6,"C"],[1,14,7,"N"],[1,15,8,"O"],[1,16,9,"F"],[1,17,10,"Ne"],[2,0,11,"Na"],[2,1,12,"Mg"],[2,12,13,"Al"],[2,13,14,"Si"],[2,14,15,"P"],[2,15,16,"S"],[2,16,17,"Cl"],[2,17,18,"Ar"],[3,0,19,"K"],[3,1,20,"Ca"],[3,2,21,"Sc"],[3,3,22,"Ti"],[3,4,23,"V"],[3,5,24,"Cr"],[3,6,25,"Mn"],[3,7,26,"Fe"],[3,8,27,"Co"],[3,9,28,"Ni"],[3,10,29,"Cu"],[3,11,30,"Zn"],[3,12,31,"Ga"],[3,13,32,"Ge"],[3,14,33,"As"],[3,15,34,"Se"],[3,16,35,"Br"],[3,17,36,"Kr"]];function wt(t,a,r,o,i=18){const e=document.createElement("div");e.className="pte-grid";for(const[c,l,n,p]of Wt){const x=document.createElement("div");x.className="pte-cell",x.style.gridRow=String(c+1),x.style.gridColumn=String(l+1);const u=n>i;u&&x.classList.add("pte-disabled");const b=n<=2||n<=4||n>=11&&n<=12||n>=19&&n<=20?"pte-s":n>=21&&n<=30?"pte-d":"pte-p";x.classList.add(b),n===o&&x.classList.add("pte-selected"),x.innerHTML=`<span class="pte-z">${n}</span><span class="pte-sym">${p}</span>`,u||x.addEventListener("click",()=>t(p,n)),e.appendChild(x)}return zt(e,{label:a,key:r,itemSelector:".pte-cell",selectedClass:"pte-selected",disabledClass:"pte-disabled",coords:c=>({row:Number(c.style.gridRow)||0,col:Number(c.style.gridColumn)||0})}),e}function Jt(t,a,r){let o=a??null,i=r??null;const e=document.createElement("div");e.className="pte-overlay";const c=document.createElement("div");c.className="pte-modal pte-dual",c.setAttribute("role","dialog"),c.setAttribute("aria-modal","true"),c.setAttribute("aria-label",F("opt.customSelect"));function l(){c.innerHTML="";const b=document.createElement("div");b.className="pte-tables";const h=document.createElement("div");h.className="pte-col";const M=document.createElement("div");M.className="pte-title",M.textContent=`${F("opt.customAtomA")}${o?": "+o.symbol:""}`,h.appendChild(M),h.appendChild(wt((A,L)=>{o={symbol:A,z:L},l()},F("opt.customAtomA"),"pte-atom-a",o?.z,20)),b.appendChild(h);const $=document.createElement("div");$.className="pte-col";const v=document.createElement("div");v.className="pte-title",v.textContent=`${F("opt.customAtomB")}${i?": "+i.symbol:""}`,$.appendChild(v),$.appendChild(wt((A,L)=>{i={symbol:A,z:L},l()},F("opt.customAtomB"),"pte-atom-b",i?.z,12)),b.appendChild($),c.appendChild(b);const S=document.createElement("div");S.className="pte-btn-row";const g=document.createElement("button");g.className="pte-btn pte-btn-cancel",g.textContent="Cancel",g.addEventListener("click",()=>p()),S.appendChild(g);const m=document.createElement("button");m.className="pte-btn pte-btn-ok",m.textContent="OK",m.disabled=!o||!i,m.addEventListener("click",()=>{o&&i&&(p(),t(o.symbol,o.z,i.symbol,i.z))}),S.appendChild(m),c.appendChild(S)}const n=document.activeElement;function p(){e.remove(),document.removeEventListener("keydown",u,!0),n?.focus?.()}function x(){return Array.from(c.querySelectorAll('[tabindex="0"], button:not([disabled])'))}function u(b){if(b.key==="Escape"){b.preventDefault(),p();return}if(b.key!=="Tab")return;const h=x();if(h.length===0)return;const M=h[0],$=h[h.length-1],v=document.activeElement;c.contains(v)?b.shiftKey&&v===M?(b.preventDefault(),$.focus()):!b.shiftKey&&v===$&&(b.preventDefault(),M.focus()):(b.preventDefault(),(b.shiftKey?$:M).focus())}l(),e.appendChild(c),e.addEventListener("click",b=>{b.target===e&&p()}),document.addEventListener("keydown",u,!0),document.body.appendChild(e),x()[0]?.focus()}const Lt=[{id:"h2",category:"dissociation",labelKey:"opt.scenH2",descKey:"opt.descH2",paramType:"bond",defaultMin:.4,defaultMax:3,defaultSteps:20,defaultCharge:0,defaultMult:1,basisOptions:["3-21G","6-31G"],generateXYZ:t=>`2
H2 R=${t.toFixed(3)}
H  0.0  0.0  0.0
H  0.0  0.0  ${t.toFixed(6)}`},{id:"hf",category:"dissociation",labelKey:"opt.scenHF",descKey:"opt.descHF",paramType:"bond",defaultMin:.5,defaultMax:3,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>`2
HF R=${t.toFixed(3)}
H  0.0  0.0  0.0
F  0.0  0.0  ${t.toFixed(6)}`},{id:"lih",category:"dissociation",labelKey:"opt.scenLiH",descKey:"opt.descLiH",paramType:"bond",defaultMin:.8,defaultMax:4,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>`2
LiH R=${t.toFixed(3)}
Li  0.0  0.0  0.0
H   0.0  0.0  ${t.toFixed(6)}`},{id:"n2",category:"dissociation",labelKey:"opt.scenN2",descKey:"opt.descN2",paramType:"bond",defaultMin:.8,defaultMax:3,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>`2
N2 R=${t.toFixed(3)}
N  0.0  0.0  0.0
N  0.0  0.0  ${t.toFixed(6)}`},{id:"f2",category:"dissociation",labelKey:"opt.scenF2",descKey:"opt.descF2",paramType:"bond",defaultMin:.8,defaultMax:3.5,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>`2
F2 R=${t.toFixed(3)}
F  0.0  0.0  0.0
F  0.0  0.0  ${t.toFixed(6)}`},{id:"heh+",category:"dissociation",labelKey:"opt.scenHeH",descKey:"opt.descHeH",paramType:"bond",defaultMin:.5,defaultMax:3,defaultSteps:20,defaultCharge:1,defaultMult:1,generateXYZ:t=>`2
HeH+ R=${t.toFixed(3)}
He  0.0  0.0  0.0
H   0.0  0.0  ${t.toFixed(6)}`},{id:"li2",category:"dissociation",labelKey:"opt.scenLi2",descKey:"opt.descLi2",paramType:"bond",defaultMin:1.5,defaultMax:5,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>`2
Li2 R=${t.toFixed(3)}
Li  0.0  0.0  0.0
Li  0.0  0.0  ${t.toFixed(6)}`},{id:"he2",category:"dissociation",labelKey:"opt.scenHe2",descKey:"opt.descHe2",paramType:"bond",defaultMin:1,defaultMax:5,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>`2
He2 R=${t.toFixed(3)}
He  0.0  0.0  0.0
He  0.0  0.0  ${t.toFixed(6)}`},{id:"c2h2",category:"dissociation",labelKey:"opt.scenC2H2",descKey:"opt.descC2H2",paramType:"bond",defaultMin:.9,defaultMax:2.5,defaultSteps:20,defaultCharge:0,defaultMult:1,basisOptions:["STO-3G"],paramAtomPair:[1,2],generateXYZ:t=>`4
C2H2 CC=${t.toFixed(3)}
H  0.0  0.0  ${(-1.06).toFixed(6)}
C  0.0  0.0  0.0
C  0.0  0.0  ${t.toFixed(6)}
H  0.0  0.0  ${(t+1.06).toFixed(6)}`},{id:"c2h4",category:"dissociation",labelKey:"opt.scenC2H4",descKey:"opt.descC2H4",paramType:"bond",defaultMin:1,defaultMax:2.5,defaultSteps:20,defaultCharge:0,defaultMult:1,basisOptions:["STO-3G"],generateXYZ:t=>{const r=121.7*Math.PI/180,o=1.08*Math.sin(r-Math.PI/2)*-1,i=1.08*Math.cos(r-Math.PI/2),e=Math.abs(o);return["6",`C2H4 CC=${t.toFixed(3)}`,"C   0.000000  0.000000  0.000000",`C   0.000000  0.000000  ${t.toFixed(6)}`,`H   ${e.toFixed(6)}  0.000000  ${(-i).toFixed(6)}`,`H   ${(-e).toFixed(6)}  0.000000  ${(-i).toFixed(6)}`,`H   ${e.toFixed(6)}  0.000000  ${(t+i).toFixed(6)}`,`H   ${(-e).toFixed(6)}  0.000000  ${(t+i).toFixed(6)}`].join(`
`)}},{id:"h2o_bond",category:"dissociation",labelKey:"opt.scenH2Obond",descKey:"opt.descH2Obond",paramType:"bond",defaultMin:.5,defaultMax:3,defaultSteps:20,defaultCharge:0,defaultMult:1,basisOptions:["3-21G","6-31G"],paramAtomPair:[0,2],generateXYZ:t=>{const a=Math.sin(52*Math.PI/180),r=Math.cos(52*Math.PI/180);return["3",`H2O R=${t.toFixed(3)}`,"O   0.000000  0.000000  0.000000",`H   ${(.96*a).toFixed(6)}  0.000000  ${(.96*r).toFixed(6)}`,`H   ${(-t*a).toFixed(6)}  0.000000  ${(t*r).toFixed(6)}`].join(`
`)}},{id:"h2o",category:"geometry",labelKey:"opt.scenH2O",descKey:"opt.descH2O",paramType:"angle",defaultMin:80,defaultMax:180,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>{const o=t*Math.PI/180/2,i=.96*Math.sin(o),e=.96*Math.cos(o);return`3
H2O angle=${t.toFixed(1)}
O  0.0  0.0  0.0
H  ${i.toFixed(6)}  0.0  ${e.toFixed(6)}
H  ${(-i).toFixed(6)}  0.0  ${e.toFixed(6)}`}},{id:"beh2",category:"geometry",labelKey:"opt.scenBeH2",descKey:"opt.descBeH2",paramType:"angle",defaultMin:90,defaultMax:270,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>{const o=t*Math.PI/180/2,i=1.33*Math.sin(o),e=1.33*Math.cos(o);return`3
BeH2 angle=${t.toFixed(1)}
Be  0.0  0.0  0.0
H   ${i.toFixed(6)}  0.0  ${e.toFixed(6)}
H   ${(-i).toFixed(6)}  0.0  ${e.toFixed(6)}`}},{id:"nh3",category:"geometry",labelKey:"opt.scenNH3",descKey:"opt.descNH3",paramType:"height",defaultMin:0,defaultMax:.5,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>{const r=Math.sqrt(Math.max(0,1.024144-t*t)),o=Math.sqrt(3)/2;return["4",`NH3 h=${t.toFixed(3)}`,`N   0.000000  0.000000  ${t.toFixed(6)}`,`H   ${r.toFixed(6)}  0.000000  0.000000`,`H   ${(-r/2).toFixed(6)}  ${(r*o).toFixed(6)}  0.000000`,`H   ${(-r/2).toFixed(6)}  ${(-r*o).toFixed(6)}  0.000000`].join(`
`)}},{id:"ch2",category:"geometry",labelKey:"opt.scenCH2",descKey:"opt.descCH2",paramType:"angle",defaultMin:90,defaultMax:180,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>{const o=t*Math.PI/180/2,i=1.08*Math.sin(o),e=1.08*Math.cos(o);return`3
CH2 angle=${t.toFixed(1)}
C  0.0  0.0  0.0
H  ${i.toFixed(6)}  0.0  ${e.toFixed(6)}
H  ${(-i).toFixed(6)}  0.0  ${e.toFixed(6)}`}},{id:"fhf",category:"geometry",labelKey:"opt.scenFHF",descKey:"opt.descFHF",paramType:"bond",defaultMin:.8,defaultMax:1.5,defaultSteps:20,defaultCharge:-1,defaultMult:1,basisOptions:["STO-3G"],generateXYZ:t=>["3",`FHF- FH=${t.toFixed(3)}`,"F   0.000000  0.000000  0.000000",`H   0.000000  0.000000  ${t.toFixed(6)}`,`F   0.000000  0.000000  ${2.3.toFixed(6)}`].join(`
`)},{id:"h3plus",category:"geometry",labelKey:"opt.scenH3plus",descKey:"opt.descH3plus",paramType:"bond",defaultMin:.5,defaultMax:2.5,defaultSteps:20,defaultCharge:1,defaultMult:1,generateXYZ:t=>{const a=Math.sqrt(3)/2;return["3",`H3+ R=${t.toFixed(3)}`,"H   0.000000  0.000000  0.000000",`H   ${t.toFixed(6)}  0.000000  0.000000`,`H   ${(t/2).toFixed(6)}  ${(t*a).toFixed(6)}  0.000000`].join(`
`)}},{id:"c2h6",category:"geometry",labelKey:"opt.scenC2H6",descKey:"opt.descC2H6",paramType:"dihedral",defaultMin:0,defaultMax:120,defaultSteps:12,defaultCharge:0,defaultMult:1,basisOptions:["STO-3G"],generateXYZ:t=>{const o=109.47*Math.PI/180,i=1.09*Math.cos(Math.PI-o),e=1.09*Math.sin(Math.PI-o),c=["8",`C2H6 dih=${t.toFixed(1)}`];c.push("C   0.000000  0.000000  0.000000"),c.push(`C   0.000000  0.000000  ${1.54.toFixed(6)}`);const l=t*Math.PI/180;for(let n=0;n<3;n++){const p=n*2*Math.PI/3;c.push(`H   ${(e*Math.cos(p)).toFixed(6)}  ${(e*Math.sin(p)).toFixed(6)}  ${(-i).toFixed(6)}`)}for(let n=0;n<3;n++){const p=n*2*Math.PI/3+l;c.push(`H   ${(e*Math.cos(p)).toFixed(6)}  ${(e*Math.sin(p)).toFixed(6)}  ${(1.54+i).toFixed(6)}`)}return c.join(`
`)}}],Qt=["STO-3G","3-21G","6-31G"];let H=Lt[0],et="STO-3G",R=!0,P=!1,D=0,ot=1,I="HF",K=H.defaultMin,Y=H.defaultMax,W=H.defaultSteps,X=!1,at=!1,C=[],G=!1,At=0,B=0,U=null,V=null;function te(t,a){const r=Ft[t.z]??1,o=Ft[a.z]??1,i=r+o,e=t.z+a.z,c=t.symbol===a.symbol?`${t.symbol}₂`:`${t.symbol}${a.symbol}`;return{id:`custom_${t.symbol}_${a.symbol}`,category:"dissociation",labelKey:"",descKey:"",paramType:"bond",defaultMin:Math.max(.4,Math.round(i*.5*10)/10),defaultMax:Math.min(5,Math.round(i*3*10)/10),defaultSteps:20,defaultCharge:0,defaultMult:e%2===0?1:2,basisOptions:["STO-3G"],generateXYZ:l=>`2
${c} R=${l.toFixed(3)}
${t.symbol}  0.0  0.0  0.0
${a.symbol}  0.0  0.0  ${l.toFixed(6)}`}}const Ht=new Map;async function ee(t){const a=Ht.get(t);if(a)return a;const r=`/GANSU-Lite/basis/${t.toLowerCase()}.gbs`,o=await fetch(r);if(!o.ok)throw new Error(`Failed to load basis set: ${t}`);const i=await o.text(),e=qt.fromGBS(i);return Ht.set(t,e),e}function J(t){return t==="rhf"?ct()?"#00d4ff":"#0077cc":ct()?"#ff8844":"#cc4400"}function pt(){const t=[];return R&&t.push({key:"rhf",label:_(I,!0),color:J("rhf")}),P&&t.push({key:"uhf",label:_(I,!1),color:J("uhf")}),t}const E=document.getElementById("app");function q(){const t=H.paramType==="angle"||H.paramType==="dihedral"?"°":"Å";E.innerHTML=`
    <div class="opt-page">
      ${Ut("optimize")}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${F("opt.scenario")}</h2>
          <div class="opt-scenario-grid" id="scenario-grid"></div>

          <h2>${F("opt.basis")}</h2>
          <div class="opt-basis-row" id="basis-row"></div>

          <h2>${F("opt.method")}</h2>
          <div class="opt-method-row">
            <button id="btn-rhf" class="opt-method-btn${R?" active":""}" style="--method-color:${J("rhf")}">${_(I,!0)}</button>
            ${H.category!=="geometry"?`<button id="btn-uhf" class="opt-method-btn${P?" active":""}" style="--method-color:${J("uhf")}">${_(I,!1)}</button>`:""}
          </div>

          <div class="opt-charge-row">
            <label>
              <span>${F("opt.charge")}</span>
              <input id="inp-charge" type="number" min="-3" max="3" value="${D}" />
            </label>
            <label>
              <span>${F("opt.mult")}</span>
              <input id="inp-mult" type="number" min="1" max="5" value="${ot}" />
            </label>
          </div>

          <h2>${F("opt.param")}</h2>
          <div class="opt-param-form">
            <label>
              <span>${F("opt.min")}</span>
              <input id="param-min" type="number" step="0.1" value="${K}" />
              <span class="opt-unit">${t}</span>
            </label>
            <label>
              <span>${F("opt.max")}</span>
              <input id="param-max" type="number" step="0.1" value="${Y}" />
              <span class="opt-unit">${t}</span>
            </label>
            <label>
              <span>${F("opt.steps")}</span>
              <input id="param-steps" type="number" min="3" max="50" value="${W}" />
            </label>
          </div>

          <div class="theory-row" style="margin:8px 0;display:flex;align-items:center;gap:8px;font-size:0.9rem;">
            <span>Theory:</span>${Zt("theory-sel",I,"",Gt)}
          </div>

          <button id="run-btn" class="opt-run-btn" ${X?"disabled":""}>
            ${X?F("opt.running"):F("opt.run")}
          </button>
          ${X?`<button id="stop-btn" class="opt-stop-btn">${F("opt.stop")}</button>`:""}

          <div id="progress-area"></div>
          <div id="result-summary"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="mol-vis"></div>
          <div id="graph-container">
            ${!G&&!X?`<p class="opt-hint">${F("opt.waiting")}</p>`:""}
          </div>
          <div id="scan-slider-area"></div>
        </div>
      </div>
    </div>`,ce();const a=E.querySelector("#scenario-grid"),r=[{key:"dissociation",labelKey:"opt.catDissociation"},{key:"geometry",labelKey:"opt.catGeometry"}];for(const l of r){const n=Lt.filter(u=>u.category===l.key);if(n.length===0)continue;const p=document.createElement("div");p.className="opt-category-header",p.textContent=F(l.labelKey),a.appendChild(p);const x=document.createElement("div");x.className="opt-category-row";for(const u of n){const b=document.createElement("div");b.className="opt-scenario-card"+(u.id===H.id?" selected":""),b.innerHTML=`<strong>${F(u.labelKey)}</strong><span>${F(u.descKey)}</span>`,b.addEventListener("click",()=>{X||(H=u,K=u.defaultMin,Y=u.defaultMax,W=u.defaultSteps,D=u.defaultCharge,ot=u.defaultMult,u.category==="geometry"&&(P=!1,R=!0),U=null,V=null,C=[],G=!1,q())}),x.appendChild(b)}a.appendChild(x)}{const l=document.createElement("div");l.className="opt-category-header",l.textContent=F("opt.catCustom"),a.appendChild(l);const n=document.createElement("div");n.className="opt-category-row";const p=H.id.startsWith("custom_"),x=document.createElement("div");x.className="opt-scenario-card"+(p?" selected":"");const u=U&&V?`${U.symbol} + ${V.symbol}`:F("opt.customSelect");x.innerHTML=`<strong>${u}</strong><span>${F("opt.catCustom")}</span>`,x.addEventListener("click",()=>{X||Jt((b,h,M,$)=>{U={symbol:b,z:h},V={symbol:M,z:$},o()},U,V)}),n.appendChild(x),a.appendChild(n)}zt(a,{label:F("opt.scenario")});function o(){if(!U||!V)return;const l=te(U,V);H=l,K=l.defaultMin,Y=l.defaultMax,W=l.defaultSteps,D=l.defaultCharge,ot=l.defaultMult,(U.z+V.z-D)%2!==0?(R=!1,P=!0):(R=!0,P=!0),C=[],G=!1,q()}const i=E.querySelector("#basis-row"),e=H.basisOptions??Qt;e.includes(et)||(et=e[0]);for(const l of e){const n=document.createElement("button");n.className="opt-basis-btn"+(l===et?" selected":""),n.textContent=l,n.addEventListener("click",()=>{X||(et=l,C=[],G=!1,q())}),i.appendChild(n)}E.querySelector("#btn-rhf").addEventListener("click",()=>{X||R&&!P||(R=!R,C=[],G=!1,q())}),E.querySelector("#btn-uhf")?.addEventListener("click",()=>{X||P&&!R||(P=!P,C=[],G=!1,q())}),E.querySelector("#inp-charge").addEventListener("change",l=>{D=parseInt(l.target.value,10)}),E.querySelector("#inp-mult").addEventListener("change",l=>{ot=parseInt(l.target.value,10)}),E.querySelector("#nav-theme").addEventListener("click",()=>{Vt(),q()}),E.querySelector("#nav-lang").addEventListener("click",()=>{_t(),q()});const c=E.querySelector("#theory-sel");c&&c.addEventListener("change",()=>{I=c.value,q()}),E.querySelector("#run-btn").addEventListener("click",ae),E.querySelector("#stop-btn")?.addEventListener("click",()=>{at=!0}),E.querySelector("#param-min").addEventListener("change",l=>{K=parseFloat(l.target.value)}),E.querySelector("#param-max").addEventListener("change",l=>{Y=parseFloat(l.target.value)}),E.querySelector("#param-steps").addEventListener("change",l=>{const n=l.target;let p=parseInt(n.value,10);(isNaN(p)||p<3)&&(p=3),p>50&&(p=50),n.value=String(p),W=p}),C.length>0&&rt(C,G),G?(ne(),se(),Kt()):st([{param:(K+Y)/2}])}function oe(t,a,r,o){if(o<=0||o>=r)return;const i=Math.PI/4,e=Math.cos(i),c=Math.sin(i),l=o-1,n=o,p=a.clone();for(let h=0;h<r;h++)p.set(h,l,e*a.get(h,l)+c*a.get(h,n));const x=a.clone();for(let h=0;h<r;h++)x.set(h,l,e*a.get(h,l)-c*a.get(h,n));const u=new Mt(r,r),b=new Mt(r,r);for(let h=0;h<r;h++)for(let M=0;M<r;M++){let $=0,v=0;for(let S=0;S<o;S++)$+=p.get(h,S)*p.get(M,S),v+=x.get(h,S)*x.get(M,S);u.set(h,M,$),b.set(h,M,v)}t.setInitialDensityGuessAlphaBeta(u,b)}async function ae(){if(X)return;if(I!=="HF"&&P&&R&&W>=15){const n=Math.round(W*2*8/60);if(!confirm(`This scan combines DFT × RHF+UHF × ${W} points → estimated ~${n} min.
Proceed?`))return}X=!0,at=!1,C=[],G=!1,q();const t=Math.max(3,Math.min(50,W)),a=[];for(let n=0;n<=t;n++)a.push(K+(Y-K)*n/t);const r=Math.floor((ot-1)/2),o=performance.now(),i=(R?a.length:0)+(P?a.length:0);let e=0;try{const n=await ee(et);let p=null;if(R)for(let x=0;x<a.length&&!at;x++){const u=a[x],b=nt(H.generateXYZ(u)),h=new it(b,n,D,r),M=await lt(h,n,I,"RHF");p&&M.setInitialDensityGuess(p);const $=await M.solve({eriBackend:"js"});p=M.density.clone(),C.push({param:u,rhf:$}),e++,kt(e,i),rt(C,!1);{const v=[{param:u,label:_(I,!0),color:J("rhf")}];P&&v.push({param:(K+Y)/2,label:_(I,!1),color:J("uhf")}),st(v)}await new Promise(v=>setTimeout(v,0))}else for(const x of a)C.push({param:x,rhf:NaN});if(P&&!at){let x=null,u=null,b=null;if((H.defaultCharge!==void 0?(a.length>0?nt(H.generateXYZ(a[0])).reduce(($,v)=>$+v.atomicNumber,0):0)-D:0)%2===0){const $=nt(H.generateXYZ(a[a.length-1])),v=new it($,n,D,r),S=await lt(v,n,I,"RHF");p&&S.setInitialDensityGuess(p),await S.solve({eriBackend:"js"}),b=S.coefficients}for(let $=a.length-1;$>=0&&!at;$--){const v=a[$],S=nt(H.generateXYZ(v));try{const g=new it(S,n,D,r),m=await lt(g,n,I,"UHF");x&&u?m.setInitialDensityGuessAlphaBeta(x,u):b&&oe(m,b,g.numBasis,g.numAlphaSpins);let A=1/0;const L=await m.solve({eriBackend:"js",onIteration:(O,T,Z)=>{A=Math.abs(Z)}});A<1e-6?(x=m.densityAlphaMatrix.clone(),u=m.densityBetaMatrix.clone(),C[$].uhf=L):console.warn(`UHF not converged at param=${v}, deltaE=${A}`)}catch(g){console.warn("UHF error at param",v,g)}e++,kt(e,i),rt(C,!1);{const g=[];if(R){const m=C.filter(L=>isFinite(L.rhf)),A=m.length>0?m.reduce((L,O)=>L.rhf<O.rhf?L:O):null;g.push({param:A?.param??(K+Y)/2,label:_(I,!0),color:J("rhf")})}g.push({param:v,label:_(I,!1),color:J("uhf")}),st(g)}await new Promise(g=>setTimeout(g,0))}}}catch(n){console.error("Scan error:",n)}At=performance.now()-o,X=!1,G=!0;const c=R?"rhf":"uhf";let l=0;for(let n=1;n<C.length;n++){const p=C[n][c],x=C[l][c];p!=null&&isFinite(p)&&(x==null||!isFinite(x)||p<x)&&(l=n)}B=l,q()}function kt(t,a){const r=E.querySelector("#progress-area");if(!r)return;const o=(t/a*100).toFixed(0),i=F("opt.progress").replace("{n}",String(t)).replace("{total}",String(a));r.innerHTML=`
    <div class="opt-progress">
      <div class="opt-progress-bar" style="width:${o}%"></div>
    </div>
    <p class="opt-progress-text">${i}</p>`}function ne(){const t=E.querySelector("#result-summary");if(!t||C.length===0)return;const a=H.paramType==="angle"||H.paramType==="dihedral"?"°":"Å",r=(At/1e3).toFixed(1),o=pt();let i="";for(const e of o){const c=C.filter(n=>n[e.key]!=null&&isFinite(n[e.key])).map(n=>({param:n.param,e:n[e.key]}));if(c.length===0)continue;let l=c[0];for(const n of c)n.e<l.e&&(l=n);i+=`<tr>
      <td><span class="opt-dot" style="background:${e.color}"></span>${e.label}</td>
      <td>${l.param.toFixed(4)} ${a}</td>
      <td><strong>${l.e.toFixed(8)} Eh</strong></td>
    </tr>`}t.innerHTML=`
    <div class="opt-summary">
      <h3>${F("opt.done")} (${r} s)</h3>
      <table>
        <tr><th></th><th>${F("opt.resultParam")}</th><th>${F("opt.resultEnergy")}</th></tr>
        ${i}
      </table>
    </div>`}function se(){const t=E.querySelector("#scan-slider-area");if(!t||C.length===0)return;B>=C.length&&(B=0);const a=H.paramType==="angle"||H.paramType==="dihedral",r=a?"°":" Å",o=a?C[B].param.toFixed(1):C[B].param.toFixed(3);t.innerHTML=`
    <div class="walsh-slider">
      <input type="range" id="scan-slider" min="0" max="${C.length-1}" value="${B}" />
      <div class="walsh-slider-label">${o}${r}</div>
    </div>`,E.querySelector("#scan-slider")?.addEventListener("input",i=>{B=parseInt(i.target.value);const e=C[B].param,c=E.querySelector(".walsh-slider-label");c&&(c.textContent=`${a?e.toFixed(1):e.toFixed(3)}${r}`),rt(C,!0),Kt()})}function Kt(){if(C.length===0)return;B>=C.length&&(B=0);const t=C[B],a=pt(),r=[];for(const o of a)r.push({param:t.param,label:o.label,color:o.color});r.length>0&&st(r)}const re={H:"#999",He:"#0CC",Li:"#C2C",Be:"#6C0",B:"#F90",C:"#555",N:"#35F",O:"#F22",F:"#9E5"},Ct={H:.31,He:.28,Li:1.28,Be:.96,B:.84,C:.76,N:.71,O:.66,F:.57};function Nt(t){if(H.paramType==="height")return{pU:"x",pV:"z",pW:"y"};const a=["x","y","z"].map(r=>({a:r,r:Math.max(...t.map(o=>o[r]))-Math.min(...t.map(o=>o[r]))})).sort((r,o)=>o.r-r.r);return{pU:a[0].a,pV:a[1].a,pW:a[2].a}}function Rt(t,a){const r=a??Nt(t);if(H.paramType==="height"){const o=Math.cos(Math.PI/6),i=Math.sin(Math.PI/6);return t.map(e=>({u:e.x*o+e.y*i,v:e.z}))}return t.map(o=>({u:o[r.pU],v:o[r.pV]}))}function dt(t){const a=t.split(`
`),r=parseInt(a[0].trim()),o=[];for(let i=2;i<2+r;i++){const e=a[i].trim().split(/\s+/);o.push({sym:e[0],x:+e[1],y:+e[2],z:+e[3]})}return o}function ie(t,a,r){let o=[];for(const g of[K,Y])o.push(...dt(H.generateXYZ(g)));const i=Nt(o);let e=[],c=[];for(const g of[K,Y]){const m=Rt(dt(H.generateXYZ(g)),i);e.push(...m.map(A=>A.u)),c.push(...m.map(A=>A.v))}const l=Math.min(...e),n=Math.max(...e),p=Math.min(...c),x=Math.max(...c),u=n-l||.01,b=x-p||.01,h=r,M=r+30,$=Math.min((t-2*r)/u,(a-h-M)/b,120),v=(p+x)/2,S=(h-M)/2;return{scale:$,uMid:(l+n)/2,vMid:v,cySvgOffset:S,axes:i}}function St(t,a,r,o,i){const e=dt(H.generateXYZ(t));if(e.length===0)return"";const c=Et(),l=35,{scale:n,uMid:p,vMid:x,cySvgOffset:u,axes:b}=ie(a,r,l),h=a/2,M=r/2+u,$=s=>h+(s-p)*n,v=s=>M-(s-x)*n,S=18,g=Rt(e,b);let m="";for(let s=0;s<e.length;s++)for(let d=s+1;d<e.length;d++){const f=e[s].x-e[d].x,y=e[s].y-e[d].y,k=e[s].z-e[d].z,w=Math.sqrt(f*f+y*y+k*k),N=1.4*((Ct[e[s].sym]??.7)+(Ct[e[d].sym]??.7));w<N&&(m+=`<line x1="${$(g[s].u)}" y1="${v(g[s].v)}" x2="${$(g[d].u)}" y2="${v(g[d].v)}" stroke="${c.grid}" stroke-width="4" stroke-linecap="round"/>`)}const A=e.map((s,d)=>d).sort((s,d)=>e[s][b.pW]-e[d][b.pW]);for(const s of A){const d=$(g[s].u),f=v(g[s].v),y=re[e[s].sym]??"#888";m+=`<circle cx="${d}" cy="${f}" r="${S}" fill="${y}" stroke="${c.axis}" stroke-width="1.2"/>`,m+=`<text x="${d}" y="${f}" text-anchor="middle" dy="0.38em" font-size="11" font-weight="bold" fill="#fff" stroke="#0003" stroke-width="0.3">${e[s].sym}</text>`}const L=H.paramType==="angle"||H.paramType==="dihedral",O=L?"°":"Å",Z=`${L?t.toFixed(1):t.toFixed(3)} ${O}`,z=i??c.accent;if(H.paramType==="bond"){const[s,d]=H.paramAtomPair??[0,1],f=$(g[s].u),y=$(g[d].u),k=Math.max(v(g[s].v),v(g[d].v))+S+10;m+=`<line x1="${f}" y1="${k}" x2="${y}" y2="${k}" stroke="${z}" stroke-width="1.2"/>`,m+=`<line x1="${f}" y1="${k-4}" x2="${f}" y2="${k+4}" stroke="${z}" stroke-width="1.2"/>`,m+=`<line x1="${y}" y1="${k-4}" x2="${y}" y2="${k+4}" stroke="${z}" stroke-width="1.2"/>`;const w=y>f?1:-1;m+=`<polygon points="${f},${k} ${f+w*6},${k-3} ${f+w*6},${k+3}" fill="${z}"/>`,m+=`<polygon points="${y},${k} ${y-w*6},${k-3} ${y-w*6},${k+3}" fill="${z}"/>`,m+=`<text x="${(f+y)/2}" y="${k+14}" text-anchor="middle" font-size="11" font-weight="600" fill="${z}">${Z}</text>`}else if(H.paramType==="angle"){const s=$(g[0].u),d=v(g[0].v),f=$(g[1].u),y=v(g[1].v),k=$(g[2].u),w=v(g[2].v),N=Math.atan2(y-d,f-s),Q=Math.atan2(w-d,k-s),j=22,ut=s+j*Math.cos(N),ft=d+j*Math.sin(N),ht=s+j*Math.cos(Q),mt=d+j*Math.sin(Q),It=(ut-s)*(mt-d)-(ft-d)*(ht-s)>0?1:0;m+=`<path d="M ${ut} ${ft} A ${j} ${j} 0 0 ${It} ${ht} ${mt}" fill="none" stroke="${z}" stroke-width="1.5"/>`;let tt=Q-N;tt>Math.PI&&(tt-=2*Math.PI),tt<-Math.PI&&(tt+=2*Math.PI);const xt=N+tt/2,gt=j+14;m+=`<text x="${s+gt*Math.cos(xt)}" y="${d+gt*Math.sin(xt)}" text-anchor="middle" dy="0.35em" font-size="11" font-weight="600" fill="${z}">${Z}</text>`;const yt=e[1].x-e[0].x,$t=e[1].y-e[0].y,bt=e[1].z-e[0].z,Bt=Math.sqrt(yt*yt+$t*$t+bt*bt),Pt=(s+f)/2,Ot=(d+y)/2,vt=Math.atan2(y-d,f-s),Xt=Pt+12*Math.cos(vt+Math.PI/2),Yt=Ot+12*Math.sin(vt+Math.PI/2);m+=`<text x="${Xt}" y="${Yt}" text-anchor="middle" dy="0.35em" font-size="9" fill="${c.dim}">${Bt.toFixed(2)} Å</text>`}else if(H.paramType==="height"){const s=$(g[0].u),d=v(g[0].v),f=v(0),y=Math.min(...g.slice(1).map(w=>$(w.u)))-20,k=Math.max(...g.slice(1).map(w=>$(w.u)))+20;if(m+=`<line x1="${y}" y1="${f}" x2="${k}" y2="${f}" stroke="${c.dim}" stroke-width="1" stroke-dasharray="5,3"/>`,m+=`<text x="${k+4}" y="${f}" dy="0.35em" font-size="9" fill="${c.dim}">H₃ plane</text>`,Math.abs(d-f)>8){const w=s>h?-24:24;m+=`<line x1="${s+w}" y1="${d}" x2="${s+w}" y2="${f}" stroke="${z}" stroke-width="1.5"/>`,m+=`<line x1="${s+w-5}" y1="${d}" x2="${s+w+5}" y2="${d}" stroke="${z}" stroke-width="1.5"/>`,m+=`<line x1="${s+w-5}" y1="${f}" x2="${s+w+5}" y2="${f}" stroke="${z}" stroke-width="1.5"/>`;const N=d<f?1:-1;m+=`<polygon points="${s+w},${d} ${s+w-3},${d+N*6} ${s+w+3},${d+N*6}" fill="${z}"/>`,m+=`<polygon points="${s+w},${f} ${s+w-3},${f-N*6} ${s+w+3},${f-N*6}" fill="${z}"/>`;const Q=s+w+(w>0?10:-10),j=w>0?"start":"end";m+=`<text x="${Q}" y="${(d+f)/2}" text-anchor="${j}" dy="0.35em" font-size="11" font-weight="600" fill="${z}">h = ${Z}</text>`}else m+=`<text x="${s}" y="${f+S+14}" text-anchor="middle" font-size="11" font-weight="600" fill="${z}">h = ${Z}</text>`}else if(H.paramType==="dihedral"){const s=Math.max(...g.map(d=>v(d.v)))+S+12;m+=`<text x="${h}" y="${s}" text-anchor="middle" font-size="11" font-weight="600" fill="${z}">φ = ${Z}</text>`}return o&&(m+=`<text x="${a/2}" y="14" text-anchor="middle" font-size="11" font-weight="700" fill="${z}">${o}</text>`),m}function st(t){const a=E.querySelector("#mol-vis");if(!a||t.length===0)return;const r=520,o=200;if(t.length===1){const i=t[0],e=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${r} ${o}" style="width:100%;max-width:${r}px;display:block;margin:0 auto;">${St(i.param,r,o,i.label,i.color)}</svg>`;a.innerHTML=e}else{const i=Math.floor(r/t.length);let e='<div style="display:flex;gap:4px;justify-content:center;">';for(const c of t)e+=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${i} ${o}" style="flex:1;max-width:${i}px;">${St(c.param,i,o,c.label,c.color)}</svg>`;e+="</div>",a.innerHTML=e}}function rt(t,a){const r=E.querySelector("#graph-container");if(!r||t.length<1)return;const o=Et(),e={bond:F("opt.xBond"),angle:F("opt.xAngle"),height:F("opt.xHeight"),dihedral:F("opt.xDihedral")}[H.paramType],c=pt(),l=520,n=360,p=72,x=24,u=36,b=44,h=l-p-x,M=n-u-b,$=[];for(const s of c)for(const d of t){const f=d[s.key];f!=null&&isFinite(f)&&$.push(f)}$.length===0&&$.push(0);const v=Math.min(...$),S=Math.max(...$),g=(S-v)*.1||.01,m=v-g,A=S+g,L=s=>p+(s-K)/(Y-K||1)*h,O=s=>u+M-(s-m)/(A-m||1)*M;let T=`<svg width="${l}" height="${n}" xmlns="http://www.w3.org/2000/svg" style="display:block;width:100%;height:auto;max-width:${l}px;" viewBox="0 0 ${l} ${n}">`;T+=`<rect x="${p}" y="${u}" width="${h}" height="${M}" fill="${o.surface}" rx="2"/>`;for(let s=0;s<=5;s++){const d=m+(A-m)*s/5,f=O(d);T+=`<line x1="${p}" y1="${f}" x2="${p+h}" y2="${f}" stroke="${o.grid}" stroke-width="0.5"/>`,T+=`<text x="${p-6}" y="${f+3}" text-anchor="end" font-size="9" font-family="monospace" fill="${o.dim}">${d.toFixed(4)}</text>`}for(let s=0;s<=5;s++){const d=K+(Y-K)*s/5,f=L(d);T+=`<line x1="${f}" y1="${u}" x2="${f}" y2="${u+M}" stroke="${o.grid}" stroke-width="0.5"/>`,T+=`<text x="${f}" y="${u+M+14}" text-anchor="middle" font-size="9" fill="${o.dim}">${d.toFixed(2)}</text>`}T+=`<line x1="${p}" y1="${u}" x2="${p}" y2="${u+M}" stroke="${o.axis}" stroke-width="1"/>`,T+=`<line x1="${p}" y1="${u+M}" x2="${p+h}" y2="${u+M}" stroke="${o.axis}" stroke-width="1"/>`;for(const s of c){const d=t.filter(y=>y[s.key]!=null&&isFinite(y[s.key])).map(y=>({x:y.param,y:y[s.key]}));if(d.length<2){for(const y of d)T+=`<circle cx="${L(y.x).toFixed(1)}" cy="${O(y.y).toFixed(1)}" r="3" fill="${s.color}"/>`;continue}let f="";for(let y=0;y<d.length;y++){const k=L(d[y].x),w=O(d[y].y);f+=y===0?`M${k.toFixed(1)},${w.toFixed(1)}`:` L${k.toFixed(1)},${w.toFixed(1)}`}T+=`<path d="${f}" fill="none" stroke="${s.color}" stroke-width="2" stroke-linejoin="round"/>`;for(const y of d)T+=`<circle cx="${L(y.x).toFixed(1)}" cy="${O(y.y).toFixed(1)}" r="3" fill="${s.color}" stroke="${o.surface}" stroke-width="0.8"/>`;if(d.length>=3&&a){let y=d[0];for(const N of d)N.y<y.y&&(y=N);const k=L(y.x),w=O(y.y);T+=`<line x1="${k}" y1="${w}" x2="${k}" y2="${u+M}" stroke="${s.color}" stroke-width="1" stroke-dasharray="3,3" opacity="0.5"/>`,T+=le(k,w,6,s.color)}}if(a&&B>=0&&B<t.length){const s=t[B],d=L(s.param),f=ct()?"#ffd700":"#cc8800";T+=`<line x1="${d}" y1="${u}" x2="${d}" y2="${u+M}" stroke="${f}" stroke-width="1" stroke-dasharray="3,2" opacity="0.6"/>`;for(const y of c){const k=s[y.key];k!=null&&isFinite(k)&&(T+=`<circle cx="${d}" cy="${O(k)}" r="5" fill="${y.color}" stroke="${f}" stroke-width="2"/>`)}}const Z=p+8;let z=u+14;for(const s of c)T+=`<rect x="${Z}" y="${z-7}" width="10" height="3" rx="1" fill="${s.color}"/>`,T+=`<text x="${Z+14}" y="${z-3}" font-size="9" font-weight="600" fill="${o.dim}">${s.label}</text>`,z+=14;T+=`<text x="${p+h/2}" y="20" text-anchor="middle" font-size="12" font-weight="600" fill="${o.titleSvg}">${F("opt.graphTitle")}</text>`,T+=`<text x="${p+h/2}" y="${n-4}" text-anchor="middle" font-size="10" fill="${o.dim}">${e}</text>`,T+=`<text x="14" y="${u+M/2}" text-anchor="middle" font-size="10" fill="${o.dim}" transform="rotate(-90,14,${u+M/2})">${F("opt.yEnergy")}</text>`,T+="</svg>",r.innerHTML=T}function le(t,a,r,o){const i=[];for(let e=0;e<10;e++){const c=Math.PI/2+e*Math.PI/5,l=e%2===0?r:r*.4;i.push(`${(t+l*Math.cos(c)).toFixed(1)},${(a-l*Math.sin(c)).toFixed(1)}`)}return`<polygon points="${i.join(" ")}" fill="${o}" stroke="${o}" stroke-width="0.5"/>`}let Tt=!1;function ce(){if(Tt)return;Tt=!0;const t=document.createElement("style");t.textContent=`
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
    .opt-category-header {
      font-size: 0.7rem; font-weight: 600; color: var(--color-text-dim);
      text-transform: uppercase; letter-spacing: 0.03em;
      margin-top: 4px;
    }
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
    .opt-scenario-card span { display: none; }

    .opt-basis-row { display: flex; gap: 6px; }
    .opt-basis-btn {
      flex: 1; padding: 6px 0; border: 1px solid var(--color-border); border-radius: 6px;
      background: none; cursor: pointer; font-size: 0.78rem; color: var(--color-text);
      transition: all 0.15s;
    }
    .opt-basis-btn:hover { background: var(--color-surface-alt); }
    .opt-basis-btn.selected {
      border-color: var(--color-accent); color: var(--color-accent);
      font-weight: 600; box-shadow: 0 0 0 1px var(--color-accent);
    }

    .opt-method-row { display: flex; gap: 6px; }
    .opt-method-btn {
      flex: 1; padding: 7px 0; border: 2px solid var(--method-color); border-radius: 6px;
      background: none; cursor: pointer; font-size: 0.82rem; font-weight: 600;
      color: var(--method-color); transition: all 0.15s;
    }
    .opt-method-btn:hover { background: color-mix(in srgb, var(--method-color) 12%, transparent); }
    .opt-method-btn.active {
      background: var(--method-color); color: #fff;
    }

    .opt-charge-row {
      display: flex; gap: 10px; margin-top: 10px;
    }
    .opt-charge-row label {
      flex: 1; display: flex; align-items: center; gap: 6px; font-size: 0.72rem;
      color: var(--color-text-secondary);
    }
    .opt-charge-row input {
      width: 50px; padding: 4px 6px; border: 1px solid var(--color-border-input); border-radius: 5px;
      background: var(--color-input); color: var(--color-text); font-size: 0.8rem;
      outline: none; font-family: monospace; text-align: center;
    }

    .opt-param-form { display: flex; flex-direction: column; gap: 6px; }
    .opt-param-form label {
      display: flex; align-items: center; gap: 6px; font-size: 0.78rem;
    }
    .opt-param-form label span:first-child {
      flex: 0 0 50px; color: var(--color-text-secondary); font-size: 0.72rem;
    }
    .opt-param-form input {
      flex: 1; padding: 5px 8px; border: 1px solid var(--color-border-input); border-radius: 5px;
      background: var(--color-input); color: var(--color-text); font-size: 0.8rem;
      outline: none; font-family: monospace;
    }
    .opt-unit { font-size: 0.72rem; color: var(--color-text-dim); flex: 0 0 16px; }

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
    .opt-summary strong { color: var(--color-text); }
    .opt-dot {
      display: inline-block; width: 8px; height: 8px; border-radius: 50%;
      margin-right: 5px; vertical-align: middle;
    }

    #mol-vis { width: 100%; margin-bottom: 8px; }
    #graph-container { width: 100%; text-align: center; }

    .walsh-slider {
      width: 100%; max-width: 520px; margin: 8px auto 0;
      display: flex; align-items: center; gap: 10px;
    }
    .walsh-slider input[type=range] { flex: 1; cursor: pointer; }
    .walsh-slider-label {
      font-size: 0.85rem; font-weight: 600; min-width: 56px; text-align: center;
      color: var(--color-text);
    }
    .opt-hint { color: var(--color-text-dim); font-size: 0.85rem; padding: 60px 20px; }

    @media (max-width: 700px) {
      .opt-content { flex-direction: column; }
      .opt-controls { flex: none; }
    }
  `,document.head.appendChild(t)}jt();Dt();q();
