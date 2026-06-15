import"./styles-B8WF9G-F.js";import{p as nt}from"./parseXYZ-ByfvMHmk.js";import{g as it,f as Yt,M as Mt}from"./ri-DT3jfd0v.js";import{t as _,a as Zt,b as lt,H as qt}from"./theoryControls-CuXHSm5X.js";import{t as H,i as Gt,c as jt,r as Ut,a as Vt,b as Dt,d as ct,g as Et}from"./nav-oJm8MoOS.js";const Ft={1:.31,2:.28,3:1.28,4:.96,5:.84,6:.76,7:.71,8:.66,9:.57,10:.58,11:1.66,12:1.41,13:1.21,14:1.11,15:1.07,16:1.05,17:1.02,18:1.06,19:2.03,20:1.76,21:1.7,22:1.6,23:1.53,24:1.39,25:1.39,26:1.32,27:1.26,28:1.24,29:1.32,30:1.22,31:1.22,32:1.2,33:1.19,34:1.2,35:1.2,36:1.16},_t=[[0,0,1,"H"],[0,17,2,"He"],[1,0,3,"Li"],[1,1,4,"Be"],[1,12,5,"B"],[1,13,6,"C"],[1,14,7,"N"],[1,15,8,"O"],[1,16,9,"F"],[1,17,10,"Ne"],[2,0,11,"Na"],[2,1,12,"Mg"],[2,12,13,"Al"],[2,13,14,"Si"],[2,14,15,"P"],[2,15,16,"S"],[2,16,17,"Cl"],[2,17,18,"Ar"],[3,0,19,"K"],[3,1,20,"Ca"],[3,2,21,"Sc"],[3,3,22,"Ti"],[3,4,23,"V"],[3,5,24,"Cr"],[3,6,25,"Mn"],[3,7,26,"Fe"],[3,8,27,"Co"],[3,9,28,"Ni"],[3,10,29,"Cu"],[3,11,30,"Zn"],[3,12,31,"Ga"],[3,13,32,"Ge"],[3,14,33,"As"],[3,15,34,"Se"],[3,16,35,"Br"],[3,17,36,"Kr"]];function wt(t,a,r=18){const o=document.createElement("div");o.className="pte-grid";for(const[i,e,l,c]of _t){const n=document.createElement("div");n.className="pte-cell",n.style.gridRow=String(i+1),n.style.gridColumn=String(e+1);const p=l>r;p&&n.classList.add("pte-disabled");const $=l<=2||l<=4||l>=11&&l<=12||l>=19&&l<=20?"pte-s":l>=21&&l<=30?"pte-d":"pte-p";n.classList.add($),l===a&&n.classList.add("pte-selected"),n.innerHTML=`<span class="pte-z">${l}</span><span class="pte-sym">${c}</span>`,p||n.addEventListener("click",()=>t(c,l)),o.appendChild(n)}return o}function Wt(t,a,r){let o=a??null,i=r??null;const e=document.createElement("div");e.className="pte-overlay";const l=document.createElement("div");l.className="pte-modal pte-dual";function c(){l.innerHTML="";const n=document.createElement("div");n.className="pte-tables";const p=document.createElement("div");p.className="pte-col";const $=document.createElement("div");$.className="pte-title",$.textContent=`${H("opt.customAtomA")}${o?": "+o.symbol:""}`,p.appendChild($),p.appendChild(wt((b,S)=>{o={symbol:b,z:S},c()},o?.z,20)),n.appendChild(p);const u=document.createElement("div");u.className="pte-col";const k=document.createElement("div");k.className="pte-title",k.textContent=`${H("opt.customAtomB")}${i?": "+i.symbol:""}`,u.appendChild(k),u.appendChild(wt((b,S)=>{i={symbol:b,z:S},c()},i?.z,12)),n.appendChild(u),l.appendChild(n);const x=document.createElement("div");x.className="pte-btn-row";const v=document.createElement("button");v.className="pte-btn pte-btn-cancel",v.textContent="Cancel",v.addEventListener("click",()=>e.remove()),x.appendChild(v);const m=document.createElement("button");m.className="pte-btn pte-btn-ok",m.textContent="OK",m.disabled=!o||!i,m.addEventListener("click",()=>{o&&i&&(e.remove(),t(o.symbol,o.z,i.symbol,i.z))}),x.appendChild(m),l.appendChild(x)}c(),e.appendChild(l),e.addEventListener("click",n=>{n.target===e&&e.remove()}),document.body.appendChild(e)}const zt=[{id:"h2",category:"dissociation",labelKey:"opt.scenH2",descKey:"opt.descH2",paramType:"bond",defaultMin:.4,defaultMax:3,defaultSteps:20,defaultCharge:0,defaultMult:1,basisOptions:["3-21G","6-31G"],generateXYZ:t=>`2
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
`)}},{id:"c2h6",category:"geometry",labelKey:"opt.scenC2H6",descKey:"opt.descC2H6",paramType:"dihedral",defaultMin:0,defaultMax:120,defaultSteps:12,defaultCharge:0,defaultMult:1,basisOptions:["STO-3G"],generateXYZ:t=>{const o=109.47*Math.PI/180,i=1.09*Math.cos(Math.PI-o),e=1.09*Math.sin(Math.PI-o),l=["8",`C2H6 dih=${t.toFixed(1)}`];l.push("C   0.000000  0.000000  0.000000"),l.push(`C   0.000000  0.000000  ${1.54.toFixed(6)}`);const c=t*Math.PI/180;for(let n=0;n<3;n++){const p=n*2*Math.PI/3;l.push(`H   ${(e*Math.cos(p)).toFixed(6)}  ${(e*Math.sin(p)).toFixed(6)}  ${(-i).toFixed(6)}`)}for(let n=0;n<3;n++){const p=n*2*Math.PI/3+c;l.push(`H   ${(e*Math.cos(p)).toFixed(6)}  ${(e*Math.sin(p)).toFixed(6)}  ${(1.54+i).toFixed(6)}`)}return l.join(`
`)}}],Jt=["STO-3G","3-21G","6-31G"];let F=zt[0],et="STO-3G",K=!0,P=!1,U=0,ot=1,A="HF",I=F.defaultMin,Y=F.defaultMax,W=F.defaultSteps,X=!1,at=!1,C=[],G=!1,Lt=0,B=0,V=null,D=null;function Qt(t,a){const r=Ft[t.z]??1,o=Ft[a.z]??1,i=r+o,e=t.z+a.z,l=t.symbol===a.symbol?`${t.symbol}₂`:`${t.symbol}${a.symbol}`;return{id:`custom_${t.symbol}_${a.symbol}`,category:"dissociation",labelKey:"",descKey:"",paramType:"bond",defaultMin:Math.max(.4,Math.round(i*.5*10)/10),defaultMax:Math.min(5,Math.round(i*3*10)/10),defaultSteps:20,defaultCharge:0,defaultMult:e%2===0?1:2,basisOptions:["STO-3G"],generateXYZ:c=>`2
${l} R=${c.toFixed(3)}
${t.symbol}  0.0  0.0  0.0
${a.symbol}  0.0  0.0  ${c.toFixed(6)}`}}const Ht=new Map;async function te(t){const a=Ht.get(t);if(a)return a;const r=`/GANSU-Lite/basis/${t.toLowerCase()}.gbs`,o=await fetch(r);if(!o.ok)throw new Error(`Failed to load basis set: ${t}`);const i=await o.text(),e=Yt.fromGBS(i);return Ht.set(t,e),e}function J(t){return t==="rhf"?ct()?"#00d4ff":"#0077cc":ct()?"#ff8844":"#cc4400"}function pt(){const t=[];return K&&t.push({key:"rhf",label:_(A,!0),color:J("rhf")}),P&&t.push({key:"uhf",label:_(A,!1),color:J("uhf")}),t}const E=document.getElementById("app");function Z(){const t=F.paramType==="angle"||F.paramType==="dihedral"?"°":"Å";E.innerHTML=`
    <div class="opt-page">
      ${Ut("optimize")}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${H("opt.scenario")}</h2>
          <div class="opt-scenario-grid" id="scenario-grid"></div>

          <h2>${H("opt.basis")}</h2>
          <div class="opt-basis-row" id="basis-row"></div>

          <h2>${H("opt.method")}</h2>
          <div class="opt-method-row">
            <button id="btn-rhf" class="opt-method-btn${K?" active":""}" style="--method-color:${J("rhf")}">${_(A,!0)}</button>
            ${F.category!=="geometry"?`<button id="btn-uhf" class="opt-method-btn${P?" active":""}" style="--method-color:${J("uhf")}">${_(A,!1)}</button>`:""}
          </div>

          <div class="opt-charge-row">
            <label>
              <span>${H("opt.charge")}</span>
              <input id="inp-charge" type="number" min="-3" max="3" value="${U}" />
            </label>
            <label>
              <span>${H("opt.mult")}</span>
              <input id="inp-mult" type="number" min="1" max="5" value="${ot}" />
            </label>
          </div>

          <h2>${H("opt.param")}</h2>
          <div class="opt-param-form">
            <label>
              <span>${H("opt.min")}</span>
              <input id="param-min" type="number" step="0.1" value="${I}" />
              <span class="opt-unit">${t}</span>
            </label>
            <label>
              <span>${H("opt.max")}</span>
              <input id="param-max" type="number" step="0.1" value="${Y}" />
              <span class="opt-unit">${t}</span>
            </label>
            <label>
              <span>${H("opt.steps")}</span>
              <input id="param-steps" type="number" min="3" max="50" value="${W}" />
            </label>
          </div>

          <div class="theory-row" style="margin:8px 0;display:flex;align-items:center;gap:8px;font-size:0.9rem;">
            <span>Theory:</span>${Zt("theory-sel",A,"",qt)}
          </div>

          <button id="run-btn" class="opt-run-btn" ${X?"disabled":""}>
            ${X?H("opt.running"):H("opt.run")}
          </button>
          ${X?`<button id="stop-btn" class="opt-stop-btn">${H("opt.stop")}</button>`:""}

          <div id="progress-area"></div>
          <div id="result-summary"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="mol-vis"></div>
          <div id="graph-container">
            ${!G&&!X?`<p class="opt-hint">${H("opt.waiting")}</p>`:""}
          </div>
          <div id="scan-slider-area"></div>
        </div>
      </div>
    </div>`,le();const a=E.querySelector("#scenario-grid"),r=[{key:"dissociation",labelKey:"opt.catDissociation"},{key:"geometry",labelKey:"opt.catGeometry"}];for(const c of r){const n=zt.filter(u=>u.category===c.key);if(n.length===0)continue;const p=document.createElement("div");p.className="opt-category-header",p.textContent=H(c.labelKey),a.appendChild(p);const $=document.createElement("div");$.className="opt-category-row";for(const u of n){const k=document.createElement("div");k.className="opt-scenario-card"+(u.id===F.id?" selected":""),k.innerHTML=`<strong>${H(u.labelKey)}</strong><span>${H(u.descKey)}</span>`,k.addEventListener("click",()=>{X||(F=u,I=u.defaultMin,Y=u.defaultMax,W=u.defaultSteps,U=u.defaultCharge,ot=u.defaultMult,u.category==="geometry"&&(P=!1,K=!0),V=null,D=null,C=[],G=!1,Z())}),$.appendChild(k)}a.appendChild($)}{const c=document.createElement("div");c.className="opt-category-header",c.textContent=H("opt.catCustom"),a.appendChild(c);const n=document.createElement("div");n.className="opt-category-row";const p=F.id.startsWith("custom_"),$=document.createElement("div");$.className="opt-scenario-card"+(p?" selected":"");const u=V&&D?`${V.symbol} + ${D.symbol}`:H("opt.customSelect");$.innerHTML=`<strong>${u}</strong><span>${H("opt.catCustom")}</span>`,$.addEventListener("click",()=>{X||Wt((k,x,v,m)=>{V={symbol:k,z:x},D={symbol:v,z:m},o()},V,D)}),n.appendChild($),a.appendChild(n)}function o(){if(!V||!D)return;const c=Qt(V,D);F=c,I=c.defaultMin,Y=c.defaultMax,W=c.defaultSteps,U=c.defaultCharge,ot=c.defaultMult,(V.z+D.z-U)%2!==0?(K=!1,P=!0):(K=!0,P=!0),C=[],G=!1,Z()}const i=E.querySelector("#basis-row"),e=F.basisOptions??Jt;e.includes(et)||(et=e[0]);for(const c of e){const n=document.createElement("button");n.className="opt-basis-btn"+(c===et?" selected":""),n.textContent=c,n.addEventListener("click",()=>{X||(et=c,C=[],G=!1,Z())}),i.appendChild(n)}E.querySelector("#btn-rhf").addEventListener("click",()=>{X||K&&!P||(K=!K,C=[],G=!1,Z())}),E.querySelector("#btn-uhf")?.addEventListener("click",()=>{X||P&&!K||(P=!P,C=[],G=!1,Z())}),E.querySelector("#inp-charge").addEventListener("change",c=>{U=parseInt(c.target.value,10)}),E.querySelector("#inp-mult").addEventListener("change",c=>{ot=parseInt(c.target.value,10)}),E.querySelector("#nav-theme").addEventListener("click",()=>{Vt(),Z()}),E.querySelector("#nav-lang").addEventListener("click",()=>{Dt(),Z()});const l=E.querySelector("#theory-sel");l&&l.addEventListener("change",()=>{A=l.value,Z()}),E.querySelector("#run-btn").addEventListener("click",oe),E.querySelector("#stop-btn")?.addEventListener("click",()=>{at=!0}),E.querySelector("#param-min").addEventListener("change",c=>{I=parseFloat(c.target.value)}),E.querySelector("#param-max").addEventListener("change",c=>{Y=parseFloat(c.target.value)}),E.querySelector("#param-steps").addEventListener("change",c=>{const n=c.target;let p=parseInt(n.value,10);(isNaN(p)||p<3)&&(p=3),p>50&&(p=50),n.value=String(p),W=p}),C.length>0&&rt(C,G),G?(ae(),ne(),It()):st([{param:(I+Y)/2}])}function ee(t,a,r,o){if(o<=0||o>=r)return;const i=Math.PI/4,e=Math.cos(i),l=Math.sin(i),c=o-1,n=o,p=a.clone();for(let x=0;x<r;x++)p.set(x,c,e*a.get(x,c)+l*a.get(x,n));const $=a.clone();for(let x=0;x<r;x++)$.set(x,c,e*a.get(x,c)-l*a.get(x,n));const u=new Mt(r,r),k=new Mt(r,r);for(let x=0;x<r;x++)for(let v=0;v<r;v++){let m=0,b=0;for(let S=0;S<o;S++)m+=p.get(x,S)*p.get(v,S),b+=$.get(x,S)*$.get(v,S);u.set(x,v,m),k.set(x,v,b)}t.setInitialDensityGuessAlphaBeta(u,k)}async function oe(){if(X)return;if(A!=="HF"&&P&&K&&W>=15){const n=Math.round(W*2*8/60);if(!confirm(`This scan combines DFT × RHF+UHF × ${W} points → estimated ~${n} min.
Proceed?`))return}X=!0,at=!1,C=[],G=!1,Z();const t=Math.max(3,Math.min(50,W)),a=[];for(let n=0;n<=t;n++)a.push(I+(Y-I)*n/t);const r=Math.floor((ot-1)/2),o=performance.now(),i=(K?a.length:0)+(P?a.length:0);let e=0;try{const n=await te(et);let p=null;if(K)for(let $=0;$<a.length&&!at;$++){const u=a[$],k=nt(F.generateXYZ(u)),x=new it(k,n,U,r),v=await lt(x,n,A,"RHF");p&&v.setInitialDensityGuess(p);const m=await v.solve({eriBackend:"js"});p=v.density.clone(),C.push({param:u,rhf:m}),e++,Ct(e,i),rt(C,!1);{const b=[{param:u,label:_(A,!0),color:J("rhf")}];P&&b.push({param:(I+Y)/2,label:_(A,!1),color:J("uhf")}),st(b)}await new Promise(b=>setTimeout(b,0))}else for(const $ of a)C.push({param:$,rhf:NaN});if(P&&!at){let $=null,u=null,k=null;if((F.defaultCharge!==void 0?(a.length>0?nt(F.generateXYZ(a[0])).reduce((m,b)=>m+b.atomicNumber,0):0)-U:0)%2===0){const m=nt(F.generateXYZ(a[a.length-1])),b=new it(m,n,U,r),S=await lt(b,n,A,"RHF");p&&S.setInitialDensityGuess(p),await S.solve({eriBackend:"js"}),k=S.coefficients}for(let m=a.length-1;m>=0&&!at;m--){const b=a[m],S=nt(F.generateXYZ(b));try{const y=new it(S,n,U,r),g=await lt(y,n,A,"UHF");$&&u?g.setInitialDensityGuessAlphaBeta($,u):k&&ee(g,k,y.numBasis,y.numAlphaSpins);let N=1/0;const L=await g.solve({eriBackend:"js",onIteration:(O,T,q)=>{N=Math.abs(q)}});N<1e-6?($=g.densityAlphaMatrix.clone(),u=g.densityBetaMatrix.clone(),C[m].uhf=L):console.warn(`UHF not converged at param=${b}, deltaE=${N}`)}catch(y){console.warn("UHF error at param",b,y)}e++,Ct(e,i),rt(C,!1);{const y=[];if(K){const g=C.filter(L=>isFinite(L.rhf)),N=g.length>0?g.reduce((L,O)=>L.rhf<O.rhf?L:O):null;y.push({param:N?.param??(I+Y)/2,label:_(A,!0),color:J("rhf")})}y.push({param:b,label:_(A,!1),color:J("uhf")}),st(y)}await new Promise(y=>setTimeout(y,0))}}}catch(n){console.error("Scan error:",n)}Lt=performance.now()-o,X=!1,G=!0;const l=K?"rhf":"uhf";let c=0;for(let n=1;n<C.length;n++){const p=C[n][l],$=C[c][l];p!=null&&isFinite(p)&&($==null||!isFinite($)||p<$)&&(c=n)}B=c,Z()}function Ct(t,a){const r=E.querySelector("#progress-area");if(!r)return;const o=(t/a*100).toFixed(0),i=H("opt.progress").replace("{n}",String(t)).replace("{total}",String(a));r.innerHTML=`
    <div class="opt-progress">
      <div class="opt-progress-bar" style="width:${o}%"></div>
    </div>
    <p class="opt-progress-text">${i}</p>`}function ae(){const t=E.querySelector("#result-summary");if(!t||C.length===0)return;const a=F.paramType==="angle"||F.paramType==="dihedral"?"°":"Å",r=(Lt/1e3).toFixed(1),o=pt();let i="";for(const e of o){const l=C.filter(n=>n[e.key]!=null&&isFinite(n[e.key])).map(n=>({param:n.param,e:n[e.key]}));if(l.length===0)continue;let c=l[0];for(const n of l)n.e<c.e&&(c=n);i+=`<tr>
      <td><span class="opt-dot" style="background:${e.color}"></span>${e.label}</td>
      <td>${c.param.toFixed(4)} ${a}</td>
      <td><strong>${c.e.toFixed(8)} Eh</strong></td>
    </tr>`}t.innerHTML=`
    <div class="opt-summary">
      <h3>${H("opt.done")} (${r} s)</h3>
      <table>
        <tr><th></th><th>${H("opt.resultParam")}</th><th>${H("opt.resultEnergy")}</th></tr>
        ${i}
      </table>
    </div>`}function ne(){const t=E.querySelector("#scan-slider-area");if(!t||C.length===0)return;B>=C.length&&(B=0);const a=F.paramType==="angle"||F.paramType==="dihedral",r=a?"°":" Å",o=a?C[B].param.toFixed(1):C[B].param.toFixed(3);t.innerHTML=`
    <div class="walsh-slider">
      <input type="range" id="scan-slider" min="0" max="${C.length-1}" value="${B}" />
      <div class="walsh-slider-label">${o}${r}</div>
    </div>`,E.querySelector("#scan-slider")?.addEventListener("input",i=>{B=parseInt(i.target.value);const e=C[B].param,l=E.querySelector(".walsh-slider-label");l&&(l.textContent=`${a?e.toFixed(1):e.toFixed(3)}${r}`),rt(C,!0),It()})}function It(){if(C.length===0)return;B>=C.length&&(B=0);const t=C[B],a=pt(),r=[];for(const o of a)r.push({param:t.param,label:o.label,color:o.color});r.length>0&&st(r)}const se={H:"#999",He:"#0CC",Li:"#C2C",Be:"#6C0",B:"#F90",C:"#555",N:"#35F",O:"#F22",F:"#9E5"},kt={H:.31,He:.28,Li:1.28,Be:.96,B:.84,C:.76,N:.71,O:.66,F:.57};function Nt(t){if(F.paramType==="height")return{pU:"x",pV:"z",pW:"y"};const a=["x","y","z"].map(r=>({a:r,r:Math.max(...t.map(o=>o[r]))-Math.min(...t.map(o=>o[r]))})).sort((r,o)=>o.r-r.r);return{pU:a[0].a,pV:a[1].a,pW:a[2].a}}function Rt(t,a){const r=a??Nt(t);if(F.paramType==="height"){const o=Math.cos(Math.PI/6),i=Math.sin(Math.PI/6);return t.map(e=>({u:e.x*o+e.y*i,v:e.z}))}return t.map(o=>({u:o[r.pU],v:o[r.pV]}))}function dt(t){const a=t.split(`
`),r=parseInt(a[0].trim()),o=[];for(let i=2;i<2+r;i++){const e=a[i].trim().split(/\s+/);o.push({sym:e[0],x:+e[1],y:+e[2],z:+e[3]})}return o}function re(t,a,r){let o=[];for(const y of[I,Y])o.push(...dt(F.generateXYZ(y)));const i=Nt(o);let e=[],l=[];for(const y of[I,Y]){const g=Rt(dt(F.generateXYZ(y)),i);e.push(...g.map(N=>N.u)),l.push(...g.map(N=>N.v))}const c=Math.min(...e),n=Math.max(...e),p=Math.min(...l),$=Math.max(...l),u=n-c||.01,k=$-p||.01,x=r,v=r+30,m=Math.min((t-2*r)/u,(a-x-v)/k,120),b=(p+$)/2,S=(x-v)/2;return{scale:m,uMid:(c+n)/2,vMid:b,cySvgOffset:S,axes:i}}function St(t,a,r,o,i){const e=dt(F.generateXYZ(t));if(e.length===0)return"";const l=Et(),c=35,{scale:n,uMid:p,vMid:$,cySvgOffset:u,axes:k}=re(a,r,c),x=a/2,v=r/2+u,m=s=>x+(s-p)*n,b=s=>v-(s-$)*n,S=18,y=Rt(e,k);let g="";for(let s=0;s<e.length;s++)for(let d=s+1;d<e.length;d++){const f=e[s].x-e[d].x,h=e[s].y-e[d].y,w=e[s].z-e[d].z,M=Math.sqrt(f*f+h*h+w*w),R=1.4*((kt[e[s].sym]??.7)+(kt[e[d].sym]??.7));M<R&&(g+=`<line x1="${m(y[s].u)}" y1="${b(y[s].v)}" x2="${m(y[d].u)}" y2="${b(y[d].v)}" stroke="${l.grid}" stroke-width="4" stroke-linecap="round"/>`)}const N=e.map((s,d)=>d).sort((s,d)=>e[s][k.pW]-e[d][k.pW]);for(const s of N){const d=m(y[s].u),f=b(y[s].v),h=se[e[s].sym]??"#888";g+=`<circle cx="${d}" cy="${f}" r="${S}" fill="${h}" stroke="${l.axis}" stroke-width="1.2"/>`,g+=`<text x="${d}" y="${f}" text-anchor="middle" dy="0.38em" font-size="11" font-weight="bold" fill="#fff" stroke="#0003" stroke-width="0.3">${e[s].sym}</text>`}const L=F.paramType==="angle"||F.paramType==="dihedral",O=L?"°":"Å",q=`${L?t.toFixed(1):t.toFixed(3)} ${O}`,z=i??l.accent;if(F.paramType==="bond"){const[s,d]=F.paramAtomPair??[0,1],f=m(y[s].u),h=m(y[d].u),w=Math.max(b(y[s].v),b(y[d].v))+S+10;g+=`<line x1="${f}" y1="${w}" x2="${h}" y2="${w}" stroke="${z}" stroke-width="1.2"/>`,g+=`<line x1="${f}" y1="${w-4}" x2="${f}" y2="${w+4}" stroke="${z}" stroke-width="1.2"/>`,g+=`<line x1="${h}" y1="${w-4}" x2="${h}" y2="${w+4}" stroke="${z}" stroke-width="1.2"/>`;const M=h>f?1:-1;g+=`<polygon points="${f},${w} ${f+M*6},${w-3} ${f+M*6},${w+3}" fill="${z}"/>`,g+=`<polygon points="${h},${w} ${h-M*6},${w-3} ${h-M*6},${w+3}" fill="${z}"/>`,g+=`<text x="${(f+h)/2}" y="${w+14}" text-anchor="middle" font-size="11" font-weight="600" fill="${z}">${q}</text>`}else if(F.paramType==="angle"){const s=m(y[0].u),d=b(y[0].v),f=m(y[1].u),h=b(y[1].v),w=m(y[2].u),M=b(y[2].v),R=Math.atan2(h-d,f-s),Q=Math.atan2(M-d,w-s),j=22,ut=s+j*Math.cos(R),ft=d+j*Math.sin(R),ht=s+j*Math.cos(Q),xt=d+j*Math.sin(Q),Kt=(ut-s)*(xt-d)-(ft-d)*(ht-s)>0?1:0;g+=`<path d="M ${ut} ${ft} A ${j} ${j} 0 0 ${Kt} ${ht} ${xt}" fill="none" stroke="${z}" stroke-width="1.5"/>`;let tt=Q-R;tt>Math.PI&&(tt-=2*Math.PI),tt<-Math.PI&&(tt+=2*Math.PI);const mt=R+tt/2,gt=j+14;g+=`<text x="${s+gt*Math.cos(mt)}" y="${d+gt*Math.sin(mt)}" text-anchor="middle" dy="0.35em" font-size="11" font-weight="600" fill="${z}">${q}</text>`;const yt=e[1].x-e[0].x,$t=e[1].y-e[0].y,bt=e[1].z-e[0].z,At=Math.sqrt(yt*yt+$t*$t+bt*bt),Bt=(s+f)/2,Pt=(d+h)/2,vt=Math.atan2(h-d,f-s),Ot=Bt+12*Math.cos(vt+Math.PI/2),Xt=Pt+12*Math.sin(vt+Math.PI/2);g+=`<text x="${Ot}" y="${Xt}" text-anchor="middle" dy="0.35em" font-size="9" fill="${l.dim}">${At.toFixed(2)} Å</text>`}else if(F.paramType==="height"){const s=m(y[0].u),d=b(y[0].v),f=b(0),h=Math.min(...y.slice(1).map(M=>m(M.u)))-20,w=Math.max(...y.slice(1).map(M=>m(M.u)))+20;if(g+=`<line x1="${h}" y1="${f}" x2="${w}" y2="${f}" stroke="${l.dim}" stroke-width="1" stroke-dasharray="5,3"/>`,g+=`<text x="${w+4}" y="${f}" dy="0.35em" font-size="9" fill="${l.dim}">H₃ plane</text>`,Math.abs(d-f)>8){const M=s>x?-24:24;g+=`<line x1="${s+M}" y1="${d}" x2="${s+M}" y2="${f}" stroke="${z}" stroke-width="1.5"/>`,g+=`<line x1="${s+M-5}" y1="${d}" x2="${s+M+5}" y2="${d}" stroke="${z}" stroke-width="1.5"/>`,g+=`<line x1="${s+M-5}" y1="${f}" x2="${s+M+5}" y2="${f}" stroke="${z}" stroke-width="1.5"/>`;const R=d<f?1:-1;g+=`<polygon points="${s+M},${d} ${s+M-3},${d+R*6} ${s+M+3},${d+R*6}" fill="${z}"/>`,g+=`<polygon points="${s+M},${f} ${s+M-3},${f-R*6} ${s+M+3},${f-R*6}" fill="${z}"/>`;const Q=s+M+(M>0?10:-10),j=M>0?"start":"end";g+=`<text x="${Q}" y="${(d+f)/2}" text-anchor="${j}" dy="0.35em" font-size="11" font-weight="600" fill="${z}">h = ${q}</text>`}else g+=`<text x="${s}" y="${f+S+14}" text-anchor="middle" font-size="11" font-weight="600" fill="${z}">h = ${q}</text>`}else if(F.paramType==="dihedral"){const s=Math.max(...y.map(d=>b(d.v)))+S+12;g+=`<text x="${x}" y="${s}" text-anchor="middle" font-size="11" font-weight="600" fill="${z}">φ = ${q}</text>`}return o&&(g+=`<text x="${a/2}" y="14" text-anchor="middle" font-size="11" font-weight="700" fill="${z}">${o}</text>`),g}function st(t){const a=E.querySelector("#mol-vis");if(!a||t.length===0)return;const r=520,o=200;if(t.length===1){const i=t[0],e=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${r} ${o}" style="width:100%;max-width:${r}px;display:block;margin:0 auto;">${St(i.param,r,o,i.label,i.color)}</svg>`;a.innerHTML=e}else{const i=Math.floor(r/t.length);let e='<div style="display:flex;gap:4px;justify-content:center;">';for(const l of t)e+=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${i} ${o}" style="flex:1;max-width:${i}px;">${St(l.param,i,o,l.label,l.color)}</svg>`;e+="</div>",a.innerHTML=e}}function rt(t,a){const r=E.querySelector("#graph-container");if(!r||t.length<1)return;const o=Et(),e={bond:H("opt.xBond"),angle:H("opt.xAngle"),height:H("opt.xHeight"),dihedral:H("opt.xDihedral")}[F.paramType],l=pt(),c=520,n=360,p=72,$=24,u=36,k=44,x=c-p-$,v=n-u-k,m=[];for(const s of l)for(const d of t){const f=d[s.key];f!=null&&isFinite(f)&&m.push(f)}m.length===0&&m.push(0);const b=Math.min(...m),S=Math.max(...m),y=(S-b)*.1||.01,g=b-y,N=S+y,L=s=>p+(s-I)/(Y-I||1)*x,O=s=>u+v-(s-g)/(N-g||1)*v;let T=`<svg width="${c}" height="${n}" xmlns="http://www.w3.org/2000/svg" style="display:block;width:100%;height:auto;max-width:${c}px;" viewBox="0 0 ${c} ${n}">`;T+=`<rect x="${p}" y="${u}" width="${x}" height="${v}" fill="${o.surface}" rx="2"/>`;for(let s=0;s<=5;s++){const d=g+(N-g)*s/5,f=O(d);T+=`<line x1="${p}" y1="${f}" x2="${p+x}" y2="${f}" stroke="${o.grid}" stroke-width="0.5"/>`,T+=`<text x="${p-6}" y="${f+3}" text-anchor="end" font-size="9" font-family="monospace" fill="${o.dim}">${d.toFixed(4)}</text>`}for(let s=0;s<=5;s++){const d=I+(Y-I)*s/5,f=L(d);T+=`<line x1="${f}" y1="${u}" x2="${f}" y2="${u+v}" stroke="${o.grid}" stroke-width="0.5"/>`,T+=`<text x="${f}" y="${u+v+14}" text-anchor="middle" font-size="9" fill="${o.dim}">${d.toFixed(2)}</text>`}T+=`<line x1="${p}" y1="${u}" x2="${p}" y2="${u+v}" stroke="${o.axis}" stroke-width="1"/>`,T+=`<line x1="${p}" y1="${u+v}" x2="${p+x}" y2="${u+v}" stroke="${o.axis}" stroke-width="1"/>`;for(const s of l){const d=t.filter(h=>h[s.key]!=null&&isFinite(h[s.key])).map(h=>({x:h.param,y:h[s.key]}));if(d.length<2){for(const h of d)T+=`<circle cx="${L(h.x).toFixed(1)}" cy="${O(h.y).toFixed(1)}" r="3" fill="${s.color}"/>`;continue}let f="";for(let h=0;h<d.length;h++){const w=L(d[h].x),M=O(d[h].y);f+=h===0?`M${w.toFixed(1)},${M.toFixed(1)}`:` L${w.toFixed(1)},${M.toFixed(1)}`}T+=`<path d="${f}" fill="none" stroke="${s.color}" stroke-width="2" stroke-linejoin="round"/>`;for(const h of d)T+=`<circle cx="${L(h.x).toFixed(1)}" cy="${O(h.y).toFixed(1)}" r="3" fill="${s.color}" stroke="${o.surface}" stroke-width="0.8"/>`;if(d.length>=3&&a){let h=d[0];for(const R of d)R.y<h.y&&(h=R);const w=L(h.x),M=O(h.y);T+=`<line x1="${w}" y1="${M}" x2="${w}" y2="${u+v}" stroke="${s.color}" stroke-width="1" stroke-dasharray="3,3" opacity="0.5"/>`,T+=ie(w,M,6,s.color)}}if(a&&B>=0&&B<t.length){const s=t[B],d=L(s.param),f=ct()?"#ffd700":"#cc8800";T+=`<line x1="${d}" y1="${u}" x2="${d}" y2="${u+v}" stroke="${f}" stroke-width="1" stroke-dasharray="3,2" opacity="0.6"/>`;for(const h of l){const w=s[h.key];w!=null&&isFinite(w)&&(T+=`<circle cx="${d}" cy="${O(w)}" r="5" fill="${h.color}" stroke="${f}" stroke-width="2"/>`)}}const q=p+8;let z=u+14;for(const s of l)T+=`<rect x="${q}" y="${z-7}" width="10" height="3" rx="1" fill="${s.color}"/>`,T+=`<text x="${q+14}" y="${z-3}" font-size="9" font-weight="600" fill="${o.dim}">${s.label}</text>`,z+=14;T+=`<text x="${p+x/2}" y="20" text-anchor="middle" font-size="12" font-weight="600" fill="${o.titleSvg}">${H("opt.graphTitle")}</text>`,T+=`<text x="${p+x/2}" y="${n-4}" text-anchor="middle" font-size="10" fill="${o.dim}">${e}</text>`,T+=`<text x="14" y="${u+v/2}" text-anchor="middle" font-size="10" fill="${o.dim}" transform="rotate(-90,14,${u+v/2})">${H("opt.yEnergy")}</text>`,T+="</svg>",r.innerHTML=T}function ie(t,a,r,o){const i=[];for(let e=0;e<10;e++){const l=Math.PI/2+e*Math.PI/5,c=e%2===0?r:r*.4;i.push(`${(t+c*Math.cos(l)).toFixed(1)},${(a-c*Math.sin(l)).toFixed(1)}`)}return`<polygon points="${i.join(" ")}" fill="${o}" stroke="${o}" stroke-width="0.5"/>`}let Tt=!1;function le(){if(Tt)return;Tt=!0;const t=document.createElement("style");t.textContent=`
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
  `,document.head.appendChild(t)}Gt();jt();Z();
