import{t as H,i as Ot,c as Xt,r as Yt,a as Zt,b as qt,d as it,g as St}from"./nav-Lc1uTV-6.js";import{p as ot}from"./parseXYZ-DTmoiCYR.js";import{t as st,u as rt,d as Gt,M as bt}from"./builder-B-MUn5dT.js";const vt={1:.31,2:.28,3:1.28,4:.96,5:.84,6:.76,7:.71,8:.66,9:.57,10:.58,11:1.66,12:1.41,13:1.21,14:1.11,15:1.07,16:1.05,17:1.02,18:1.06,19:2.03,20:1.76,21:1.7,22:1.6,23:1.53,24:1.39,25:1.39,26:1.32,27:1.26,28:1.24,29:1.32,30:1.22,31:1.22,32:1.2,33:1.19,34:1.2,35:1.2,36:1.16},jt=[[0,0,1,"H"],[0,17,2,"He"],[1,0,3,"Li"],[1,1,4,"Be"],[1,12,5,"B"],[1,13,6,"C"],[1,14,7,"N"],[1,15,8,"O"],[1,16,9,"F"],[1,17,10,"Ne"],[2,0,11,"Na"],[2,1,12,"Mg"],[2,12,13,"Al"],[2,13,14,"Si"],[2,14,15,"P"],[2,15,16,"S"],[2,16,17,"Cl"],[2,17,18,"Ar"],[3,0,19,"K"],[3,1,20,"Ca"],[3,2,21,"Sc"],[3,3,22,"Ti"],[3,4,23,"V"],[3,5,24,"Cr"],[3,6,25,"Mn"],[3,7,26,"Fe"],[3,8,27,"Co"],[3,9,28,"Ni"],[3,10,29,"Cu"],[3,11,30,"Zn"],[3,12,31,"Ga"],[3,13,32,"Ge"],[3,14,33,"As"],[3,15,34,"Se"],[3,16,35,"Br"],[3,17,36,"Kr"]];function Mt(t,n,r=18){const o=document.createElement("div");o.className="pte-grid";for(const[l,e,a,p]of jt){const i=document.createElement("div");i.className="pte-cell",i.style.gridRow=String(l+1),i.style.gridColumn=String(e+1);const d=a>r;d&&i.classList.add("pte-disabled");const f=a<=2||a<=4||a>=11&&a<=12||a>=19&&a<=20?"pte-s":a>=21&&a<=30?"pte-d":"pte-p";i.classList.add(f),a===n&&i.classList.add("pte-selected"),i.innerHTML=`<span class="pte-z">${a}</span><span class="pte-sym">${p}</span>`,d||i.addEventListener("click",()=>t(p,a)),o.appendChild(i)}return o}function Ut(t,n,r){let o=n??null,l=r??null;const e=document.createElement("div");e.className="pte-overlay";const a=document.createElement("div");a.className="pte-modal pte-dual";function p(){a.innerHTML="";const i=document.createElement("div");i.className="pte-tables";const d=document.createElement("div");d.className="pte-col";const f=document.createElement("div");f.className="pte-title",f.textContent=`${H("opt.customAtomA")}${o?": "+o.symbol:""}`,d.appendChild(f),d.appendChild(Mt((b,C)=>{o={symbol:b,z:C},p()},o?.z,20)),i.appendChild(d);const x=document.createElement("div");x.className="pte-col";const T=document.createElement("div");T.className="pte-title",T.textContent=`${H("opt.customAtomB")}${l?": "+l.symbol:""}`,x.appendChild(T),x.appendChild(Mt((b,C)=>{l={symbol:b,z:C},p()},l?.z,12)),i.appendChild(x),a.appendChild(i);const m=document.createElement("div");m.className="pte-btn-row";const v=document.createElement("button");v.className="pte-btn pte-btn-cancel",v.textContent="Cancel",v.addEventListener("click",()=>e.remove()),m.appendChild(v);const $=document.createElement("button");$.className="pte-btn pte-btn-ok",$.textContent="OK",$.disabled=!o||!l,$.addEventListener("click",()=>{o&&l&&(e.remove(),t(o.symbol,o.z,l.symbol,l.z))}),m.appendChild($),a.appendChild(m)}p(),e.appendChild(a),e.addEventListener("click",i=>{i.target===e&&e.remove()}),document.body.appendChild(e)}const Tt=[{id:"h2",category:"dissociation",labelKey:"opt.scenH2",descKey:"opt.descH2",paramType:"bond",defaultMin:.4,defaultMax:3,defaultSteps:20,defaultCharge:0,defaultMult:1,basisOptions:["3-21G","6-31G"],generateXYZ:t=>`2
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
H  0.0  0.0  ${(t+1.06).toFixed(6)}`},{id:"c2h4",category:"dissociation",labelKey:"opt.scenC2H4",descKey:"opt.descC2H4",paramType:"bond",defaultMin:1,defaultMax:2.5,defaultSteps:20,defaultCharge:0,defaultMult:1,basisOptions:["STO-3G"],generateXYZ:t=>{const r=121.7*Math.PI/180,o=1.08*Math.sin(r-Math.PI/2)*-1,l=1.08*Math.cos(r-Math.PI/2),e=Math.abs(o);return["6",`C2H4 CC=${t.toFixed(3)}`,"C   0.000000  0.000000  0.000000",`C   0.000000  0.000000  ${t.toFixed(6)}`,`H   ${e.toFixed(6)}  0.000000  ${(-l).toFixed(6)}`,`H   ${(-e).toFixed(6)}  0.000000  ${(-l).toFixed(6)}`,`H   ${e.toFixed(6)}  0.000000  ${(t+l).toFixed(6)}`,`H   ${(-e).toFixed(6)}  0.000000  ${(t+l).toFixed(6)}`].join(`
`)}},{id:"h2o_bond",category:"dissociation",labelKey:"opt.scenH2Obond",descKey:"opt.descH2Obond",paramType:"bond",defaultMin:.5,defaultMax:3,defaultSteps:20,defaultCharge:0,defaultMult:1,basisOptions:["3-21G","6-31G"],paramAtomPair:[0,2],generateXYZ:t=>{const n=Math.sin(52*Math.PI/180),r=Math.cos(52*Math.PI/180);return["3",`H2O R=${t.toFixed(3)}`,"O   0.000000  0.000000  0.000000",`H   ${(.96*n).toFixed(6)}  0.000000  ${(.96*r).toFixed(6)}`,`H   ${(-t*n).toFixed(6)}  0.000000  ${(t*r).toFixed(6)}`].join(`
`)}},{id:"h2o",category:"geometry",labelKey:"opt.scenH2O",descKey:"opt.descH2O",paramType:"angle",defaultMin:80,defaultMax:180,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>{const o=t*Math.PI/180/2,l=.96*Math.sin(o),e=.96*Math.cos(o);return`3
H2O angle=${t.toFixed(1)}
O  0.0  0.0  0.0
H  ${l.toFixed(6)}  0.0  ${e.toFixed(6)}
H  ${(-l).toFixed(6)}  0.0  ${e.toFixed(6)}`}},{id:"beh2",category:"geometry",labelKey:"opt.scenBeH2",descKey:"opt.descBeH2",paramType:"angle",defaultMin:90,defaultMax:270,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>{const o=t*Math.PI/180/2,l=1.33*Math.sin(o),e=1.33*Math.cos(o);return`3
BeH2 angle=${t.toFixed(1)}
Be  0.0  0.0  0.0
H   ${l.toFixed(6)}  0.0  ${e.toFixed(6)}
H   ${(-l).toFixed(6)}  0.0  ${e.toFixed(6)}`}},{id:"nh3",category:"geometry",labelKey:"opt.scenNH3",descKey:"opt.descNH3",paramType:"height",defaultMin:0,defaultMax:.5,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>{const r=Math.sqrt(Math.max(0,1.024144-t*t)),o=Math.sqrt(3)/2;return["4",`NH3 h=${t.toFixed(3)}`,`N   0.000000  0.000000  ${t.toFixed(6)}`,`H   ${r.toFixed(6)}  0.000000  0.000000`,`H   ${(-r/2).toFixed(6)}  ${(r*o).toFixed(6)}  0.000000`,`H   ${(-r/2).toFixed(6)}  ${(-r*o).toFixed(6)}  0.000000`].join(`
`)}},{id:"ch2",category:"geometry",labelKey:"opt.scenCH2",descKey:"opt.descCH2",paramType:"angle",defaultMin:90,defaultMax:180,defaultSteps:20,defaultCharge:0,defaultMult:1,generateXYZ:t=>{const o=t*Math.PI/180/2,l=1.08*Math.sin(o),e=1.08*Math.cos(o);return`3
CH2 angle=${t.toFixed(1)}
C  0.0  0.0  0.0
H  ${l.toFixed(6)}  0.0  ${e.toFixed(6)}
H  ${(-l).toFixed(6)}  0.0  ${e.toFixed(6)}`}},{id:"fhf",category:"geometry",labelKey:"opt.scenFHF",descKey:"opt.descFHF",paramType:"bond",defaultMin:.8,defaultMax:1.5,defaultSteps:20,defaultCharge:-1,defaultMult:1,basisOptions:["STO-3G"],generateXYZ:t=>["3",`FHF- FH=${t.toFixed(3)}`,"F   0.000000  0.000000  0.000000",`H   0.000000  0.000000  ${t.toFixed(6)}`,`F   0.000000  0.000000  ${2.3.toFixed(6)}`].join(`
`)},{id:"h3plus",category:"geometry",labelKey:"opt.scenH3plus",descKey:"opt.descH3plus",paramType:"bond",defaultMin:.5,defaultMax:2.5,defaultSteps:20,defaultCharge:1,defaultMult:1,generateXYZ:t=>{const n=Math.sqrt(3)/2;return["3",`H3+ R=${t.toFixed(3)}`,"H   0.000000  0.000000  0.000000",`H   ${t.toFixed(6)}  0.000000  0.000000`,`H   ${(t/2).toFixed(6)}  ${(t*n).toFixed(6)}  0.000000`].join(`
`)}},{id:"c2h6",category:"geometry",labelKey:"opt.scenC2H6",descKey:"opt.descC2H6",paramType:"dihedral",defaultMin:0,defaultMax:120,defaultSteps:12,defaultCharge:0,defaultMult:1,basisOptions:["STO-3G"],generateXYZ:t=>{const o=109.47*Math.PI/180,l=1.09*Math.cos(Math.PI-o),e=1.09*Math.sin(Math.PI-o),a=["8",`C2H6 dih=${t.toFixed(1)}`];a.push("C   0.000000  0.000000  0.000000"),a.push(`C   0.000000  0.000000  ${1.54.toFixed(6)}`);const p=t*Math.PI/180;for(let i=0;i<3;i++){const d=i*2*Math.PI/3;a.push(`H   ${(e*Math.cos(d)).toFixed(6)}  ${(e*Math.sin(d)).toFixed(6)}  ${(-l).toFixed(6)}`)}for(let i=0;i<3;i++){const d=i*2*Math.PI/3+p;a.push(`H   ${(e*Math.cos(d)).toFixed(6)}  ${(e*Math.sin(d)).toFixed(6)}  ${(1.54+l).toFixed(6)}`)}return a.join(`
`)}}],Vt=["STO-3G","3-21G","6-31G"];let F=Tt[0],J="STO-3G",A=!0,O=!1,j=0,Q=1,R=F.defaultMin,X=F.defaultMax,tt=F.defaultSteps,P=!1,et=!1,k=[],Z=!1,zt=0,K=0,U=null,V=null;function _t(t,n){const r=vt[t.z]??1,o=vt[n.z]??1,l=r+o,e=t.z+n.z,a=t.symbol===n.symbol?`${t.symbol}₂`:`${t.symbol}${n.symbol}`;return{id:`custom_${t.symbol}_${n.symbol}`,category:"dissociation",labelKey:"",descKey:"",paramType:"bond",defaultMin:Math.max(.4,Math.round(l*.5*10)/10),defaultMax:Math.min(5,Math.round(l*3*10)/10),defaultSteps:20,defaultCharge:0,defaultMult:e%2===0?1:2,basisOptions:["STO-3G"],generateXYZ:p=>`2
${a} R=${p.toFixed(3)}
${t.symbol}  0.0  0.0  0.0
${n.symbol}  0.0  0.0  ${p.toFixed(6)}`}}const Ft=new Map;async function Dt(t){const n=Ft.get(t);if(n)return n;const r=`/GANSU-Lite/basis/${t.toLowerCase()}.gbs`,o=await fetch(r);if(!o.ok)throw new Error(`Failed to load basis set: ${t}`);const l=await o.text(),e=Gt.fromGBS(l);return Ft.set(t,e),e}function _(t){return t==="rhf"?it()?"#00d4ff":"#0077cc":it()?"#ff8844":"#cc4400"}function ct(){const t=[];return A&&t.push({key:"rhf",label:"RHF",color:_("rhf")}),O&&t.push({key:"uhf",label:"UHF",color:_("uhf")}),t}const E=document.getElementById("app");function q(){const t=F.paramType==="angle"||F.paramType==="dihedral"?"°":"Å";E.innerHTML=`
    <div class="opt-page">
      ${Yt("optimize")}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${H("opt.scenario")}</h2>
          <div class="opt-scenario-grid" id="scenario-grid"></div>

          <h2>${H("opt.basis")}</h2>
          <div class="opt-basis-row" id="basis-row"></div>

          <h2>${H("opt.method")}</h2>
          <div class="opt-method-row">
            <button id="btn-rhf" class="opt-method-btn${A?" active":""}" style="--method-color:${_("rhf")}">RHF</button>
            ${F.category!=="geometry"?`<button id="btn-uhf" class="opt-method-btn${O?" active":""}" style="--method-color:${_("uhf")}">UHF</button>`:""}
          </div>

          <div class="opt-charge-row">
            <label>
              <span>${H("opt.charge")}</span>
              <input id="inp-charge" type="number" min="-3" max="3" value="${j}" />
            </label>
            <label>
              <span>${H("opt.mult")}</span>
              <input id="inp-mult" type="number" min="1" max="5" value="${Q}" />
            </label>
          </div>

          <h2>${H("opt.param")}</h2>
          <div class="opt-param-form">
            <label>
              <span>${H("opt.min")}</span>
              <input id="param-min" type="number" step="0.1" value="${R}" />
              <span class="opt-unit">${t}</span>
            </label>
            <label>
              <span>${H("opt.max")}</span>
              <input id="param-max" type="number" step="0.1" value="${X}" />
              <span class="opt-unit">${t}</span>
            </label>
            <label>
              <span>${H("opt.steps")}</span>
              <input id="param-steps" type="number" min="3" max="50" value="${tt}" />
            </label>
          </div>

          <button id="run-btn" class="opt-run-btn" ${P?"disabled":""}>
            ${P?H("opt.running"):H("opt.run")}
          </button>
          ${P?`<button id="stop-btn" class="opt-stop-btn">${H("opt.stop")}</button>`:""}

          <div id="progress-area"></div>
          <div id="result-summary"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="mol-vis"></div>
          <div id="graph-container">
            ${!Z&&!P?`<p class="opt-hint">${H("opt.waiting")}</p>`:""}
          </div>
          <div id="scan-slider-area"></div>
        </div>
      </div>
    </div>`,ne();const n=E.querySelector("#scenario-grid"),r=[{key:"dissociation",labelKey:"opt.catDissociation"},{key:"geometry",labelKey:"opt.catGeometry"}];for(const a of r){const p=Tt.filter(f=>f.category===a.key);if(p.length===0)continue;const i=document.createElement("div");i.className="opt-category-header",i.textContent=H(a.labelKey),n.appendChild(i);const d=document.createElement("div");d.className="opt-category-row";for(const f of p){const x=document.createElement("div");x.className="opt-scenario-card"+(f.id===F.id?" selected":""),x.innerHTML=`<strong>${H(f.labelKey)}</strong><span>${H(f.descKey)}</span>`,x.addEventListener("click",()=>{P||(F=f,R=f.defaultMin,X=f.defaultMax,tt=f.defaultSteps,j=f.defaultCharge,Q=f.defaultMult,f.category==="geometry"&&(O=!1,A=!0),U=null,V=null,k=[],Z=!1,q())}),d.appendChild(x)}n.appendChild(d)}{const a=document.createElement("div");a.className="opt-category-header",a.textContent=H("opt.catCustom"),n.appendChild(a);const p=document.createElement("div");p.className="opt-category-row";const i=F.id.startsWith("custom_"),d=document.createElement("div");d.className="opt-scenario-card"+(i?" selected":"");const f=U&&V?`${U.symbol} + ${V.symbol}`:H("opt.customSelect");d.innerHTML=`<strong>${f}</strong><span>${H("opt.catCustom")}</span>`,d.addEventListener("click",()=>{P||Ut((x,T,m,v)=>{U={symbol:x,z:T},V={symbol:m,z:v},o()},U,V)}),p.appendChild(d),n.appendChild(p)}function o(){if(!U||!V)return;const a=_t(U,V);F=a,R=a.defaultMin,X=a.defaultMax,tt=a.defaultSteps,j=a.defaultCharge,Q=a.defaultMult,(U.z+V.z-j)%2!==0?(A=!1,O=!0):(A=!0,O=!0),k=[],Z=!1,q()}const l=E.querySelector("#basis-row"),e=F.basisOptions??Vt;e.includes(J)||(J=e[0]);for(const a of e){const p=document.createElement("button");p.className="opt-basis-btn"+(a===J?" selected":""),p.textContent=a,p.addEventListener("click",()=>{P||(J=a,k=[],Z=!1,q())}),l.appendChild(p)}E.querySelector("#btn-rhf").addEventListener("click",()=>{P||A&&!O||(A=!A,k=[],Z=!1,q())}),E.querySelector("#btn-uhf")?.addEventListener("click",()=>{P||O&&!A||(O=!O,k=[],Z=!1,q())}),E.querySelector("#inp-charge").addEventListener("change",a=>{j=parseInt(a.target.value,10)}),E.querySelector("#inp-mult").addEventListener("change",a=>{Q=parseInt(a.target.value,10)}),E.querySelector("#nav-theme").addEventListener("click",()=>{Zt(),q()}),E.querySelector("#nav-lang").addEventListener("click",()=>{qt(),q()}),E.querySelector("#run-btn").addEventListener("click",Jt),E.querySelector("#stop-btn")?.addEventListener("click",()=>{et=!0}),E.querySelector("#param-min").addEventListener("change",a=>{R=parseFloat(a.target.value)}),E.querySelector("#param-max").addEventListener("change",a=>{X=parseFloat(a.target.value)}),E.querySelector("#param-steps").addEventListener("change",a=>{const p=a.target;let i=parseInt(p.value,10);(isNaN(i)||i<3)&&(i=3),i>50&&(i=50),p.value=String(i),tt=i}),k.length>0&&nt(k,Z),Z?(Qt(),te(),Et()):at([{param:(R+X)/2}])}function Wt(t,n,r,o){if(o<=0||o>=r)return;const l=Math.PI/4,e=Math.cos(l),a=Math.sin(l),p=o-1,i=o,d=n.clone();for(let m=0;m<r;m++)d.set(m,p,e*n.get(m,p)+a*n.get(m,i));const f=n.clone();for(let m=0;m<r;m++)f.set(m,p,e*n.get(m,p)-a*n.get(m,i));const x=new bt(r,r),T=new bt(r,r);for(let m=0;m<r;m++)for(let v=0;v<r;v++){let $=0,b=0;for(let C=0;C<o;C++)$+=d.get(m,C)*d.get(v,C),b+=f.get(m,C)*f.get(v,C);x.set(m,v,$),T.set(m,v,b)}t.setInitialDensityGuessAlphaBeta(x,T)}async function Jt(){if(P)return;P=!0,et=!1,k=[],Z=!1,q();const t=Math.max(3,Math.min(50,tt)),n=[];for(let i=0;i<=t;i++)n.push(R+(X-R)*i/t);const r=Math.floor((Q-1)/2),o=performance.now(),l=(A?n.length:0)+(O?n.length:0);let e=0;try{const i=await Dt(J);let d=null;if(A)for(let f=0;f<n.length&&!et;f++){const x=n[f],T=ot(F.generateXYZ(x)),m=new st(T,i,j,r),v=rt(m,"RHF");d&&v.setInitialDensityGuess(d);const $=await v.solve({eriBackend:"js"});d=v.density.clone(),k.push({param:x,rhf:$}),e++,wt(e,l),nt(k,!1);{const b=[{param:x,label:"RHF",color:_("rhf")}];O&&b.push({param:(R+X)/2,label:"UHF",color:_("uhf")}),at(b)}await new Promise(b=>setTimeout(b,0))}else for(const f of n)k.push({param:f,rhf:NaN});if(O&&!et){let f=null,x=null,T=null;if((F.defaultCharge!==void 0?(n.length>0?ot(F.generateXYZ(n[0])).reduce(($,b)=>$+b.atomicNumber,0):0)-j:0)%2===0){const $=ot(F.generateXYZ(n[n.length-1])),b=new st($,i,j,r),C=rt(b,"RHF");d&&C.setInitialDensityGuess(d),await C.solve({eriBackend:"js"}),T=C.coefficients}for(let $=n.length-1;$>=0&&!et;$--){const b=n[$],C=ot(F.generateXYZ(b));try{const y=new st(C,i,j,r),g=rt(y,"UHF");f&&x?g.setInitialDensityGuessAlphaBeta(f,x):T&&Wt(g,T,y.numBasis,y.numAlphaSpins);let N=1/0;const L=await g.solve({eriBackend:"js",onIteration:(B,S,Y)=>{N=Math.abs(Y)}});N<1e-6?(f=g.densityAlphaMatrix.clone(),x=g.densityBetaMatrix.clone(),k[$].uhf=L):console.warn(`UHF not converged at param=${b}, deltaE=${N}`)}catch(y){console.warn("UHF error at param",b,y)}e++,wt(e,l),nt(k,!1);{const y=[];if(A){const g=k.filter(L=>isFinite(L.rhf)),N=g.length>0?g.reduce((L,B)=>L.rhf<B.rhf?L:B):null;y.push({param:N?.param??(R+X)/2,label:"RHF",color:_("rhf")})}y.push({param:b,label:"UHF",color:_("uhf")}),at(y)}await new Promise(y=>setTimeout(y,0))}}}catch(i){console.error("Scan error:",i)}zt=performance.now()-o,P=!1,Z=!0;const a=A?"rhf":"uhf";let p=0;for(let i=1;i<k.length;i++){const d=k[i][a],f=k[p][a];d!=null&&isFinite(d)&&(f==null||!isFinite(f)||d<f)&&(p=i)}K=p,q()}function wt(t,n){const r=E.querySelector("#progress-area");if(!r)return;const o=(t/n*100).toFixed(0),l=H("opt.progress").replace("{n}",String(t)).replace("{total}",String(n));r.innerHTML=`
    <div class="opt-progress">
      <div class="opt-progress-bar" style="width:${o}%"></div>
    </div>
    <p class="opt-progress-text">${l}</p>`}function Qt(){const t=E.querySelector("#result-summary");if(!t||k.length===0)return;const n=F.paramType==="angle"||F.paramType==="dihedral"?"°":"Å",r=(zt/1e3).toFixed(1),o=ct();let l="";for(const e of o){const a=k.filter(i=>i[e.key]!=null&&isFinite(i[e.key])).map(i=>({param:i.param,e:i[e.key]}));if(a.length===0)continue;let p=a[0];for(const i of a)i.e<p.e&&(p=i);l+=`<tr>
      <td><span class="opt-dot" style="background:${e.color}"></span>${e.label}</td>
      <td>${p.param.toFixed(4)} ${n}</td>
      <td><strong>${p.e.toFixed(8)} Eh</strong></td>
    </tr>`}t.innerHTML=`
    <div class="opt-summary">
      <h3>${H("opt.done")} (${r} s)</h3>
      <table>
        <tr><th></th><th>${H("opt.resultParam")}</th><th>${H("opt.resultEnergy")}</th></tr>
        ${l}
      </table>
    </div>`}function te(){const t=E.querySelector("#scan-slider-area");if(!t||k.length===0)return;K>=k.length&&(K=0);const n=F.paramType==="angle"||F.paramType==="dihedral",r=n?"°":" Å",o=n?k[K].param.toFixed(1):k[K].param.toFixed(3);t.innerHTML=`
    <div class="walsh-slider">
      <input type="range" id="scan-slider" min="0" max="${k.length-1}" value="${K}" />
      <div class="walsh-slider-label">${o}${r}</div>
    </div>`,E.querySelector("#scan-slider")?.addEventListener("input",l=>{K=parseInt(l.target.value);const e=k[K].param,a=E.querySelector(".walsh-slider-label");a&&(a.textContent=`${n?e.toFixed(1):e.toFixed(3)}${r}`),nt(k,!0),Et()})}function Et(){if(k.length===0)return;K>=k.length&&(K=0);const t=k[K],n=ct(),r=[];for(const o of n)r.push({param:t.param,label:o.label,color:o.color});r.length>0&&at(r)}const ee={H:"#999",He:"#0CC",Li:"#C2C",Be:"#6C0",B:"#F90",C:"#555",N:"#35F",O:"#F22",F:"#9E5"},Ht={H:.31,He:.28,Li:1.28,Be:.96,B:.84,C:.76,N:.71,O:.66,F:.57};function Lt(t){if(F.paramType==="height")return{pU:"x",pV:"z",pW:"y"};const n=["x","y","z"].map(r=>({a:r,r:Math.max(...t.map(o=>o[r]))-Math.min(...t.map(o=>o[r]))})).sort((r,o)=>o.r-r.r);return{pU:n[0].a,pV:n[1].a,pW:n[2].a}}function Rt(t,n){const r=n??Lt(t);if(F.paramType==="height"){const o=Math.cos(Math.PI/6),l=Math.sin(Math.PI/6);return t.map(e=>({u:e.x*o+e.y*l,v:e.z}))}return t.map(o=>({u:o[r.pU],v:o[r.pV]}))}function lt(t){const n=t.split(`
`),r=parseInt(n[0].trim()),o=[];for(let l=2;l<2+r;l++){const e=n[l].trim().split(/\s+/);o.push({sym:e[0],x:+e[1],y:+e[2],z:+e[3]})}return o}function oe(t,n,r){let o=[];for(const y of[R,X])o.push(...lt(F.generateXYZ(y)));const l=Lt(o);let e=[],a=[];for(const y of[R,X]){const g=Rt(lt(F.generateXYZ(y)),l);e.push(...g.map(N=>N.u)),a.push(...g.map(N=>N.v))}const p=Math.min(...e),i=Math.max(...e),d=Math.min(...a),f=Math.max(...a),x=i-p||.01,T=f-d||.01,m=r,v=r+30,$=Math.min((t-2*r)/x,(n-m-v)/T,120),b=(d+f)/2,C=(m-v)/2;return{scale:$,uMid:(p+i)/2,vMid:b,cySvgOffset:C,axes:l}}function kt(t,n,r,o,l){const e=lt(F.generateXYZ(t));if(e.length===0)return"";const a=St(),p=35,{scale:i,uMid:d,vMid:f,cySvgOffset:x,axes:T}=oe(n,r,p),m=n/2,v=r/2+x,$=s=>m+(s-d)*i,b=s=>v-(s-f)*i,C=18,y=Rt(e,T);let g="";for(let s=0;s<e.length;s++)for(let c=s+1;c<e.length;c++){const u=e[s].x-e[c].x,h=e[s].y-e[c].y,w=e[s].z-e[c].z,M=Math.sqrt(u*u+h*h+w*w),I=1.4*((Ht[e[s].sym]??.7)+(Ht[e[c].sym]??.7));M<I&&(g+=`<line x1="${$(y[s].u)}" y1="${b(y[s].v)}" x2="${$(y[c].u)}" y2="${b(y[c].v)}" stroke="${a.grid}" stroke-width="4" stroke-linecap="round"/>`)}const N=e.map((s,c)=>c).sort((s,c)=>e[s][T.pW]-e[c][T.pW]);for(const s of N){const c=$(y[s].u),u=b(y[s].v),h=ee[e[s].sym]??"#888";g+=`<circle cx="${c}" cy="${u}" r="${C}" fill="${h}" stroke="${a.axis}" stroke-width="1.2"/>`,g+=`<text x="${c}" y="${u}" text-anchor="middle" dy="0.38em" font-size="11" font-weight="bold" fill="#fff" stroke="#0003" stroke-width="0.3">${e[s].sym}</text>`}const L=F.paramType==="angle"||F.paramType==="dihedral",B=L?"°":"Å",Y=`${L?t.toFixed(1):t.toFixed(3)} ${B}`,z=l??a.accent;if(F.paramType==="bond"){const[s,c]=F.paramAtomPair??[0,1],u=$(y[s].u),h=$(y[c].u),w=Math.max(b(y[s].v),b(y[c].v))+C+10;g+=`<line x1="${u}" y1="${w}" x2="${h}" y2="${w}" stroke="${z}" stroke-width="1.2"/>`,g+=`<line x1="${u}" y1="${w-4}" x2="${u}" y2="${w+4}" stroke="${z}" stroke-width="1.2"/>`,g+=`<line x1="${h}" y1="${w-4}" x2="${h}" y2="${w+4}" stroke="${z}" stroke-width="1.2"/>`;const M=h>u?1:-1;g+=`<polygon points="${u},${w} ${u+M*6},${w-3} ${u+M*6},${w+3}" fill="${z}"/>`,g+=`<polygon points="${h},${w} ${h-M*6},${w-3} ${h-M*6},${w+3}" fill="${z}"/>`,g+=`<text x="${(u+h)/2}" y="${w+14}" text-anchor="middle" font-size="11" font-weight="600" fill="${z}">${Y}</text>`}else if(F.paramType==="angle"){const s=$(y[0].u),c=b(y[0].v),u=$(y[1].u),h=b(y[1].v),w=$(y[2].u),M=b(y[2].v),I=Math.atan2(h-c,u-s),D=Math.atan2(M-c,w-s),G=22,dt=s+G*Math.cos(I),pt=c+G*Math.sin(I),ut=s+G*Math.cos(D),ft=c+G*Math.sin(D),Nt=(dt-s)*(ft-c)-(pt-c)*(ut-s)>0?1:0;g+=`<path d="M ${dt} ${pt} A ${G} ${G} 0 0 ${Nt} ${ut} ${ft}" fill="none" stroke="${z}" stroke-width="1.5"/>`;let W=D-I;W>Math.PI&&(W-=2*Math.PI),W<-Math.PI&&(W+=2*Math.PI);const ht=I+W/2,xt=G+14;g+=`<text x="${s+xt*Math.cos(ht)}" y="${c+xt*Math.sin(ht)}" text-anchor="middle" dy="0.35em" font-size="11" font-weight="600" fill="${z}">${Y}</text>`;const mt=e[1].x-e[0].x,gt=e[1].y-e[0].y,yt=e[1].z-e[0].z,It=Math.sqrt(mt*mt+gt*gt+yt*yt),Kt=(s+u)/2,At=(c+h)/2,$t=Math.atan2(h-c,u-s),Bt=Kt+12*Math.cos($t+Math.PI/2),Pt=At+12*Math.sin($t+Math.PI/2);g+=`<text x="${Bt}" y="${Pt}" text-anchor="middle" dy="0.35em" font-size="9" fill="${a.dim}">${It.toFixed(2)} Å</text>`}else if(F.paramType==="height"){const s=$(y[0].u),c=b(y[0].v),u=b(0),h=Math.min(...y.slice(1).map(M=>$(M.u)))-20,w=Math.max(...y.slice(1).map(M=>$(M.u)))+20;if(g+=`<line x1="${h}" y1="${u}" x2="${w}" y2="${u}" stroke="${a.dim}" stroke-width="1" stroke-dasharray="5,3"/>`,g+=`<text x="${w+4}" y="${u}" dy="0.35em" font-size="9" fill="${a.dim}">H₃ plane</text>`,Math.abs(c-u)>8){const M=s>m?-24:24;g+=`<line x1="${s+M}" y1="${c}" x2="${s+M}" y2="${u}" stroke="${z}" stroke-width="1.5"/>`,g+=`<line x1="${s+M-5}" y1="${c}" x2="${s+M+5}" y2="${c}" stroke="${z}" stroke-width="1.5"/>`,g+=`<line x1="${s+M-5}" y1="${u}" x2="${s+M+5}" y2="${u}" stroke="${z}" stroke-width="1.5"/>`;const I=c<u?1:-1;g+=`<polygon points="${s+M},${c} ${s+M-3},${c+I*6} ${s+M+3},${c+I*6}" fill="${z}"/>`,g+=`<polygon points="${s+M},${u} ${s+M-3},${u-I*6} ${s+M+3},${u-I*6}" fill="${z}"/>`;const D=s+M+(M>0?10:-10),G=M>0?"start":"end";g+=`<text x="${D}" y="${(c+u)/2}" text-anchor="${G}" dy="0.35em" font-size="11" font-weight="600" fill="${z}">h = ${Y}</text>`}else g+=`<text x="${s}" y="${u+C+14}" text-anchor="middle" font-size="11" font-weight="600" fill="${z}">h = ${Y}</text>`}else if(F.paramType==="dihedral"){const s=Math.max(...y.map(c=>b(c.v)))+C+12;g+=`<text x="${m}" y="${s}" text-anchor="middle" font-size="11" font-weight="600" fill="${z}">φ = ${Y}</text>`}return o&&(g+=`<text x="${n/2}" y="14" text-anchor="middle" font-size="11" font-weight="700" fill="${z}">${o}</text>`),g}function at(t){const n=E.querySelector("#mol-vis");if(!n||t.length===0)return;const r=520,o=200;if(t.length===1){const l=t[0],e=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${r} ${o}" style="width:100%;max-width:${r}px;display:block;margin:0 auto;">${kt(l.param,r,o,l.label,l.color)}</svg>`;n.innerHTML=e}else{const l=Math.floor(r/t.length);let e='<div style="display:flex;gap:4px;justify-content:center;">';for(const a of t)e+=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${l} ${o}" style="flex:1;max-width:${l}px;">${kt(a.param,l,o,a.label,a.color)}</svg>`;e+="</div>",n.innerHTML=e}}function nt(t,n){const r=E.querySelector("#graph-container");if(!r||t.length<1)return;const o=St(),e={bond:H("opt.xBond"),angle:H("opt.xAngle"),height:H("opt.xHeight"),dihedral:H("opt.xDihedral")}[F.paramType],a=ct(),p=520,i=360,d=72,f=24,x=36,T=44,m=p-d-f,v=i-x-T,$=[];for(const s of a)for(const c of t){const u=c[s.key];u!=null&&isFinite(u)&&$.push(u)}$.length===0&&$.push(0);const b=Math.min(...$),C=Math.max(...$),y=(C-b)*.1||.01,g=b-y,N=C+y,L=s=>d+(s-R)/(X-R||1)*m,B=s=>x+v-(s-g)/(N-g||1)*v;let S=`<svg width="${p}" height="${i}" xmlns="http://www.w3.org/2000/svg" style="display:block;width:100%;height:auto;max-width:${p}px;" viewBox="0 0 ${p} ${i}">`;S+=`<rect x="${d}" y="${x}" width="${m}" height="${v}" fill="${o.surface}" rx="2"/>`;for(let s=0;s<=5;s++){const c=g+(N-g)*s/5,u=B(c);S+=`<line x1="${d}" y1="${u}" x2="${d+m}" y2="${u}" stroke="${o.grid}" stroke-width="0.5"/>`,S+=`<text x="${d-6}" y="${u+3}" text-anchor="end" font-size="9" font-family="monospace" fill="${o.dim}">${c.toFixed(4)}</text>`}for(let s=0;s<=5;s++){const c=R+(X-R)*s/5,u=L(c);S+=`<line x1="${u}" y1="${x}" x2="${u}" y2="${x+v}" stroke="${o.grid}" stroke-width="0.5"/>`,S+=`<text x="${u}" y="${x+v+14}" text-anchor="middle" font-size="9" fill="${o.dim}">${c.toFixed(2)}</text>`}S+=`<line x1="${d}" y1="${x}" x2="${d}" y2="${x+v}" stroke="${o.axis}" stroke-width="1"/>`,S+=`<line x1="${d}" y1="${x+v}" x2="${d+m}" y2="${x+v}" stroke="${o.axis}" stroke-width="1"/>`;for(const s of a){const c=t.filter(h=>h[s.key]!=null&&isFinite(h[s.key])).map(h=>({x:h.param,y:h[s.key]}));if(c.length<2){for(const h of c)S+=`<circle cx="${L(h.x).toFixed(1)}" cy="${B(h.y).toFixed(1)}" r="3" fill="${s.color}"/>`;continue}let u="";for(let h=0;h<c.length;h++){const w=L(c[h].x),M=B(c[h].y);u+=h===0?`M${w.toFixed(1)},${M.toFixed(1)}`:` L${w.toFixed(1)},${M.toFixed(1)}`}S+=`<path d="${u}" fill="none" stroke="${s.color}" stroke-width="2" stroke-linejoin="round"/>`;for(const h of c)S+=`<circle cx="${L(h.x).toFixed(1)}" cy="${B(h.y).toFixed(1)}" r="3" fill="${s.color}" stroke="${o.surface}" stroke-width="0.8"/>`;if(c.length>=3&&n){let h=c[0];for(const I of c)I.y<h.y&&(h=I);const w=L(h.x),M=B(h.y);S+=`<line x1="${w}" y1="${M}" x2="${w}" y2="${x+v}" stroke="${s.color}" stroke-width="1" stroke-dasharray="3,3" opacity="0.5"/>`,S+=ae(w,M,6,s.color)}}if(n&&K>=0&&K<t.length){const s=t[K],c=L(s.param),u=it()?"#ffd700":"#cc8800";S+=`<line x1="${c}" y1="${x}" x2="${c}" y2="${x+v}" stroke="${u}" stroke-width="1" stroke-dasharray="3,2" opacity="0.6"/>`;for(const h of a){const w=s[h.key];w!=null&&isFinite(w)&&(S+=`<circle cx="${c}" cy="${B(w)}" r="5" fill="${h.color}" stroke="${u}" stroke-width="2"/>`)}}const Y=d+8;let z=x+14;for(const s of a)S+=`<rect x="${Y}" y="${z-7}" width="10" height="3" rx="1" fill="${s.color}"/>`,S+=`<text x="${Y+14}" y="${z-3}" font-size="9" font-weight="600" fill="${o.dim}">${s.label}</text>`,z+=14;S+=`<text x="${d+m/2}" y="20" text-anchor="middle" font-size="12" font-weight="600" fill="${o.titleSvg}">${H("opt.graphTitle")}</text>`,S+=`<text x="${d+m/2}" y="${i-4}" text-anchor="middle" font-size="10" fill="${o.dim}">${e}</text>`,S+=`<text x="14" y="${x+v/2}" text-anchor="middle" font-size="10" fill="${o.dim}" transform="rotate(-90,14,${x+v/2})">${H("opt.yEnergy")}</text>`,S+="</svg>",r.innerHTML=S}function ae(t,n,r,o){const l=[];for(let e=0;e<10;e++){const a=Math.PI/2+e*Math.PI/5,p=e%2===0?r:r*.4;l.push(`${(t+p*Math.cos(a)).toFixed(1)},${(n-p*Math.sin(a)).toFixed(1)}`)}return`<polygon points="${l.join(" ")}" fill="${o}" stroke="${o}" stroke-width="0.5"/>`}let Ct=!1;function ne(){if(Ct)return;Ct=!0;const t=document.createElement("style");t.textContent=`
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
  `,document.head.appendChild(t)}Ot();Xt();q();
