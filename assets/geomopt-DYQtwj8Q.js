import"./styles-B8WF9G-F.js";import{p as et}from"./parseXYZ-BE3rbSP_.js";import{z as xt,g as ut,f as ht,b as yt,n as bt}from"./ri-Cuj4t-H2.js";import{a as vt,H as $t,b as wt,c as Ft}from"./theoryControls-AWp8GxPQ.js";import{c as Mt}from"./gradient-gHuNZ-_b.js";import{c as St,i as zt,r as Ht,t as $,a as kt,b as At,g as Et}from"./nav-BKPoqOVV.js";import{c as Nt}from"./moleculeViewer3D-DIiPA6tx.js";const k=1.8897259886;function T(c){return`${c.length}
geomopt
${c.map(([a,s,r,e])=>`${a}  ${s.toFixed(6)}  ${r.toFixed(6)}  ${e.toFixed(6)}`).join(`
`)}`}const dt=[{id:"h2",labelKey:"gopt.scenH2",descKey:"gopt.descH2",charge:0,mult:1,basisName:"sto-3g",xyz:T([["H",0,0,0],["H",0,0,2*k]])},{id:"hf",labelKey:"gopt.scenHF",descKey:"gopt.descHF",charge:0,mult:1,basisName:"sto-3g",xyz:T([["H",0,0,0],["F",0,0,2.2*k]])},{id:"lih",labelKey:"gopt.scenLiH",descKey:"gopt.descLiH",charge:0,mult:1,basisName:"sto-3g",xyz:T([["Li",0,0,0],["H",0,0,2*k]])},{id:"h2o",labelKey:"gopt.scenH2O",descKey:"gopt.descH2O",charge:0,mult:1,basisName:"sto-3g",xyz:T([["O",0,0,0],["H",1.1*k,0,0],["H",0,1.1*k,0]])},{id:"nh3",labelKey:"gopt.scenNH3",descKey:"gopt.descNH3",charge:0,mult:1,basisName:"sto-3g",xyz:(()=>{const c=1.05*k,a=95*Math.PI/180,s=c*Math.cos(a),r=c*Math.sin(a);return T([["N",0,0,0],["H",r,0,s],["H",-r*.5,r*Math.sqrt(3)/2,s],["H",-r*.5,-r*Math.sqrt(3)/2,s]])})()},{id:"beh2",labelKey:"gopt.scenBeH2",descKey:"gopt.descBeH2",charge:0,mult:1,basisName:"sto-3g",xyz:T([["Be",0,0,0],["H",1.4*k,0,1*k],["H",-1.4*k,0,1*k]])},{id:"hcn",labelKey:"gopt.scenHCN",descKey:"gopt.descHCN",charge:0,mult:1,basisName:"sto-3g",xyz:T([["H",0,0,0],["C",0,0,1.2*k],["N",0,0,1.2*k+1.3*k]])},{id:"ch4",labelKey:"gopt.scenCH4",descKey:"gopt.descCH4",charge:0,mult:1,basisName:"sto-3g",xyz:(()=>{const c=1.15*k,a=100*Math.PI/180,s=Math.sin(a),r=Math.cos(a);return T([["C",0,0,0],["H",c*s,0,c*r],["H",-c*s*.5,c*s*Math.sqrt(3)/2,c*r],["H",-c*s*.5,-c*s*Math.sqrt(3)/2,c*r],["H",0,0,-c]])})()},{id:"bh3",labelKey:"gopt.scenBH3",descKey:"gopt.descBH3",charge:0,mult:1,basisName:"sto-3g",xyz:(()=>{const c=1.25*k;return T([["B",0,0,0],["H",c,0,.5*k],["H",-c*.5,c*Math.sqrt(3)/2,.5*k],["H",-c*.5,-c*Math.sqrt(3)/2,.5*k]])})()}],pt=[{label:"Steepest Descent",items:[{id:"sd",name:"SD"}]},{label:"Conjugate Gradient",items:[{id:"cg-fr",name:"CG-FR"},{id:"cg-pr",name:"CG-PR"},{id:"cg-hs",name:"CG-HS"},{id:"cg-dy",name:"CG-DY"}]},{label:"Quasi-Newton",items:[{id:"bfgs",name:"BFGS"},{id:"dfp",name:"DFP"},{id:"sr1",name:"SR1"}]},{label:"GDIIS",items:[{id:"gdiis",name:"GDIIS"}]}];let P=dt[0],q="bfgs",b=[],I=!1,tt=!1,K=!1,_=!1,gt=0,W=50,C=0,Q="HF";const V=45e-5;let O=null,it="auto";const st=new Map;async function qt(c){const a=st.get(c);if(a)return a;const s=`/GANSU-Lite/basis/${c.toLowerCase()}.gbs`,r=ht.fromGBS(await(await fetch(s)).text());return st.set(c,r),r}function lt(c,a,s){const r=3*s,e=[];for(let t=0;t<3;t++){const i=new Float64Array(r);for(let n=0;n<s;n++)i[3*n+t]=1;e.push(i)}let g=0,x=0,l=0;for(let t=0;t<s;t++)g+=a[3*t],x+=a[3*t+1],l+=a[3*t+2];g/=s,x/=s,l/=s;const o=[];for(let t=0;t<3;t++){const i=new Float64Array(r);for(let n=0;n<s;n++){const d=a[3*n]-g,u=a[3*n+1]-x,w=a[3*n+2]-l;t===0?(i[3*n+1]=-w,i[3*n+2]=u):t===1?(i[3*n]=w,i[3*n+2]=-d):(i[3*n]=-u,i[3*n+1]=d)}o.push(i)}const m=[];for(const t of[...e,...o]){const i=new Float64Array(t);for(const d of m){let u=0;for(let w=0;w<r;w++)u+=i[w]*d[w];for(let w=0;w<r;w++)i[w]-=u*d[w]}let n=0;for(let d=0;d<r;d++)n+=i[d]*i[d];if(n=Math.sqrt(n),n>1e-10){for(let d=0;d<r;d++)i[d]/=n;m.push(i)}}for(const t of m){let i=0;for(let n=0;n<r;n++)i+=c[n]*t[n];for(let n=0;n<r;n++)c[n]-=i*t[n]}}function Ct(c,a,s,r){const e=a.length;let g=0,x=0,l=0,o=0,m=0;for(let t=0;t<e;t++){const i=a[t]-s[t];g+=a[t]*a[t],x+=s[t]*s[t],l+=a[t]*i,o+=i*r[t],m+=a[t]*r[t]}if(x<1e-30)return 0;switch(c){case"cg-fr":return g/x;case"cg-pr":return Math.max(0,l/x);case"cg-hs":return o>1e-30?l/o:0;case"cg-dy":return o>1e-30?g/o:0;default:return 0}}function Lt(c,a,s,r,e){let g=0;for(let t=0;t<e;t++)g+=s[t]*r[t];if(c==="sr1"){const t=new Float64Array(e);for(let d=0;d<e;d++){let u=s[d];for(let w=0;w<e;w++)u-=a[d*e+w]*r[w];t[d]=u}let i=0;for(let d=0;d<e;d++)i+=t[d]*r[d];if(Math.abs(i)<1e-12*Math.sqrt(g>0?g:1))return a;const n=new Float64Array(a);for(let d=0;d<e;d++)for(let u=0;u<e;u++)n[d*e+u]+=t[d]*t[u]/i;return n}if(g<1e-14)return a;const x=1/g,l=new Float64Array(e);for(let t=0;t<e;t++){let i=0;for(let n=0;n<e;n++)i+=a[t*e+n]*r[n];l[t]=i}let o=0;for(let t=0;t<e;t++)o+=r[t]*l[t];const m=new Float64Array(e*e);if(c==="bfgs")for(let t=0;t<e;t++)for(let i=0;i<e;i++)m[t*e+i]=a[t*e+i]-x*(s[t]*l[i]+l[t]*s[i])+x*(x*o+1)*s[t]*s[i];else for(let t=0;t<e;t++)for(let i=0;i<e;i++)m[t*e+i]=a[t*e+i]-l[t]*l[i]/(o||1)+x*s[t]*s[i];return m}const j=[],Tt=6;function It(c){const a=j.length;if(a<2)return null;const s=new Float64Array((a+1)*(a+1));for(let o=0;o<a;o++)for(let m=0;m<a;m++){let t=0;for(let i=0;i<c;i++)t+=j[o].grad[i]*j[m].grad[i];s[o*(a+1)+m]=t}for(let o=0;o<a;o++)s[a*(a+1)+o]=-1,s[o*(a+1)+a]=-1;s[a*(a+1)+a]=0;const r=new Float64Array(a+1);r[a]=-1;const e=a+1,g=new Float64Array(e*(e+1));for(let o=0;o<e;o++){for(let m=0;m<e;m++)g[o*(e+1)+m]=s[o*e+m];g[o*(e+1)+e]=r[o]}for(let o=0;o<e;o++){let m=o,t=Math.abs(g[o*(e+1)+o]);for(let n=o+1;n<e;n++){const d=Math.abs(g[n*(e+1)+o]);d>t&&(t=d,m=n)}if(t<1e-14)return null;if(m!==o)for(let n=0;n<=e;n++){const d=g[o*(e+1)+n];g[o*(e+1)+n]=g[m*(e+1)+n],g[m*(e+1)+n]=d}const i=g[o*(e+1)+o];for(let n=o;n<=e;n++)g[o*(e+1)+n]/=i;for(let n=0;n<e;n++){if(n===o)continue;const d=g[n*(e+1)+o];for(let u=o;u<=e;u++)g[n*(e+1)+u]-=d*g[o*(e+1)+u]}}const x=new Float64Array(a);for(let o=0;o<a;o++)x[o]=g[o*(e+1)+e];const l=new Float64Array(c);for(let o=0;o<a;o++)for(let m=0;m<c;m++)l[m]+=x[o]*j[o].coords[m];return l}async function Gt(){I=!0,tt=!1,b=[],K=!1,_=!1,C=0,j.length=0,B();const c=P,a=performance.now();await xt("/GANSU-Lite/"),it=bt()!=="js"?"wasm":"auto";const e=await qt(c.basisName),g=et(c.xyz),x=g.length,l=3*x;let o=new Float64Array(l);for(let S=0;S<x;S++)o[3*S]=g[S].coordinate.x,o[3*S+1]=g[S].coordinate.y,o[3*S+2]=g[S].coordinate.z;const m=g.map(S=>S.atomicNumber);let t=.3,i=null,n=null,d=null,u=!1;const w=q.startsWith("cg-"),G=q==="bfgs"||q==="dfp"||q==="sr1",U=q==="gdiis";for(let S=0;S<W&&!tt;S++){const E=m.map((h,p)=>({atomicNumber:h,coordinate:{x:o[3*p],y:o[3*p+1],z:o[3*p+2]},atomIndex:p})),v=new ut(E,e,c.charge),F=await wt(v,e,Q,"RHF"),L=await F.solve({eriBackend:it}),R=v.numAlphaSpins,ft=F.xcFunctional&&F.grid?{functional:F.xcFunctional,grid:F.grid}:void 0,M=Mt(v.primitiveShells,v.atoms,v.cgtoNormalizationFactors,v.numBasis,R,F.density,F.coefficients,F.orbitalEnergies,void 0,ft).total;lt(M,o,x);let D=0,at=0;for(let h=0;h<l;h++)D=Math.max(D,Math.abs(M[h])),at+=M[h]*M[h];if(b.push({energy:L,maxForce:D,rmsForce:Math.sqrt(at/l),coords:Array.from(o),gradient:Array.from(M)}),C=b.length-1,Bt(S+1,D),ot(),mt(),rt(),D<V){_=!0;break}if(b.length>=2&&L>b[b.length-2].energy&&(!G||!u)){const h=b[b.length-2];for(let p=0;p<l;p++)o[p]=h.coords[p];t*=.5,n=null,await new Promise(p=>setTimeout(p,0));continue}let H=new Float64Array(l);if(q==="sd")for(let h=0;h<l;h++)H[h]=-M[h];else if(w){let h=0;i&&n&&(h=Ct(q,M,i,n));for(let f=0;f<l;f++)H[f]=-M[f]+h*(n?n[f]:0);let p=0;for(let f=0;f<l;f++)p+=H[f]*M[f];if(p>0)for(let f=0;f<l;f++)H[f]=-M[f]}else if(G){if(!d){d=new Float64Array(l*l);for(let p=0;p<l;p++)d[p*l+p]=1}if(i&&b.length>=2){const p=b[b.length-2].coords,f=new Float64Array(l),y=new Float64Array(l);for(let z=0;z<l;z++)f[z]=o[z]-p[z],y[z]=M[z]-i[z];lt(f,o,x);let N=0,J=0,nt=0;for(let z=0;z<l;z++)N+=f[z]*y[z],J+=y[z]*y[z],nt+=f[z]*f[z];if(N>1e-14&&nt>1e-20){if(!u&&J>1e-14){const z=N/J;for(let X=0;X<l;X++)d[X*l+X]=z;u=!0}d=Lt(q,d,f,y,l)}}for(let p=0;p<l;p++){let f=0;for(let y=0;y<l;y++)f+=d[p*l+y]*-M[y];H[p]=f}let h=0;for(let p=0;p<l;p++)h+=H[p]*M[p];if(h>0){d=new Float64Array(l*l);for(let p=0;p<l;p++)d[p*l+p]=1;u=!1;for(let p=0;p<l;p++)H[p]=-M[p]}}else if(U){j.push({coords:Float64Array.from(o),grad:Float64Array.from(M)}),j.length>Tt&&j.shift();const h=It(l);if(h){for(let y=0;y<l;y++)H[y]=h[y]-o[y];const p=Math.sqrt(H.reduce((y,N)=>y+N*N,0));let f=0;for(let y=0;y<l;y++)f+=H[y]*M[y];if(p<1e-12||f>0)for(let y=0;y<l;y++)H[y]=-M[y];else{const y=new Float64Array(l);for(let N=0;N<l;N++)y[N]=o[N]+H[N];o=y,i=Float64Array.from(M),n=H,await new Promise(N=>setTimeout(N,0));continue}}else for(let p=0;p<l;p++)H[p]=-M[p]}const Y=.5,Z=Math.sqrt(H.reduce((h,p)=>h+p*p,0));if(G&&u){const h=Z>Y?Y/Z:1,p=new Float64Array(l);for(let f=0;f<l;f++)p[f]=o[f]+h*H[f];o=p}else{const h=Math.min(t,Y)/(Z||1),p=new Float64Array(l);for(let f=0;f<l;f++)p[f]=o[f]+h*H[f];o=p}(!G||!u)&&(t=Math.min(t*1.2,Y)),i=Float64Array.from(M),n=H,await new Promise(h=>setTimeout(h,0))}gt=performance.now()-a,I=!1,K=!0,B()}let ct=!1;function jt(){if(ct)return;ct=!0;const c=document.createElement("style");c.textContent=`
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body { background: var(--color-bg); color: var(--color-text);
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; min-height: 100vh; }
    .opt-page { max-width: 980px; margin: 0 auto; padding: 16px 20px; }
    .opt-content { display: flex; gap: 20px; }
    .opt-panel { background: var(--color-surface); border: 1px solid var(--color-border); border-radius: 10px; padding: 16px 18px; }
    .opt-controls { flex: 0 0 280px; }
    .opt-graph-panel { flex: 1; min-width: 0; display: flex; flex-direction: column; align-items: center; padding-top: 8px; }
    .opt-controls h2 {
      font-size: 0.72rem; font-weight: 600; color: var(--color-text-secondary);
      text-transform: uppercase; letter-spacing: 0.04em; margin: 10px 0 5px; }
    .opt-controls h2:first-child { margin-top: 0; }
    .opt-scenario-grid { display: flex; flex-direction: column; gap: 5px; }
    .opt-scenario-card {
      padding: 7px 12px; border: 2px solid var(--color-border); border-radius: 8px;
      cursor: pointer; transition: all 0.15s; display: flex; align-items: baseline; gap: 8px; }
    .opt-scenario-card:hover { background: var(--color-surface-alt); border-color: var(--color-text-dim); }
    .opt-scenario-card.selected { border-color: var(--color-accent); background: var(--color-surface-alt); box-shadow: 0 0 0 1px var(--color-accent); }
    .opt-scenario-card strong { font-size: 0.85rem; white-space: nowrap; }
    .conv-desc { font-size: 0.7rem; color: var(--color-text-dim); }

    .gopt-algo-group { margin-bottom: 4px; }
    .gopt-algo-label { font-size: 0.65rem; color: var(--color-text-dim); margin-bottom: 2px; }
    .gopt-algo-row { display: flex; gap: 3px; margin-bottom: 3px; }
    .gopt-algo-btn {
      flex: 1; padding: 4px 2px; border: 1px solid var(--color-border); border-radius: 5px;
      background: none; color: var(--color-text); font-size: 0.7rem; cursor: pointer;
      transition: all 0.12s; text-align: center; }
    .gopt-algo-btn:hover { background: var(--color-surface-alt); }
    .gopt-algo-btn.active { border-color: var(--color-accent); background: var(--color-surface-alt); font-weight: 600; box-shadow: 0 0 0 1px var(--color-accent); }

    .opt-run-btn {
      width: 100%; margin-top: 10px; padding: 10px; border: none; border-radius: 8px;
      font-size: 0.85rem; font-weight: 600; cursor: pointer;
      background: var(--color-accent); color: var(--color-accent-on); transition: background 0.15s; }
    .opt-run-btn:hover:not([disabled]) { background: var(--color-accent-hover); }
    .opt-run-btn[disabled] { opacity: 0.6; cursor: not-allowed; }
    .opt-stop-btn {
      width: 100%; margin-top: 6px; padding: 8px; border: 1px solid var(--color-error); border-radius: 8px;
      font-size: 0.8rem; font-weight: 600; cursor: pointer; background: none; color: var(--color-error); }

    .opt-progress { height: 6px; background: var(--color-progress-bg, #e0e4ea); border-radius: 3px; margin-top: 10px; overflow: hidden; }
    .opt-progress-bar { height: 100%; background: var(--color-accent); border-radius: 3px; transition: width 0.2s ease; }
    .opt-progress-text { font-size: 0.72rem; color: var(--color-text-dim); text-align: center; margin-top: 4px; }

    .opt-summary { margin-top: 10px; padding: 10px; background: var(--color-surface-alt); border-radius: 8px; border: 1px solid var(--color-border); }
    .opt-summary h3 { font-size: 0.8rem; margin-bottom: 6px; }
    .opt-summary table { width: 100%; font-size: 0.72rem; border-collapse: collapse; }
    .opt-summary td { padding: 2px 4px; }

    #graph-container { width: 100%; text-align: center; }
    #mol-vis { width: 100%; max-width: 400px; margin: 6px auto 0; }
    .opt-hint { color: var(--color-text-dim); font-size: 0.85rem; padding: 60px 20px; }
    .gopt-slider { width: 100%; max-width: 500px; margin: 6px auto 0; display: flex; align-items: center; gap: 10px; }
    .gopt-slider input[type=range] { flex: 1; cursor: pointer; }
    .gopt-slider-label { font-size: 0.78rem; font-weight: 600; min-width: 90px; text-align: center; }
    .gopt-step-info { font-size: 0.72rem; color: var(--color-text-dim); text-align: center; margin-top: 2px; }
    .gopt-xyz-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; margin-top: 10px; }
    .gopt-xyz-block { background: var(--color-surface-alt); border: 1px solid var(--color-border); border-radius: 6px; padding: 6px 8px; overflow: hidden; }
    .gopt-xyz-label { font-size: 0.65rem; color: var(--color-text-dim); text-transform: uppercase; letter-spacing: 0.04em; margin-bottom: 3px; font-weight: 600; }
    .gopt-xyz-pre { font-family: ui-monospace, 'SF Mono', Consolas, monospace; font-size: 0.65rem; line-height: 1.35; color: var(--color-text); white-space: pre; overflow-x: auto; margin: 0; }
    @media (max-width: 700px) { .gopt-xyz-grid { grid-template-columns: 1fr; } }
    @media (max-width: 700px) { .opt-content { flex-direction: column; } .opt-controls { flex: none; } }
  `,document.head.appendChild(c)}const A=document.getElementById("app");function B(){O=null,A.innerHTML=`
    <div class="opt-page">
      ${Ht("geomopt")}
      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${$("gopt.molecule")}</h2>
          <div class="opt-scenario-grid" id="scen-grid"></div>

          <h2>${$("gopt.algorithm")}</h2>
          <div id="algo-area"></div>

          <div class="theory-row" style="margin:8px 0;display:flex;align-items:center;gap:8px;font-size:0.9rem;">
            <span>Theory:</span>${vt("theory-sel",Q,"",$t)}
          </div>

          <div style="margin:6px 0 0;">
            <label style="font-size:0.72rem;">${$("gopt.maxIter")}
              <input type="number" id="max-iter" min="5" max="200" value="${W}" style="width:50px;margin-left:4px;">
            </label>
          </div>

          <button id="run-btn" class="opt-run-btn" ${I?"disabled":""}>
            ${I?$("gopt.running"):$("gopt.run")}
          </button>
          ${I?`<button id="stop-btn" class="opt-stop-btn">${$("gopt.stop")}</button>`:""}
          <div id="progress-area"></div>
        </div>
        <div class="opt-panel opt-graph-panel">
          <div id="graph-container">
            ${!K&&!I&&b.length===0?`<p class="opt-hint">${$("gopt.waiting")}</p>`:""}
          </div>
          <div id="slider-area"></div>
          <div id="mol-vis"></div>
          <div id="summary-area"></div>
        </div>
      </div>
    </div>`,jt();const c=A.querySelector("#scen-grid");for(const r of dt){const e=document.createElement("div");e.className="opt-scenario-card"+(r.id===P.id?" selected":""),e.innerHTML=`<strong>${$(r.labelKey)}</strong><span class="conv-desc">${$(r.descKey)}</span>`,e.addEventListener("click",()=>{I||(P=r,b=[],K=!1,_=!1,C=0,O=null,B())}),c.appendChild(e)}const a=A.querySelector("#algo-area");for(const r of pt){const e=document.createElement("div");e.className="gopt-algo-group",e.innerHTML=`<div class="gopt-algo-label">${r.label}</div><div class="gopt-algo-row">${r.items.map(g=>`<button class="gopt-algo-btn${g.id===q?" active":""}" data-algo="${g.id}">${g.name}</button>`).join("")}</div>`,a.appendChild(e)}for(const r of A.querySelectorAll(".gopt-algo-btn"))r.addEventListener("click",()=>{if(!I){q=r.dataset.algo;for(const e of A.querySelectorAll(".gopt-algo-btn"))e.classList.remove("active");r.classList.add("active")}});A.querySelector("#nav-theme").addEventListener("click",()=>{kt(),B()}),A.querySelector("#nav-lang").addEventListener("click",()=>{At(),B()});const s=A.querySelector("#theory-sel");s&&s.addEventListener("change",()=>{Q=s.value}),A.querySelector("#run-btn").addEventListener("click",()=>{I||(W=+(A.querySelector("#max-iter")?.value??50),Gt())}),A.querySelector("#stop-btn")?.addEventListener("click",()=>{tt=!0}),b.length>0&&(ot(),Kt(),mt(),rt())}function Bt(c,a,s){const r=A.querySelector("#progress-area");if(!r)return;const e=`Step ${c} | Max|F| = ${a.toExponential(2)}`;r.innerHTML=`
    <div class="opt-progress"><div class="opt-progress-bar" style="width:${Math.min(100,c/W*100).toFixed(0)}%"></div></div>
    <div class="opt-progress-text">${e}</div>`}function Kt(){const c=A.querySelector("#summary-area");if(!c||b.length===0)return;const a=b[0],s=b[b.length-1],r=_?"var(--color-accent)":"var(--color-error)",e=_?`✅ ${$("gopt.converged")}`:K?`⚠️ ${$("gopt.notConverged")}`:$("gopt.running"),g=pt.flatMap(t=>t.items).find(t=>t.id===q),x=(t,i)=>`<tr><td style="color:var(--color-text-dim)">${t}</td><td>${i}</td></tr>`,l=.529177210903,o=et(P.xyz).map(t=>t.atomicNumber),m=t=>{const i=[];for(let n=0;n<o.length;n++){const d=yt(o[n]),u=(t[3*n]*l).toFixed(6).padStart(11),w=(t[3*n+1]*l).toFixed(6).padStart(11),G=(t[3*n+2]*l).toFixed(6).padStart(11);i.push(`${d.padEnd(2)}${u}${w}${G}`)}return i.join(`
`)};c.innerHTML=`<div class="opt-summary"><h3 style="color:${r}">${e}</h3><table>
    ${x("Theory",Ft(Q))}
    ${x($("gopt.algorithm"),g?.name??q)}
    ${x($("gopt.steps"),String(b.length))}
    ${x($("gopt.initEnergy"),a.energy.toFixed(8)+" Eh")}
    ${x($("gopt.finalEnergy"),s.energy.toFixed(8)+" Eh")}
    ${x($("gopt.energyChange"),(s.energy-a.energy).toFixed(8)+" Eh")}
    ${x($("gopt.maxForce"),s.maxForce.toExponential(4)+" Eh/bohr")}
    ${x($("gopt.threshold"),V.toExponential(1)+" Eh/bohr")}
    ${K?x($("gopt.time"),(gt/1e3).toFixed(1)+" s"):""}
  </table>
  <div class="gopt-xyz-grid">
    <div class="gopt-xyz-block">
      <div class="gopt-xyz-label">Initial (Å)</div>
      <pre class="gopt-xyz-pre">${m(Array.from(a.coords))}</pre>
    </div>
    <div class="gopt-xyz-block">
      <div class="gopt-xyz-label">Final (Å)</div>
      <pre class="gopt-xyz-pre">${m(Array.from(s.coords))}</pre>
    </div>
  </div>
  </div>`}function mt(){const c=A.querySelector("#slider-area");if(!c||b.length<2)return;const a=b[C];c.innerHTML=`
    <div class="gopt-slider">
      <span class="gopt-slider-label">Step ${C+1}/${b.length}</span>
      <input type="range" id="step-slider" min="0" max="${b.length-1}" value="${C}" step="1">
    </div>
    <div class="gopt-step-info">E = ${a.energy.toFixed(8)} Eh &nbsp; Max|F| = ${a.maxForce.toExponential(3)}</div>`,A.querySelector("#step-slider").addEventListener("input",s=>{C=+s.target.value,ot(),rt();const r=b[C];c.querySelector(".gopt-slider-label").textContent=`Step ${C+1}/${b.length}`,c.querySelector(".gopt-step-info").innerHTML=`E = ${r.energy.toFixed(8)} Eh &nbsp; Max|F| = ${r.maxForce.toExponential(3)}`})}function ot(){const c=A.querySelector("#graph-container");if(!c||b.length<1)return;const a=500,s=155,r={l:70,r:20,t:24,b:28},e=Et();function g(x,l,o,m){const t=o.length,i=1,n=Math.max(t,2);let d,u;if(m){const v=o.map(F=>Math.log10(Math.max(F,1e-15)));d=Math.min(...v),u=Math.max(...v),u-d<.5&&(d-=.25,u+=.25)}else{d=Math.min(...o),u=Math.max(...o);const v=(u-d)*.1||.001;d-=v,u+=v}const w=a-r.l-r.r,G=s-r.t-r.b,U=v=>r.l+(v-i)/(n-i)*w,S=v=>r.t+(1-(v-d)/(u-d))*G;let E=`<svg viewBox="0 0 ${a} ${s}" style="width:100%;max-width:${a}px;background:${e.surface};border-radius:6px;margin-bottom:4px;">`;E+=`<text x="${a/2}" y="14" text-anchor="middle" font-size="10" fill="${e.titleSvg}" font-weight="600">${x}</text>`,E+=`<line x1="${r.l}" y1="${r.t}" x2="${r.l}" y2="${s-r.b}" stroke="${e.grid}"/>`,E+=`<line x1="${r.l}" y1="${s-r.b}" x2="${a-r.r}" y2="${s-r.b}" stroke="${e.grid}"/>`,E+=`<text x="12" y="${s/2}" text-anchor="middle" font-size="7" fill="${e.dim}" transform="rotate(-90,12,${s/2})">${l}</text>`,E+=`<text x="${r.l+w/2}" y="${s-3}" text-anchor="middle" font-size="7" fill="${e.dim}">${$("gopt.xStep")}</text>`;for(let v=0;v<=4;v++){const F=d+(u-d)*v/4,L=S(F);E+=`<line x1="${r.l-3}" y1="${L}" x2="${a-r.r}" y2="${L}" stroke="${e.grid}" stroke-dasharray="2,3"/>`,E+=`<text x="${r.l-5}" y="${L+3}" text-anchor="end" font-size="6.5" fill="${e.dim}">${m?`1e${F.toFixed(1)}`:F.toPrecision(6)}</text>`}if(m){const v=S(Math.log10(V));v>=r.t&&v<=s-r.b&&(E+=`<line x1="${r.l}" y1="${v}" x2="${a-r.r}" y2="${v}" stroke="#ef4444" stroke-dasharray="4,3"/><text x="${a-r.r-2}" y="${v-3}" text-anchor="end" font-size="7" fill="#ef4444">threshold</text>`)}if(t>=1){const v=o.map((F,L)=>`${U(L+1).toFixed(1)},${S(m?Math.log10(Math.max(F,1e-15)):F).toFixed(1)}`);E+=`<polyline points="${v.join(" ")}" fill="none" stroke="${e.occupied}" stroke-width="1.5"/>`;for(let F=0;F<t;F++){const L=m?Math.log10(Math.max(o[F],1e-15)):o[F],R=F===C;E+=`<circle cx="${U(F+1).toFixed(1)}" cy="${S(L).toFixed(1)}" r="${R?5:2.5}" fill="${R?e.accent:e.occupied}" ${R?`stroke="${e.titleSvg}" stroke-width="1.5"`:""} />`}}return E+="</svg>",E}c.innerHTML=g($("gopt.graphEnergy"),$("gopt.yEnergy"),b.map(x=>x.energy),!1)+g($("gopt.graphGrad"),$("gopt.yForce"),b.map(x=>x.maxForce),!0)}function rt(){const c=A.querySelector("#mol-vis");if(!c||b.length===0)return;const a=b[C],s=et(P.xyz).map(r=>r.atomicNumber);O?O.update(a.coords,a.gradient,V):O=Nt(c,s,a.coords,a.gradient,V)}St();zt();B();
