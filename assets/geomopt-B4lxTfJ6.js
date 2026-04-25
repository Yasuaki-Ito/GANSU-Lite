import{c as mt,i as ft,r as xt,t as $,a as ut,b as ht,g as yt}from"./nav-Lc1uTV-6.js";import{p as tt}from"./parseXYZ-DTmoiCYR.js";import{O as bt,t as vt,u as $t,d as wt,b as Ft,x as Mt}from"./builder-B-MUn5dT.js";import{c as zt}from"./gradient-CXo3J1Uj.js";import{c as St}from"./moleculeViewer3D-B1aGFjc3.js";const k=1.8897259886;function I(c){return`${c.length}
geomopt
${c.map(([a,i,r,o])=>`${a}  ${i.toFixed(6)}  ${r.toFixed(6)}  ${o.toFixed(6)}`).join(`
`)}`}const ct=[{id:"h2",labelKey:"gopt.scenH2",descKey:"gopt.descH2",charge:0,mult:1,basisName:"sto-3g",xyz:I([["H",0,0,0],["H",0,0,2*k]])},{id:"hf",labelKey:"gopt.scenHF",descKey:"gopt.descHF",charge:0,mult:1,basisName:"sto-3g",xyz:I([["H",0,0,0],["F",0,0,2.2*k]])},{id:"lih",labelKey:"gopt.scenLiH",descKey:"gopt.descLiH",charge:0,mult:1,basisName:"sto-3g",xyz:I([["Li",0,0,0],["H",0,0,2*k]])},{id:"h2o",labelKey:"gopt.scenH2O",descKey:"gopt.descH2O",charge:0,mult:1,basisName:"sto-3g",xyz:I([["O",0,0,0],["H",1.1*k,0,0],["H",0,1.1*k,0]])},{id:"nh3",labelKey:"gopt.scenNH3",descKey:"gopt.descNH3",charge:0,mult:1,basisName:"sto-3g",xyz:(()=>{const c=1.05*k,a=95*Math.PI/180,i=c*Math.cos(a),r=c*Math.sin(a);return I([["N",0,0,0],["H",r,0,i],["H",-r*.5,r*Math.sqrt(3)/2,i],["H",-r*.5,-r*Math.sqrt(3)/2,i]])})()},{id:"beh2",labelKey:"gopt.scenBeH2",descKey:"gopt.descBeH2",charge:0,mult:1,basisName:"sto-3g",xyz:I([["Be",0,0,0],["H",1.4*k,0,1*k],["H",-1.4*k,0,1*k]])},{id:"hcn",labelKey:"gopt.scenHCN",descKey:"gopt.descHCN",charge:0,mult:1,basisName:"sto-3g",xyz:I([["H",0,0,0],["C",0,0,1.2*k],["N",0,0,1.2*k+1.3*k]])},{id:"ch4",labelKey:"gopt.scenCH4",descKey:"gopt.descCH4",charge:0,mult:1,basisName:"sto-3g",xyz:(()=>{const c=1.15*k,a=100*Math.PI/180,i=Math.sin(a),r=Math.cos(a);return I([["C",0,0,0],["H",c*i,0,c*r],["H",-c*i*.5,c*i*Math.sqrt(3)/2,c*r],["H",-c*i*.5,-c*i*Math.sqrt(3)/2,c*r],["H",0,0,-c]])})()},{id:"bh3",labelKey:"gopt.scenBH3",descKey:"gopt.descBH3",charge:0,mult:1,basisName:"sto-3g",xyz:(()=>{const c=1.25*k;return I([["B",0,0,0],["H",c,0,.5*k],["H",-c*.5,c*Math.sqrt(3)/2,.5*k],["H",-c*.5,-c*Math.sqrt(3)/2,.5*k]])})()}],dt=[{label:"Steepest Descent",items:[{id:"sd",name:"SD"}]},{label:"Conjugate Gradient",items:[{id:"cg-fr",name:"CG-FR"},{id:"cg-pr",name:"CG-PR"},{id:"cg-hs",name:"CG-HS"},{id:"cg-dy",name:"CG-DY"}]},{label:"Quasi-Newton",items:[{id:"bfgs",name:"BFGS"},{id:"dfp",name:"DFP"},{id:"sr1",name:"SR1"}]},{label:"GDIIS",items:[{id:"gdiis",name:"GDIIS"}]}];let P=ct[0],q="bfgs",b=[],L=!1,J=!1,T=!1,_=!1,pt=0,X=50,C=0;const V=45e-5;let D=null,nt="auto";const it=new Map;async function Ht(c){const a=it.get(c);if(a)return a;const i=`/GANSU-Lite/basis/${c.toLowerCase()}.gbs`,r=wt.fromGBS(await(await fetch(i)).text());return it.set(c,r),r}function st(c,a,i){const r=3*i,o=[];for(let t=0;t<3;t++){const s=new Float64Array(r);for(let n=0;n<i;n++)s[3*n+t]=1;o.push(s)}let f=0,x=0,l=0;for(let t=0;t<i;t++)f+=a[3*t],x+=a[3*t+1],l+=a[3*t+2];f/=i,x/=i,l/=i;const e=[];for(let t=0;t<3;t++){const s=new Float64Array(r);for(let n=0;n<i;n++){const d=a[3*n]-f,u=a[3*n+1]-x,w=a[3*n+2]-l;t===0?(s[3*n+1]=-w,s[3*n+2]=u):t===1?(s[3*n]=w,s[3*n+2]=-d):(s[3*n]=-u,s[3*n+1]=d)}e.push(s)}const g=[];for(const t of[...o,...e]){const s=new Float64Array(t);for(const d of g){let u=0;for(let w=0;w<r;w++)u+=s[w]*d[w];for(let w=0;w<r;w++)s[w]-=u*d[w]}let n=0;for(let d=0;d<r;d++)n+=s[d]*s[d];if(n=Math.sqrt(n),n>1e-10){for(let d=0;d<r;d++)s[d]/=n;g.push(s)}}for(const t of g){let s=0;for(let n=0;n<r;n++)s+=c[n]*t[n];for(let n=0;n<r;n++)c[n]-=s*t[n]}}function kt(c,a,i,r){const o=a.length;let f=0,x=0,l=0,e=0,g=0;for(let t=0;t<o;t++){const s=a[t]-i[t];f+=a[t]*a[t],x+=i[t]*i[t],l+=a[t]*s,e+=s*r[t],g+=a[t]*r[t]}if(x<1e-30)return 0;switch(c){case"cg-fr":return f/x;case"cg-pr":return Math.max(0,l/x);case"cg-hs":return e>1e-30?l/e:0;case"cg-dy":return e>1e-30?f/e:0;default:return 0}}function At(c,a,i,r,o){let f=0;for(let t=0;t<o;t++)f+=i[t]*r[t];if(c==="sr1"){const t=new Float64Array(o);for(let d=0;d<o;d++){let u=i[d];for(let w=0;w<o;w++)u-=a[d*o+w]*r[w];t[d]=u}let s=0;for(let d=0;d<o;d++)s+=t[d]*r[d];if(Math.abs(s)<1e-12*Math.sqrt(f>0?f:1))return a;const n=new Float64Array(a);for(let d=0;d<o;d++)for(let u=0;u<o;u++)n[d*o+u]+=t[d]*t[u]/s;return n}if(f<1e-14)return a;const x=1/f,l=new Float64Array(o);for(let t=0;t<o;t++){let s=0;for(let n=0;n<o;n++)s+=a[t*o+n]*r[n];l[t]=s}let e=0;for(let t=0;t<o;t++)e+=r[t]*l[t];const g=new Float64Array(o*o);if(c==="bfgs")for(let t=0;t<o;t++)for(let s=0;s<o;s++)g[t*o+s]=a[t*o+s]-x*(i[t]*l[s]+l[t]*i[s])+x*(x*e+1)*i[t]*i[s];else for(let t=0;t<o;t++)for(let s=0;s<o;s++)g[t*o+s]=a[t*o+s]-l[t]*l[s]/(e||1)+x*i[t]*i[s];return g}const B=[],Nt=6;function Et(c){const a=B.length;if(a<2)return null;const i=new Float64Array((a+1)*(a+1));for(let e=0;e<a;e++)for(let g=0;g<a;g++){let t=0;for(let s=0;s<c;s++)t+=B[e].grad[s]*B[g].grad[s];i[e*(a+1)+g]=t}for(let e=0;e<a;e++)i[a*(a+1)+e]=-1,i[e*(a+1)+a]=-1;i[a*(a+1)+a]=0;const r=new Float64Array(a+1);r[a]=-1;const o=a+1,f=new Float64Array(o*(o+1));for(let e=0;e<o;e++){for(let g=0;g<o;g++)f[e*(o+1)+g]=i[e*o+g];f[e*(o+1)+o]=r[e]}for(let e=0;e<o;e++){let g=e,t=Math.abs(f[e*(o+1)+e]);for(let n=e+1;n<o;n++){const d=Math.abs(f[n*(o+1)+e]);d>t&&(t=d,g=n)}if(t<1e-14)return null;if(g!==e)for(let n=0;n<=o;n++){const d=f[e*(o+1)+n];f[e*(o+1)+n]=f[g*(o+1)+n],f[g*(o+1)+n]=d}const s=f[e*(o+1)+e];for(let n=e;n<=o;n++)f[e*(o+1)+n]/=s;for(let n=0;n<o;n++){if(n===e)continue;const d=f[n*(o+1)+e];for(let u=e;u<=o;u++)f[n*(o+1)+u]-=d*f[e*(o+1)+u]}}const x=new Float64Array(a);for(let e=0;e<a;e++)x[e]=f[e*(o+1)+o];const l=new Float64Array(c);for(let e=0;e<a;e++)for(let g=0;g<c;g++)l[g]+=x[e]*B[e].coords[g];return l}async function qt(){L=!0,J=!1,b=[],T=!1,_=!1,C=0,B.length=0,K();const c=P,a=performance.now();await bt("/GANSU-Lite/"),nt=Mt()!=="js"?"wasm":"auto";const o=await Ht(c.basisName),f=tt(c.xyz),x=f.length,l=3*x;let e=new Float64Array(l);for(let z=0;z<x;z++)e[3*z]=f[z].coordinate.x,e[3*z+1]=f[z].coordinate.y,e[3*z+2]=f[z].coordinate.z;const g=f.map(z=>z.atomicNumber);let t=.3,s=null,n=null,d=null,u=!1;const w=q.startsWith("cg-"),j=q==="bfgs"||q==="dfp"||q==="sr1",U=q==="gdiis";for(let z=0;z<X&&!J;z++){const N=g.map((h,p)=>({atomicNumber:h,coordinate:{x:e[3*p],y:e[3*p+1],z:e[3*p+2]},atomIndex:p})),v=new vt(N,o,c.charge),F=$t(v,"RHF"),G=await F.solve({eriBackend:nt}),R=v.numAlphaSpins,M=zt(v.primitiveShells,v.atoms,v.cgtoNormalizationFactors,v.numBasis,R,F.density,F.coefficients,F.orbitalEnergies).total;st(M,e,x);let O=0,rt=0;for(let h=0;h<l;h++)O=Math.max(O,Math.abs(M[h])),rt+=M[h]*M[h];if(b.push({energy:G,maxForce:O,rmsForce:Math.sqrt(rt/l),coords:Array.from(e),gradient:Array.from(M)}),C=b.length-1,Gt(z+1,O),et(),gt(),ot(),O<V){_=!0;break}if(b.length>=2&&G>b[b.length-2].energy&&(!j||!u)){const h=b[b.length-2];for(let p=0;p<l;p++)e[p]=h.coords[p];t*=.5,n=null,await new Promise(p=>setTimeout(p,0));continue}let H=new Float64Array(l);if(q==="sd")for(let h=0;h<l;h++)H[h]=-M[h];else if(w){let h=0;s&&n&&(h=kt(q,M,s,n));for(let m=0;m<l;m++)H[m]=-M[m]+h*(n?n[m]:0);let p=0;for(let m=0;m<l;m++)p+=H[m]*M[m];if(p>0)for(let m=0;m<l;m++)H[m]=-M[m]}else if(j){if(!d){d=new Float64Array(l*l);for(let p=0;p<l;p++)d[p*l+p]=1}if(s&&b.length>=2){const p=b[b.length-2].coords,m=new Float64Array(l),y=new Float64Array(l);for(let S=0;S<l;S++)m[S]=e[S]-p[S],y[S]=M[S]-s[S];st(m,e,x);let E=0,Z=0,at=0;for(let S=0;S<l;S++)E+=m[S]*y[S],Z+=y[S]*y[S],at+=m[S]*m[S];if(E>1e-14&&at>1e-20){if(!u&&Z>1e-14){const S=E/Z;for(let W=0;W<l;W++)d[W*l+W]=S;u=!0}d=At(q,d,m,y,l)}}for(let p=0;p<l;p++){let m=0;for(let y=0;y<l;y++)m+=d[p*l+y]*-M[y];H[p]=m}let h=0;for(let p=0;p<l;p++)h+=H[p]*M[p];if(h>0){d=new Float64Array(l*l);for(let p=0;p<l;p++)d[p*l+p]=1;u=!1;for(let p=0;p<l;p++)H[p]=-M[p]}}else if(U){B.push({coords:Float64Array.from(e),grad:Float64Array.from(M)}),B.length>Nt&&B.shift();const h=Et(l);if(h){for(let y=0;y<l;y++)H[y]=h[y]-e[y];const p=Math.sqrt(H.reduce((y,E)=>y+E*E,0));let m=0;for(let y=0;y<l;y++)m+=H[y]*M[y];if(p<1e-12||m>0)for(let y=0;y<l;y++)H[y]=-M[y];else{const y=new Float64Array(l);for(let E=0;E<l;E++)y[E]=e[E]+H[E];e=y,s=Float64Array.from(M),n=H,await new Promise(E=>setTimeout(E,0));continue}}else for(let p=0;p<l;p++)H[p]=-M[p]}const Y=.5,Q=Math.sqrt(H.reduce((h,p)=>h+p*p,0));if(j&&u){const h=Q>Y?Y/Q:1,p=new Float64Array(l);for(let m=0;m<l;m++)p[m]=e[m]+h*H[m];e=p}else{const h=Math.min(t,Y)/(Q||1),p=new Float64Array(l);for(let m=0;m<l;m++)p[m]=e[m]+h*H[m];e=p}(!j||!u)&&(t=Math.min(t*1.2,Y)),s=Float64Array.from(M),n=H,await new Promise(h=>setTimeout(h,0))}pt=performance.now()-a,L=!1,T=!0,K()}let lt=!1;function Ct(){if(lt)return;lt=!0;const c=document.createElement("style");c.textContent=`
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
  `,document.head.appendChild(c)}const A=document.getElementById("app");function K(){D=null,A.innerHTML=`
    <div class="opt-page">
      ${xt("geomopt")}
      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${$("gopt.molecule")}</h2>
          <div class="opt-scenario-grid" id="scen-grid"></div>

          <h2>${$("gopt.algorithm")}</h2>
          <div id="algo-area"></div>

          <div style="margin:6px 0 0;">
            <label style="font-size:0.72rem;">${$("gopt.maxIter")}
              <input type="number" id="max-iter" min="5" max="200" value="${X}" style="width:50px;margin-left:4px;">
            </label>
          </div>

          <button id="run-btn" class="opt-run-btn" ${L?"disabled":""}>
            ${L?$("gopt.running"):$("gopt.run")}
          </button>
          ${L?`<button id="stop-btn" class="opt-stop-btn">${$("gopt.stop")}</button>`:""}
          <div id="progress-area"></div>
        </div>
        <div class="opt-panel opt-graph-panel">
          <div id="graph-container">
            ${!T&&!L&&b.length===0?`<p class="opt-hint">${$("gopt.waiting")}</p>`:""}
          </div>
          <div id="slider-area"></div>
          <div id="mol-vis"></div>
          <div id="summary-area"></div>
        </div>
      </div>
    </div>`,Ct();const c=A.querySelector("#scen-grid");for(const i of ct){const r=document.createElement("div");r.className="opt-scenario-card"+(i.id===P.id?" selected":""),r.innerHTML=`<strong>${$(i.labelKey)}</strong><span class="conv-desc">${$(i.descKey)}</span>`,r.addEventListener("click",()=>{L||(P=i,b=[],T=!1,_=!1,C=0,D=null,K())}),c.appendChild(r)}const a=A.querySelector("#algo-area");for(const i of dt){const r=document.createElement("div");r.className="gopt-algo-group",r.innerHTML=`<div class="gopt-algo-label">${i.label}</div><div class="gopt-algo-row">${i.items.map(o=>`<button class="gopt-algo-btn${o.id===q?" active":""}" data-algo="${o.id}">${o.name}</button>`).join("")}</div>`,a.appendChild(r)}for(const i of A.querySelectorAll(".gopt-algo-btn"))i.addEventListener("click",()=>{if(!L){q=i.dataset.algo;for(const r of A.querySelectorAll(".gopt-algo-btn"))r.classList.remove("active");i.classList.add("active")}});A.querySelector("#nav-theme").addEventListener("click",()=>{ut(),K()}),A.querySelector("#nav-lang").addEventListener("click",()=>{ht(),K()}),A.querySelector("#run-btn").addEventListener("click",()=>{L||(X=+(A.querySelector("#max-iter")?.value??50),qt())}),A.querySelector("#stop-btn")?.addEventListener("click",()=>{J=!0}),b.length>0&&(et(),It(),gt(),ot())}function Gt(c,a,i){const r=A.querySelector("#progress-area");if(!r)return;const o=`Step ${c} | Max|F| = ${a.toExponential(2)}`;r.innerHTML=`
    <div class="opt-progress"><div class="opt-progress-bar" style="width:${Math.min(100,c/X*100).toFixed(0)}%"></div></div>
    <div class="opt-progress-text">${o}</div>`}function It(){const c=A.querySelector("#summary-area");if(!c||b.length===0)return;const a=b[0],i=b[b.length-1],r=_?"var(--color-accent)":"var(--color-error)",o=_?`✅ ${$("gopt.converged")}`:T?`⚠️ ${$("gopt.notConverged")}`:$("gopt.running"),f=dt.flatMap(t=>t.items).find(t=>t.id===q),x=(t,s)=>`<tr><td style="color:var(--color-text-dim)">${t}</td><td>${s}</td></tr>`,l=.529177210903,e=tt(P.xyz).map(t=>t.atomicNumber),g=t=>{const s=[];for(let n=0;n<e.length;n++){const d=Ft(e[n]),u=(t[3*n]*l).toFixed(6).padStart(11),w=(t[3*n+1]*l).toFixed(6).padStart(11),j=(t[3*n+2]*l).toFixed(6).padStart(11);s.push(`${d.padEnd(2)}${u}${w}${j}`)}return s.join(`
`)};c.innerHTML=`<div class="opt-summary"><h3 style="color:${r}">${o}</h3><table>
    ${x($("gopt.algorithm"),f?.name??q)}
    ${x($("gopt.steps"),String(b.length))}
    ${x($("gopt.initEnergy"),a.energy.toFixed(8)+" Eh")}
    ${x($("gopt.finalEnergy"),i.energy.toFixed(8)+" Eh")}
    ${x($("gopt.energyChange"),(i.energy-a.energy).toFixed(8)+" Eh")}
    ${x($("gopt.maxForce"),i.maxForce.toExponential(4)+" Eh/bohr")}
    ${x($("gopt.threshold"),V.toExponential(1)+" Eh/bohr")}
    ${T?x($("gopt.time"),(pt/1e3).toFixed(1)+" s"):""}
  </table>
  <div class="gopt-xyz-grid">
    <div class="gopt-xyz-block">
      <div class="gopt-xyz-label">Initial (Å)</div>
      <pre class="gopt-xyz-pre">${g(Array.from(a.coords))}</pre>
    </div>
    <div class="gopt-xyz-block">
      <div class="gopt-xyz-label">Final (Å)</div>
      <pre class="gopt-xyz-pre">${g(Array.from(i.coords))}</pre>
    </div>
  </div>
  </div>`}function gt(){const c=A.querySelector("#slider-area");if(!c||b.length<2)return;const a=b[C];c.innerHTML=`
    <div class="gopt-slider">
      <span class="gopt-slider-label">Step ${C+1}/${b.length}</span>
      <input type="range" id="step-slider" min="0" max="${b.length-1}" value="${C}" step="1">
    </div>
    <div class="gopt-step-info">E = ${a.energy.toFixed(8)} Eh &nbsp; Max|F| = ${a.maxForce.toExponential(3)}</div>`,A.querySelector("#step-slider").addEventListener("input",i=>{C=+i.target.value,et(),ot();const r=b[C];c.querySelector(".gopt-slider-label").textContent=`Step ${C+1}/${b.length}`,c.querySelector(".gopt-step-info").innerHTML=`E = ${r.energy.toFixed(8)} Eh &nbsp; Max|F| = ${r.maxForce.toExponential(3)}`})}function et(){const c=A.querySelector("#graph-container");if(!c||b.length<1)return;const a=500,i=155,r={l:70,r:20,t:24,b:28},o=yt();function f(x,l,e,g){const t=e.length,s=1,n=Math.max(t,2);let d,u;if(g){const v=e.map(F=>Math.log10(Math.max(F,1e-15)));d=Math.min(...v),u=Math.max(...v),u-d<.5&&(d-=.25,u+=.25)}else{d=Math.min(...e),u=Math.max(...e);const v=(u-d)*.1||.001;d-=v,u+=v}const w=a-r.l-r.r,j=i-r.t-r.b,U=v=>r.l+(v-s)/(n-s)*w,z=v=>r.t+(1-(v-d)/(u-d))*j;let N=`<svg viewBox="0 0 ${a} ${i}" style="width:100%;max-width:${a}px;background:${o.surface};border-radius:6px;margin-bottom:4px;">`;N+=`<text x="${a/2}" y="14" text-anchor="middle" font-size="10" fill="${o.titleSvg}" font-weight="600">${x}</text>`,N+=`<line x1="${r.l}" y1="${r.t}" x2="${r.l}" y2="${i-r.b}" stroke="${o.grid}"/>`,N+=`<line x1="${r.l}" y1="${i-r.b}" x2="${a-r.r}" y2="${i-r.b}" stroke="${o.grid}"/>`,N+=`<text x="12" y="${i/2}" text-anchor="middle" font-size="7" fill="${o.dim}" transform="rotate(-90,12,${i/2})">${l}</text>`,N+=`<text x="${r.l+w/2}" y="${i-3}" text-anchor="middle" font-size="7" fill="${o.dim}">${$("gopt.xStep")}</text>`;for(let v=0;v<=4;v++){const F=d+(u-d)*v/4,G=z(F);N+=`<line x1="${r.l-3}" y1="${G}" x2="${a-r.r}" y2="${G}" stroke="${o.grid}" stroke-dasharray="2,3"/>`,N+=`<text x="${r.l-5}" y="${G+3}" text-anchor="end" font-size="6.5" fill="${o.dim}">${g?`1e${F.toFixed(1)}`:F.toPrecision(6)}</text>`}if(g){const v=z(Math.log10(V));v>=r.t&&v<=i-r.b&&(N+=`<line x1="${r.l}" y1="${v}" x2="${a-r.r}" y2="${v}" stroke="#ef4444" stroke-dasharray="4,3"/><text x="${a-r.r-2}" y="${v-3}" text-anchor="end" font-size="7" fill="#ef4444">threshold</text>`)}if(t>=1){const v=e.map((F,G)=>`${U(G+1).toFixed(1)},${z(g?Math.log10(Math.max(F,1e-15)):F).toFixed(1)}`);N+=`<polyline points="${v.join(" ")}" fill="none" stroke="${o.occupied}" stroke-width="1.5"/>`;for(let F=0;F<t;F++){const G=g?Math.log10(Math.max(e[F],1e-15)):e[F],R=F===C;N+=`<circle cx="${U(F+1).toFixed(1)}" cy="${z(G).toFixed(1)}" r="${R?5:2.5}" fill="${R?o.accent:o.occupied}" ${R?`stroke="${o.titleSvg}" stroke-width="1.5"`:""} />`}}return N+="</svg>",N}c.innerHTML=f($("gopt.graphEnergy"),$("gopt.yEnergy"),b.map(x=>x.energy),!1)+f($("gopt.graphGrad"),$("gopt.yForce"),b.map(x=>x.maxForce),!0)}function ot(){const c=A.querySelector("#mol-vis");if(!c||b.length===0)return;const a=b[C],i=tt(P.xyz).map(r=>r.atomicNumber);D?D.update(a.coords,a.gradient,V):D=St(c,i,a.coords,a.gradient,V)}mt();ft();K();
