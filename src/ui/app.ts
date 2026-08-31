/** Main UI controller */

import { parseXYZ } from '../core/parseXYZ';
import { BasisSet } from '../core/basisSet';
import { Molecular } from '../core/molecular';
import { buildHF } from '../core/builder';
import type { HFMethod } from '../core/builder';
import { RHF } from '../core/rhf';
import { UHF } from '../core/uhf';
import { ROHF } from '../core/rohf';
import { computeMullikenCharges, computeDipoleMoment, computeS2, computeLowdinCharges, computeWibergBondOrder, computeEnergyComponents } from '../core/properties';
import { atomicNumberToElementName, ELEMENT_NAME_TO_ATOMIC_NUMBER, ATOMIC_GROUND_MULT } from '../core/constants';
import { renderOrbitalDiagram, renderROHFOrbitalDiagram, renderUHFOrbitalDiagram } from './orbitalDiagram';
import { renderConvergenceGraph } from './convergenceGraph';
import { renderDensityHeatmap } from './densityHeatmap';
import { Matrix } from '../linalg/matrix';
import { computeMP2Energy, computeUMP2Energy, computeROMP2Energy } from '../core/mp2';
import { computeCCSDEnergy, computeUCCSDEnergy, computeROCCSDEnergy } from '../core/ccsd';
import { computeMP3Energy, computeUMP3Energy, computeROMP3Energy } from '../core/mp3';
import { computeCIS, type CISExcitedState } from '../core/cis';
import { renderSpectrumChart } from './spectrumChart';
import { computeSOADDensity } from '../core/initialGuess';
import { writeMolden } from '../core/moldenWriter';
import type { WorkerRequest, WorkerResponse, WorkerResult, EriBackend, SCFAccelMethod, SCFAccelParams } from '../core/scfWorker';
import { getCalibration, getWasmPostHFThreshold } from '../core/calibration';
import { isWasmAvailable, getActiveBackend, computeMP2EnergyWasm, computeMP3EnergyWasm, computeCCSDEnergyWasm, computeUMP2EnergyWasm, computeUMP3EnergyWasm, computeUCCSDEnergyWasm, computeROMP2EnergyWasm, computeROMP3EnergyWasm, computeROCCSDEnergyWasm } from '../core/eriWasm';
import { toggleTheme } from './theme';
import { renderMoleculePreview } from './moleculeViewer3D';
import { ProgressTracker, type StepDef } from './progressTracker';
import { t, toggleLang, initLang } from './i18n';
import { renderHeader } from './nav';
import { RIData, generateAutoAuxBasis } from '../core/ri';
import { computeRIMP2Energy, computeRIUMP2Energy } from '../core/riMP2';
import type { DFTConfig } from '../core/builder';
import type { FunctionalName } from '../core/xcFunctional';
import type { GridLevel } from '../core/grid';

interface SampleMolecule {
  xyz: string;
  charge: number;
  /** Optional spin multiplicity (2S+1). If unset, resolved from ATOMIC_GROUND_MULT for single atoms, else kept at user's current selection. */
  mult?: number;
}

const SAMPLE_MOLECULES: SampleMolecule[] = [
  { charge: 0, xyz: `2
Hydrogen (H2)
H  0.0  0.0  0.0
H  0.0  0.0  0.74` },
  { charge: 1, xyz: `2
Helium Hydride Ion (HeH+)
He  0.0  0.0  0.0
H   0.0  0.0  0.77` },
  { charge: 1, xyz: `5
Ammonium (NH4+)
N   0.000   0.000   0.000
H   0.592   0.592   0.592
H  -0.592  -0.592   0.592
H  -0.592   0.592  -0.592
H   0.592  -0.592  -0.592` },
  { charge: 1, xyz: `4
Hydronium (H3O+)
O   0.000   0.000   0.069
H   0.000   0.957  -0.276
H   0.829  -0.478  -0.276
H  -0.829  -0.478  -0.276` },
  { charge: 1, xyz: `2
Nitrosonium (NO+)
N  0.0  0.0  0.0
O  0.0  0.0  1.063` },
  { charge: -1, xyz: `2
Hydroxide (OH-)
O  0.0  0.0  0.0
H  0.0  0.0  0.964` },
  { charge: -1, xyz: `2
Cyanide (CN-)
C  0.0  0.0  0.0
N  0.0  0.0  1.172` },
  { charge: -1, xyz: `5
Borohydride (BH4-)
B   0.000   0.000   0.000
H   0.629   0.629   0.629
H  -0.629  -0.629   0.629
H  -0.629   0.629  -0.629
H   0.629  -0.629  -0.629` },
  { charge: 0, xyz: `2
Lithium Hydride (LiH)
Li  0.0  0.0  0.0
H   0.0  0.0  1.595` },
  { charge: 0, xyz: `2
Hydrogen Fluoride (HF)
H  0.0  0.0  0.0
F  0.0  0.0  0.917` },
  { charge: 0, xyz: `3
Water (H2O)
O   0.00000   0.00000   0.11779
H   0.00000   0.75545  -0.47116
H   0.00000  -0.75545  -0.47116` },
  { charge: 0, xyz: `4
Ammonia (NH3)
N   0.000   0.000   0.128
H   0.000   0.941  -0.298
H   0.815  -0.470  -0.298
H  -0.815  -0.470  -0.298` },
  { charge: 0, xyz: `5
Methane (CH4)
C   0.000   0.000   0.000
H   0.629   0.629   0.629
H  -0.629  -0.629   0.629
H  -0.629   0.629  -0.629
H   0.629  -0.629  -0.629` },
  { charge: 0, xyz: `2
Nitrogen (N2)
N  0.0  0.0  0.0
N  0.0  0.0  1.098` },
  { charge: 0, xyz: `2
Oxygen (O2)
O  0.0  0.0  0.0
O  0.0  0.0  1.208` },
  { charge: 0, xyz: `3
Carbon Dioxide (CO2)
C  0.000  0.000  0.000
O  0.000  0.000  1.162
O  0.000  0.000 -1.162` },
  { charge: 0, xyz: `2
Bromine (Br2)
Br  0.0  0.0  0.0
Br  0.0  0.0  2.281` },
  { charge: 0, xyz: `3
Ozone (O3)
O  0.000  0.000  0.000
O  0.000  1.090  0.671
O  0.000 -1.090  0.671` },
  { charge: 0, xyz: `6
Ethylene (C2H4)
C  0.000  0.000  0.6695
C  0.000  0.000 -0.6695
H  0.923  0.000  1.2425
H -0.923  0.000  1.2425
H  0.923  0.000 -1.2425
H -0.923  0.000 -1.2425` },
  { charge: 0, xyz: `4
Formaldehyde (H2CO)
O   0.000   0.000   0.683
C   0.000   0.000  -0.534
H   0.000   0.926  -1.129
H   0.000  -0.926  -1.129` },
  { charge: 0, xyz: `9
Ethanol (C2H6O)
C   1.187  -0.429   0.000
C   0.000   0.555   0.000
O  -1.218  -0.204   0.000
H  -1.930   0.485   0.000
H   2.126   0.116   0.000
H   1.151  -1.062   0.881
H   1.151  -1.062  -0.881
H   0.063   1.203   0.883
H   0.063   1.203  -0.883` },
  { charge: 0, xyz: `10
Acetone (C3H6O)
C   0.000   0.000   0.198
O   0.000   0.000   1.417
C   0.000   1.304  -0.627
C   0.000  -1.304  -0.627
H  -0.002   2.166   0.032
H   0.002  -2.166   0.032
H   0.881   1.344  -1.262
H  -0.879   1.342  -1.265
H  -0.881  -1.344  -1.262
H   0.879  -1.342  -1.265` },
  { charge: 0, xyz: `11
Pyridine (C5H5N)
C  -0.1802   0.3609  -1.1203
C  -0.1802   1.5593  -0.4079
C  -0.1802   1.5032   0.9869
N  -0.1802   0.3609   1.2902
C  -0.1802  -0.7813   0.9869
C  -0.1802  -0.8374  -0.4079
H  -0.1802   0.3609  -2.2065
H  -0.1802   2.5180  -0.9171
H  -0.1802   2.4213   1.5721
H  -0.1802  -1.6994   1.5721
H  -0.1802  -1.7961  -0.9171` },
  { charge: 0, xyz: `12
Benzene (C6H6)
C   1.3862   0.0000   0.0000
C   0.6931   1.2004   0.0000
C  -0.6931   1.2004   0.0000
C  -1.3862   0.0000   0.0000
C  -0.6931  -1.2004   0.0000
C   0.6931  -1.2004   0.0000
H   2.4618   0.0000   0.0000
H   1.2309   2.1319   0.0000
H  -1.2309   2.1319   0.0000
H  -2.4618   0.0000   0.0000
H  -1.2309  -2.1319   0.0000
H   1.2309  -2.1319   0.0000` },
  { charge: 0, xyz: `13
Alanine (C3H7NO2)
O   1.4573  -1.0438   0.2682
O   1.2492   1.1165  -0.4047
N  -1.4105   1.1507   0.1821
C  -0.7085  -0.1136   0.3937
C  -1.3345  -1.2000  -0.4702
C   0.7470   0.0903   0.0308
H  -0.7666  -0.3737   1.4558
H  -0.8580  -2.1695  -0.2878
H  -2.4023  -1.3127  -0.2521
H  -1.2248  -0.9797  -1.5384
H  -2.3916   1.0420   0.4376
H  -1.4071   1.3875  -0.8099
H   2.4062  -0.9341   0.0447` },
  { charge: 0, xyz: `18
Naphthalene (C10H8)
C   0.000   1.253   1.395
C   0.000   2.421   0.713
C   0.000   2.421  -0.713
C   0.000   1.253  -1.395
C   0.000  -1.253  -1.395
C   0.000  -2.421  -0.713
C   0.000  -2.421   0.713
C   0.000  -1.253   1.395
C   0.000   0.000   0.702
C   0.000   0.000  -0.702
H   0.000   1.244   2.478
H   0.000   3.367   1.240
H   0.000   3.367  -1.240
H   0.000   1.244  -2.478
H   0.000  -1.244  -2.478
H   0.000  -3.367  -1.240
H   0.000  -3.367   1.240
H   0.000  -1.244   2.478` },
  { charge: 0, xyz: `20
Vitamin C (C6H8O6)
O  -0.016   1.342   0.009
C   0.002  -0.004   0.002
C   1.337   1.809  -0.001
O  -0.985  -0.715   0.007
C   1.379  -0.484  -0.014
C   2.210   0.576  -0.016
C   1.593   2.650  -1.253
O   1.762  -1.792  -0.024
O   3.561   0.543  -0.029
O   1.443   1.832  -2.415
C   0.589   3.803  -1.308
O   0.900   4.654  -2.412
H   1.536   2.398   0.894
H   2.606   3.051  -1.220
H   1.025  -2.418  -0.021
H   3.932  -0.349  -0.037
H   0.564   1.440  -2.510
H  -0.418   3.402  -1.430
H   0.642   4.375  -0.382
H   0.303   5.408  -2.508` },
  // ── Single atoms (H–Ar) — multiplicity auto-resolved via ATOMIC_GROUND_MULT ──
  { charge: 0, xyz: `1\nHydrogen (H)\nH   0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nHelium (He)\nHe  0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nLithium (Li)\nLi  0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nBeryllium (Be)\nBe  0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nBoron (B)\nB   0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nCarbon (C)\nC   0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nNitrogen (N)\nN   0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nOxygen (O)\nO   0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nFluorine (F)\nF   0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nNeon (Ne)\nNe  0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nSodium (Na)\nNa  0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nMagnesium (Mg)\nMg  0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nAluminum (Al)\nAl  0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nSilicon (Si)\nSi  0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nPhosphorus (P)\nP   0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nSulfur (S)\nS   0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nChlorine (Cl)\nCl  0.0  0.0  0.0` },
  { charge: 0, xyz: `1\nArgon (Ar)\nAr  0.0  0.0  0.0` },
];

export function initApp(root: HTMLElement) {
  initLang();
  root.innerHTML = `
    ${renderHeader('calc')}
    <div class="controls">
      <div class="panel">
        <h2>${t('panel.molecule')}</h2>
        <label for="sample-select">${t('mol.sample')}</label>
        <select id="sample-select" title="${t('mol.sampleTip')}">
          <option value="">${t('mol.selectDefault')}</option>
        </select>
        <label>${t('mol.xyzLabel')}</label>
        <textarea id="xyz-input" placeholder="${t('mol.xyzPlaceholder')}" title="${t('mol.xyzTip')}"></textarea>
        <div id="xyz-error" class="xyz-error"></div>
        <div id="mol-preview"></div>
        <div id="mol-info" class="mol-info"></div>
      </div>
      <div class="panel">
        <h2>${t('panel.settings')}</h2>
        <label>Theory</label>
        <div class="backend-toggles" id="theory-toggles">
          <button type="button" class="backend-btn active" data-theory="hf">HF</button>
          <button type="button" class="backend-btn" data-theory="dft">DFT</button>
        </div>
        <label>${t('set.basisSet')}</label>
        <div class="backend-toggles" id="basis-toggles">
          <button type="button" class="backend-btn active" data-basis="sto-3g" title="${t('tip.sto3g')}">STO-3G</button>
          <button type="button" class="backend-btn" data-basis="3-21g" title="${t('tip.321g')}">3-21G</button>
          <button type="button" class="backend-btn" data-basis="6-31g" title="${t('tip.631g')}">6-31G</button>
          <button type="button" class="backend-btn" data-basis="cc-pvdz" title="${t('tip.ccpvdz')}">cc-pVDZ</button>
          <button type="button" class="backend-btn" data-basis="aug-cc-pvdz" title="aug-cc-pVDZ (cc-pVDZ + diffuse, for Rydberg / CT / anions)">aug-cc-pVDZ</button>
          <button type="button" class="backend-btn" data-basis="def2-svp" title="def2-SVP (Karlsruhe split-valence + polarisation, ≈ cc-pVDZ)">def2-SVP</button>
          <button type="button" class="backend-btn" data-basis="def2-tzvp" title="def2-TZVP (Karlsruhe triple-zeta + polarisation, common DFT standard)">def2-TZVP</button>
        </div>
        <label>${t('set.method')}</label>
        <div class="backend-toggles" id="method-toggles">
          <button type="button" class="backend-btn active" data-method="RHF" title="${t('tip.rhf')}">RHF</button>
          <button type="button" class="backend-btn" data-method="UHF" title="${t('tip.uhf')}">UHF</button>
          <button type="button" class="backend-btn" data-method="ROHF" title="${t('tip.rohf')}">ROHF</button>
        </div>
        <div id="functional-row" style="display:none;">
          <label>Functional</label>
          <div class="backend-toggles" id="functional-toggles">
            <button type="button" class="backend-btn active" data-functional="SVWN">SVWN</button>
            <button type="button" class="backend-btn" data-functional="BLYP">BLYP</button>
            <button type="button" class="backend-btn" data-functional="PBE">PBE</button>
            <button type="button" class="backend-btn" data-functional="PBE0">PBE0</button>
            <button type="button" class="backend-btn" data-functional="TPSS" title="TPSS meta-GGA (τ-dependent, Tao-Perdew-Staroverov-Scuseria 2003)">TPSS</button>
            <button type="button" class="backend-btn" data-functional="B3LYP">B3LYP</button>
            <button type="button" class="backend-btn" data-functional="CAM-B3LYP" title="CAM-B3LYP range-separated hybrid (α=0.19, β=0.46, ω=0.33). RSH-lite approximation: HF range-split, B3LYP DFT exchange unchanged.">CAM-B3LYP</button>
            <button type="button" class="backend-btn" data-functional="ωB97X-D" title="ωB97X-D long-range corrected (α+β=1, ω=0.30). Use D2 dispersion together. RSH-lite simplification.">ωB97X-D</button>
          </div>
          <label>Grid</label>
          <div class="backend-toggles" id="grid-toggles">
            <button type="button" class="backend-btn" data-grid="coarse">Coarse</button>
            <button type="button" class="backend-btn active" data-grid="medium">Medium</button>
            <button type="button" class="backend-btn" data-grid="fine">Fine</button>
          </div>
        </div>
        <label>${t('set.charge')}</label>
        <div class="backend-toggles" id="charge-toggles">
          <button type="button" class="backend-btn" data-charge="-2" title="${t('tip.charge.-2')}">-2</button>
          <button type="button" class="backend-btn" data-charge="-1" title="${t('tip.charge.-1')}">-1</button>
          <button type="button" class="backend-btn active" data-charge="0" title="${t('tip.charge.0')}">0</button>
          <button type="button" class="backend-btn" data-charge="1" title="${t('tip.charge.1')}">+1</button>
          <button type="button" class="backend-btn" data-charge="2" title="${t('tip.charge.2')}">+2</button>
        </div>
        <label>${t('set.multiplicity')}</label>
        <div class="backend-toggles" id="mult-toggles">
          <button type="button" class="backend-btn active" data-mult="1" title="${t('tip.singlet')}">Singlet</button>
          <button type="button" class="backend-btn" data-mult="2" title="${t('tip.doublet')}">Doublet</button>
          <button type="button" class="backend-btn" data-mult="3" title="${t('tip.triplet')}">Triplet</button>
          <button type="button" class="backend-btn" data-mult="4" title="Quartet (S=3/2)">Quartet</button>
        </div>
        <div id="mult-preview" class="mult-preview"></div>
        <label>${t('set.eriMethod')}</label>
        <div class="backend-toggles" id="eri-method-toggles">
          <button type="button" class="backend-btn active" data-erimethod="auto" title="${t('tip.eriAuto')}">Auto</button>
          <button type="button" class="backend-btn" data-erimethod="stored" title="${t('tip.eriStored')}">Stored</button>
          <button type="button" class="backend-btn" data-erimethod="ri" title="${t('tip.eriRI')}">RI</button>
        </div>
        <div id="aux-basis-row" style="display:none;">
          <label>${t('set.auxBasis')}</label>
          <div class="backend-toggles" id="aux-basis-toggles">
          </div>
        </div>
        <label>${t('set.eriBackend')}</label>
        <div class="backend-toggles" id="backend-toggles">
          <button type="button" class="backend-btn active" data-backend="auto" title="${t('tip.backendAuto')}">Auto</button>
          <button type="button" class="backend-btn" data-backend="wasm" title="${t('tip.wasm')}">WASM</button>
          <button type="button" class="backend-btn" data-backend="js" title="${t('tip.js')}">JS</button>
        </div>
        <span id="cal-info" style="font-size:0.75rem;color:var(--color-text-dim);"></span>
        <label>${t('set.initialGuess')}</label>
        <div class="backend-toggles" id="guess-toggles">
          <button type="button" class="backend-btn active" data-guess="sad" title="${t('tip.sad')}">SAD</button>
          <button type="button" class="backend-btn" data-guess="gwh" title="${t('tip.gwh')}">GWH</button>
          <button type="button" class="backend-btn" data-guess="core" title="${t('tip.coreH')}">Core H</button>
        </div>
        <label>${t('set.scfConvergence')}</label>
        <div class="backend-toggles" id="scf-accel-toggles">
          <button type="button" class="backend-btn active" data-accel="diis" title="${t('tip.diis')}">DIIS</button>
          <button type="button" class="backend-btn" data-accel="damping" title="${t('tip.damp')}">Damp</button>
          <button type="button" class="backend-btn" data-accel="level-shift" title="${t('tip.lshift')}">L-Shift</button>
          <button type="button" class="backend-btn" data-accel="ediis" title="${t('tip.ediis')}">EDIIS</button>
          <button type="button" class="backend-btn" data-accel="adiis-diis" title="${t('tip.adiis')}">ADIIS</button>
        </div>
        <div id="scf-params" class="scf-params">
          <div id="param-diis" class="param-row">
            <label>${t('param.histSize')} <input type="number" id="p-diis-hist" min="2" max="20" value="8" step="1"></label>
            <span class="param-hint">${t('param.diisHint')}</span>
          </div>
          <div id="param-damping" class="param-row" style="display:none">
            <label>${t('param.newRatio')} <input type="number" id="p-damp-alpha" min="0.05" max="0.95" value="0.50" step="0.05"></label>
            <span class="param-hint">${t('param.dampHint')}</span>
          </div>
          <div id="param-level-shift" class="param-row" style="display:none">
            <label>${t('param.shift')} <input type="number" id="p-ls-shift" min="0.1" max="5.0" value="0.50" step="0.1"></label>
            <span class="param-hint">${t('param.shiftHint')}</span>
          </div>
          <div id="param-ediis" class="param-row" style="display:none">
            <label>${t('param.histSize')} <input type="number" id="p-ediis-hist" min="2" max="20" value="8" step="1"></label>
            <span class="param-hint">${t('param.ediisHint')}</span>
          </div>
          <div id="param-adiis-diis" class="param-row" style="display:none">
            <label>${t('param.histSize')} <input type="number" id="p-adiis-hist" min="2" max="20" value="8" step="1"></label>
            <label>${t('param.diisSwitch')} <input type="number" id="p-adiis-thresh" min="0.01" max="1.0" value="0.10" step="0.01"></label>
            <span class="param-hint">${t('param.adiisHint')}</span>
          </div>
        </div>
        <label>${t('set.postHF')}</label>
        <div class="backend-toggles" id="posthf-toggles">
          <button type="button" class="backend-btn active" data-posthf="none" title="${t('tip.none')}">None</button>
          <span class="posthf-group-label" data-group="corr">Correlation:</span>
          <button type="button" class="backend-btn" data-posthf="mp2" data-group="corr" title="${t('tip.mp2')}">MP2</button>
          <button type="button" class="backend-btn" data-posthf="mp3" data-group="corr" title="${t('tip.mp3')}">MP3</button>
          <button type="button" class="backend-btn" data-posthf="ccsd" data-group="corr" title="${t('tip.ccsd')}">CCSD</button>
          <span class="posthf-group-label" data-group="excited">Excited:</span>
          <button type="button" class="backend-btn" data-posthf="cis" data-group="excited" title="${t('tip.cis')}">CIS</button>
          <button type="button" class="backend-btn" data-posthf="adc2" data-group="excited" title="${t('tip.adc2')}">ADC(2)</button>
        </div>
        <div id="cis-options" style="display:none;margin-top:4px;">
          <label>${t('set.spinType')}
            <select id="cis-spin-type">
              <option value="singlet">${t('set.cisSinglet')}</option>
              <option value="triplet">${t('set.excitedTriplet')}</option>
            </select>
          </label>
          <label>${t('set.nStates')} <input type="number" id="cis-nstates" min="1" max="20" value="5" step="1" style="width:48px"></label>
          <label id="full-casida-row" title="Solve the full Casida (RPA) eigenproblem with B coupling instead of TDA. Slightly lower excitation energies; can show triplet instabilities (negative ω²).">
            <input type="checkbox" id="full-casida-check"> Full Casida
          </label>
        </div>
        <div style="margin-top:6px;">
          <label title="${t('set.gradientTip')}"><input type="checkbox" id="gradient-check"> ${t('set.gradient')}</label>
        </div>
        <div style="margin-top:4px;">
          <label title="${t('set.hessianTip')}"><input type="checkbox" id="hessian-check"> ${t('set.hessian')}</label>
        </div>
        <div style="margin-top:4px;">
          <span class="opt-sublabel">${t('set.dispersion')}:</span>
          <div class="backend-toggles" id="dispersion-toggles" style="display:inline-flex;vertical-align:middle;margin-bottom:0;">
            <button type="button" class="backend-btn active" data-disp="none" title="${t('set.dispersionNone')}" style="flex:0 0 auto;white-space:nowrap;">${t('set.dispersionNone')}</button>
            <button type="button" class="backend-btn" data-disp="d2" title="${t('set.dispersionD2Tip')}" style="flex:0 0 auto;white-space:nowrap;">D2</button>
            <button type="button" class="backend-btn" data-disp="d3bj" title="${t('set.dispersionD3BJTip')}" style="flex:0 0 auto;white-space:nowrap;">D3(BJ)</button>
          </div>
        </div>
      </div>
    </div>
    <button id="run-btn" title="${t('btn.runTip')}">${t('btn.run')}</button>
    <progress id="eri-progress" max="100" value="0" style="display:none;"></progress>
    <div class="output-panel" id="output"></div>
    <div id="results-panel" class="results-panel" style="display:none;">
      <h2>${t('panel.results')}</h2>
      <div class="results-grid">
        <div class="panel" id="energy-summary"></div>
        <div class="panel" id="orbital-diagram"></div>
        <div class="panel" id="convergence-graph"></div>
        <div class="panel" id="density-panel"></div>
        <div class="panel" id="mulliken-panel"></div>
        <div class="panel" id="bondorder-panel"></div>
        <div class="panel" id="dipole-panel"></div>
        <div class="panel" id="gradient-panel"></div>
        <div class="panel" id="frequency-panel"></div>
        <div class="panel" id="spectrum-panel"></div>
        <div class="panel" id="molden-panel"></div>
      </div>
    </div>
  `;

  const sampleSelect = root.querySelector<HTMLSelectElement>('#sample-select')!;
  // Populate sample dropdown — molecules first, single atoms grouped in an <optgroup>
  const molGroup = document.createElement('optgroup');
  molGroup.label = 'Molecules';
  const atomGroup = document.createElement('optgroup');
  atomGroup.label = 'Atoms';
  SAMPLE_MOLECULES.forEach((mol, i) => {
    const rawName = mol.xyz.split(/\r?\n/)[1]?.trim() || `Sample ${i}`;
    const atomCount = parseInt(mol.xyz.split(/\r?\n/)[0].trim(), 10) || 0;
    const opt = document.createElement('option');
    opt.value = String(i);
    opt.textContent = t(`sample.${rawName}`) !== `sample.${rawName}` ? t(`sample.${rawName}`) : rawName;
    (atomCount === 1 ? atomGroup : molGroup).appendChild(opt);
  });
  sampleSelect.appendChild(molGroup);
  sampleSelect.appendChild(atomGroup);
  // Default to Water (index 10)
  sampleSelect.value = '10';
  const xyzInput = root.querySelector<HTMLTextAreaElement>('#xyz-input')!;
  const molInfo = root.querySelector<HTMLDivElement>('#mol-info')!;
  const molPreview = root.querySelector<HTMLDivElement>('#mol-preview')!;
  const xyzError = root.querySelector<HTMLDivElement>('#xyz-error')!;

  function validateXYZ(text: string): string[] {
    const errors: string[] = [];
    const lines = text.split(/\r?\n/);
    if (lines.length < 3) { errors.push('Need at least 3 lines (count, comment, atoms)'); return errors; }
    const count = parseInt(lines[0].trim(), 10);
    if (isNaN(count) || count < 1) { errors.push('Line 1: atom count must be a positive integer'); return errors; }
    if (lines.length - 2 < count) errors.push(`Need ${count} atom lines but only ${lines.length - 2} given`);
    const elements = ELEMENT_NAME_TO_ATOMIC_NUMBER;
    for (let i = 2; i < 2 + count && i < lines.length; i++) {
      const line = lines[i];
      if (!line || !line.trim()) { errors.push(`Line ${i+1}: empty`); continue; }
      const parts = line.trim().split(/\s+/);
      if (parts.length < 4) { errors.push(`Line ${i+1}: expected "Symbol x y z"`); continue; }
      if (!(parts[0] in elements)) errors.push(`Line ${i+1}: unknown element "${parts[0]}"`);
      for (let k = 1; k <= 3; k++) {
        if (isNaN(parseFloat(parts[k]))) errors.push(`Line ${i+1}: invalid coordinate "${parts[k]}"`);
      }
    }
    return errors;
  }

  function updateMolPreview() {
    const text = xyzInput.value.trim();
    if (!text) {
      molPreview.innerHTML = '';
      xyzError.textContent = '';
      xyzError.style.display = 'none';
      return;
    }
    const errors = validateXYZ(text);
    if (errors.length > 0) {
      xyzError.innerHTML = errors.slice(0, 5).map(e => `\u26a0 ${e}`).join('<br>');
      xyzError.style.display = 'block';
    } else {
      xyzError.textContent = '';
      xyzError.style.display = 'none';
    }
    renderMoleculePreview(molPreview, text);
    updateMultPreview();
  }

  /** Update the α/β preview + validate (electrons, mult) consistency. */
  function updateMultPreview() {
    const previewEl = root.querySelector<HTMLDivElement>('#mult-preview');
    if (!previewEl) return;
    const text = xyzInput.value.trim();
    if (!text) { previewEl.textContent = ''; previewEl.className = 'mult-preview'; return; }
    let nElec: number;
    try {
      const atoms = parseXYZ(text);
      let totalZ = 0;
      for (const a of atoms) totalZ += a.atomicNumber;
      const charge = parseInt(chargeGroup.value, 10) || 0;
      nElec = totalZ - charge;
    } catch {
      previewEl.textContent = '';
      previewEl.className = 'mult-preview';
      return;
    }
    if (nElec < 1) {
      previewEl.textContent = `\u26a0 ${nElec} electrons (charge too high)`;
      previewEl.className = 'mult-preview error';
      return;
    }
    const mult = parseInt(multGroup.value, 10) || 1;
    const nUnpaired = mult - 1;
    // Validity: (nElec - nUnpaired) must be even and >= 0
    if (nElec < nUnpaired || (nElec - nUnpaired) % 2 !== 0) {
      previewEl.textContent = `\u26a0 ${nElec} electrons cannot be a ${multName(mult)} (mult ${mult})`;
      previewEl.className = 'mult-preview error';
      return;
    }
    const nAlpha = (nElec + nUnpaired) / 2;
    const nBeta = (nElec - nUnpaired) / 2;
    previewEl.textContent = `\u2192 \u03b1${nAlpha} \u03b2${nBeta} (${nElec} electrons)`;
    previewEl.className = 'mult-preview';
  }

  function multName(m: number): string {
    return ['', 'singlet', 'doublet', 'triplet', 'quartet'][m] ?? `mult-${m}`;
  }

  // ── Generic toggle-group helper ──
  function initToggleGroup(containerId: string, dataAttr: string, initial: string, onChange?: (val: string) => void) {
    const container = root.querySelector<HTMLDivElement>(`#${containerId}`)!;
    const btns = container.querySelectorAll<HTMLButtonElement>('.backend-btn');
    const state = { value: initial };
    for (const btn of btns) {
      btn.addEventListener('click', () => {
        if (btn.disabled) return;
        for (const b of btns) b.classList.remove('active');
        btn.classList.add('active');
        state.value = btn.dataset[dataAttr] ?? initial;
        if (onChange) onChange(state.value);
      });
    }
    return { get value() { return state.value; }, set value(v: string) { state.value = v; },
             container, btns };
  }

  const basisGroup = initToggleGroup('basis-toggles', 'basis', 'sto-3g', () => { updatePostHFButtons(); updateAuxBasisRow(); });
  const theoryGroup = initToggleGroup('theory-toggles', 'theory', 'hf', () => { updateTheoryUI(); updatePostHFButtons(); updateCasidaVisibility(postHFGroup.value); });
  const methodGroup = initToggleGroup('method-toggles', 'method', 'RHF', () => updatePostHFButtons());
  const chargeGroup = initToggleGroup('charge-toggles', 'charge', '0', () => { updatePostHFButtons(); updateMultPreview(); });
  const multGroup = initToggleGroup('mult-toggles', 'mult', '1', () => { updatePostHFButtons(); updateMultPreview(); });
  const dispersionGroup = initToggleGroup('dispersion-toggles', 'disp', 'none');
  const auxBasisRow = root.querySelector<HTMLDivElement>('#aux-basis-row')!;
  const auxBasisContainer = root.querySelector<HTMLDivElement>('#aux-basis-toggles')!;
  let auxBasisValue = 'auto';

  function updateAuxBasisRow() {
    const isRI = eriMethodGroup.value === 'ri';
    auxBasisRow.style.display = isRI ? '' : 'none';
    if (!isRI) return;

    const basis = basisGroup.value;
    const hasOptimized = basis in AUX_BASIS_MAP;

    // Rebuild buttons
    auxBasisContainer.innerHTML = '';
    const makeBtn = (value: string, label: string, title: string) => {
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'backend-btn';
      btn.dataset.auxbasis = value;
      btn.title = title;
      btn.textContent = label;
      btn.addEventListener('click', () => {
        for (const b of auxBasisContainer.querySelectorAll<HTMLButtonElement>('.backend-btn')) b.classList.remove('active');
        btn.classList.add('active');
        auxBasisValue = value;
      });
      return btn;
    };

    const autoBtn = makeBtn('auto', 'Auto', t('tip.auxAuto'));
    auxBasisContainer.appendChild(autoBtn);

    if (hasOptimized) {
      const name = AUX_BASIS_MAP[basis];
      const optBtn = makeBtn('optimized', name, t('tip.auxOptimized'));
      auxBasisContainer.appendChild(optBtn);
      // Default to optimized when available
      optBtn.classList.add('active');
      auxBasisValue = 'optimized';
    } else {
      autoBtn.classList.add('active');
      auxBasisValue = 'auto';
    }
  }

  const eriMethodGroup = initToggleGroup('eri-method-toggles', 'erimethod', 'auto', () => { updateAuxBasisRow(); updatePostHFButtons(); });
  const backendGroup = initToggleGroup('backend-toggles', 'backend', 'auto');
  const guessGroup = initToggleGroup('guess-toggles', 'guess', 'sad');
  const accelGroup = initToggleGroup('scf-accel-toggles', 'accel', 'diis', (val) => {
    const paramIds = ['param-diis', 'param-damping', 'param-level-shift', 'param-ediis', 'param-adiis-diis'];
    for (const id of paramIds) {
      const el = root.querySelector<HTMLDivElement>(`#${id}`);
      if (el) el.style.display = id === `param-${val}` ? '' : 'none';
    }
  });
  const cisOptionsDiv = root.querySelector<HTMLDivElement>('#cis-options')!;
  const fullCasidaRow = root.querySelector<HTMLLabelElement>('#full-casida-row')!;
  const fullCasidaCheck = root.querySelector<HTMLInputElement>('#full-casida-check')!;
  function updateCasidaVisibility(val: string) {
    // Full Casida only meaningful for TDDFT (DFT excited states), not CIS or ADC(2)
    const isDFT = theoryGroup.value === 'dft';
    fullCasidaRow.style.display = (val === 'cis' && isDFT) ? '' : 'none';
    if (!(val === 'cis' && isDFT)) fullCasidaCheck.checked = false;
  }
  const postHFGroup = initToggleGroup('posthf-toggles', 'posthf', 'none', (val) => {
    cisOptionsDiv.style.display = (val === 'cis' || val === 'adc2') ? '' : 'none';
    updateCasidaVisibility(val);
  });
  const functionalGroup = initToggleGroup('functional-toggles', 'functional', 'SVWN');
  const gridGroup = initToggleGroup('grid-toggles', 'grid', 'medium');

  const functionalRow = root.querySelector<HTMLDivElement>('#functional-row')!;
  const postHFLabel = root.querySelector<HTMLDivElement>('#posthf-toggles')?.parentElement;
  // Find Post-HF label
  const postHFToggleContainer = root.querySelector<HTMLDivElement>('#posthf-toggles')!;

  function updateTheoryUI() {
    const isDFT = theoryGroup.value === 'dft';
    functionalRow.style.display = isDFT ? '' : 'none';
    // Dynamic label: Post-HF (HF reference) ↔ Excited states (DFT reference)
    const postHFLabelEl = postHFToggleContainer.previousElementSibling;
    if (postHFLabelEl && postHFLabelEl.tagName === 'LABEL') {
      (postHFLabelEl as HTMLElement).textContent = isDFT ? 'Excited states' : t('set.postHF');
    }
    // For DFT: show only CIS (becomes TDDFT) + None; hide correlation group entirely
    // For HF: show everything
    const postHFBtns = postHFToggleContainer.querySelectorAll<HTMLButtonElement>('.backend-btn');
    for (const b of postHFBtns) {
      const v = b.dataset.posthf;
      const grp = b.dataset.group;
      if (isDFT) {
        // Hide MP2/MP3/CCSD (correlation) and ADC(2) — only TDDFT applies for DFT
        if (grp === 'corr' || v === 'adc2') {
          b.style.display = 'none';
        } else {
          b.style.display = '';
        }
        b.disabled = false;
        b.classList.remove('opt-disabled');
        if (v === 'cis') b.textContent = 'TDDFT';
      } else {
        b.style.display = '';
        b.disabled = false;
        b.classList.remove('opt-disabled');
        if (v === 'cis') b.textContent = 'CIS';
      }
    }
    // Hide both group sub-labels for DFT (section label is already "Excited states")
    const groupLabels = postHFToggleContainer.querySelectorAll<HTMLSpanElement>('.posthf-group-label');
    for (const lbl of groupLabels) {
      lbl.style.display = isDFT ? 'none' : '';
    }
    // Update method labels for DFT
    const methodBtns = root.querySelectorAll<HTMLButtonElement>('#method-toggles .backend-btn');
    for (const btn of methodBtns) {
      const m = btn.dataset.method;
      if (isDFT) {
        if (m === 'RHF') btn.textContent = 'RKS';
        else if (m === 'UHF') btn.textContent = 'UKS';
        else if (m === 'ROHF') { btn.textContent = 'ROKS'; btn.disabled = false; }
      } else {
        if (m === 'RHF') btn.textContent = 'RHF';
        else if (m === 'UHF') btn.textContent = 'UHF';
        else if (m === 'ROHF') { btn.textContent = 'ROHF'; btn.disabled = false; }
      }
    }
  }

  // ── Backend support detection ──
  {
    const wasmSupported = typeof WebAssembly !== 'undefined';
    for (const btn of backendGroup.btns) {
      const backend = btn.dataset.backend!;
      const supported = backend === 'auto' || backend === 'js' ||
        (backend === 'wasm' && wasmSupported);
      if (!supported) { btn.disabled = true; btn.classList.add('unsupported'); }
    }
  }

  sampleSelect.addEventListener('change', () => {
    const idx = parseInt(sampleSelect.value, 10);
    const mol = isNaN(idx) ? undefined : SAMPLE_MOLECULES[idx];
    if (mol) {
      xyzInput.value = mol.xyz;
      for (const b of chargeGroup.btns) b.classList.remove('active');
      const chargeBtn = chargeGroup.container.querySelector<HTMLButtonElement>(`[data-charge="${mol.charge}"]`);
      if (chargeBtn) { chargeBtn.classList.add('active'); chargeGroup.value = String(mol.charge); }

      // Resolve multiplicity: explicit > single-atom ground state > keep current
      let targetMult: number | null = mol.mult ?? null;
      if (targetMult === null) {
        const parsed = parseXYZ(mol.xyz);
        if (parsed.length === 1) {
          targetMult = ATOMIC_GROUND_MULT[parsed[0].atomicNumber] ?? null;
        }
      }
      if (targetMult !== null) {
        for (const b of multGroup.btns) b.classList.remove('active');
        const multBtn = multGroup.container.querySelector<HTMLButtonElement>(`[data-mult="${targetMult}"]`);
        if (multBtn) { multBtn.classList.add('active'); multGroup.value = String(targetMult); }
        // Auto-switch RHF → UHF for open-shell
        if (targetMult > 1 && methodGroup.value === 'RHF') {
          for (const b of methodGroup.btns) b.classList.remove('active');
          const uhfBtn = methodGroup.container.querySelector<HTMLButtonElement>('[data-method="UHF"]');
          if (uhfBtn) { uhfBtn.classList.add('active'); methodGroup.value = 'UHF'; }
        }
      }
    }
    updatePostHFButtons();
    updateMolPreview();
  });

  // GBS file cache for quick numBasis lookups
  const gbsCache: Record<string, string> = {};
  async function loadGBS(basisName: string): Promise<string> {
    if (!gbsCache[basisName]) {
      const r = await fetch(`${import.meta.env.BASE_URL}basis/${basisName}.gbs`);
      if (!r.ok) throw new Error(`Failed to load ${basisName}`);
      gbsCache[basisName] = await r.text();
    }
    return gbsCache[basisName];
  }

  const AUX_BASIS_MAP: Record<string, string> = {
    'cc-pvdz': 'cc-pvdz-rifit',
    'cc-pvtz': 'cc-pvtz-rifit',
  };

  async function loadAuxGBS(basisName: string): Promise<string | null> {
    const auxName = AUX_BASIS_MAP[basisName];
    if (!auxName) return null;
    const key = `aux/${auxName}`;
    if (!gbsCache[key]) {
      const r = await fetch(`${import.meta.env.BASE_URL}basis/auxiliary/${auxName}.gbs`);
      if (!r.ok) return null;
      gbsCache[key] = await r.text();
    }
    return gbsCache[key];
  }

  const SHELL_LABELS = ['s', 'p', 'd', 'f', 'g'];
  const SHELL_NBASIS = [1, 3, 6, 10, 15]; // Cartesian basis functions per shell type

  /** Count contracted shells and basis functions by angular momentum type */
  function shellBreakdownDetail(mol: Molecular): { label: string; nBasis: number }[] {
    const seen = new Set<number>();
    const counts: number[] = [];
    for (const ps of mol.primitiveShells) {
      if (!seen.has(ps.basisIndex)) {
        seen.add(ps.basisIndex);
        while (counts.length <= ps.shellType) counts.push(0);
        counts[ps.shellType]++;
      }
    }
    return counts
      .map((c, i) => c > 0 ? { label: SHELL_LABELS[i] ?? `L${i}`, nBasis: c * (SHELL_NBASIS[i] ?? 1), count: c } : null)
      .filter((x): x is { label: string; nBasis: number; count: number } => x !== null)
      .map(x => ({ label: `${x.count}${x.label}`, nBasis: x.nBasis }));
  }

  /** Recompute numBasis and enable/disable post-HF buttons */
  async function updatePostHFButtons() {
    const xyzText = xyzInput.value.trim();
    if (!xyzText) {
      for (const btn of postHFGroup.btns) {
        const m = btn.dataset.posthf ?? 'none';
        if (m === 'none') continue;
        btn.disabled = true;
        btn.classList.add('unsupported');
      }
      molInfo.innerHTML = '';
      return;
    }
    try {
      const gbsText = await loadGBS(basisGroup.value);
      const atoms = parseXYZ(xyzText);
      const basisSet = BasisSet.fromGBS(gbsText);
      const charge = parseInt(chargeGroup.value, 10);
      const mult = parseInt(multGroup.value, 10);
      const betaToAlpha = Math.floor((mult - 1) / 2);
      const mol = new Molecular(atoms, basisSet, charge, betaToAlpha);
      const n = mol.numBasis;

      const isRI = eriMethodGroup.value === 'ri';
      const method = methodGroup.value;
      for (const btn of postHFGroup.btns) {
        const m = btn.dataset.posthf ?? 'none';
        // RI availability: MP2 supported for RHF/UHF, everything else unsupported
        const riUnsupported = isRI && m !== 'none' && (
          m === 'mp3' || m === 'ccsd' ||
          (m === 'mp2' && method === 'ROHF')
        );
        // CIS only available for RHF
        const cisUnsupported = (m === 'cis' || m === 'adc2') && method !== 'RHF';
        const unsupported = riUnsupported || cisUnsupported;
        btn.disabled = unsupported;
        btn.classList.toggle('unsupported', unsupported);
        // If currently selected button becomes disabled, fall back to 'none'
        if (unsupported && btn.classList.contains('active')) {
          btn.classList.remove('active');
          const noneBtn = Array.from(postHFGroup.btns).find((b: HTMLButtonElement) => b.dataset.posthf === 'none');
          if (noneBtn) noneBtn.classList.add('active');
          postHFGroup.value = 'none';
          cisOptionsDiv.style.display = 'none';
        }
      }

      // Gradient & Hessian: available for RHF (HF) and RKS (closed-shell DFT)
      // For DFT, Hessian is numerical (FD of analytic KS gradient) — slower but correct
      const gradCheck = root.querySelector<HTMLInputElement>('#gradient-check')!;
      const hessCheck = root.querySelector<HTMLInputElement>('#hessian-check')!;
      const gradientSupported = method === 'RHF';
      const hessianSupported = method === 'RHF';
      gradCheck.disabled = !gradientSupported;
      hessCheck.disabled = !hessianSupported;
      if (!gradientSupported) gradCheck.checked = false;
      if (!hessianSupported) hessCheck.checked = false;
      // Visually grey out the labels when unsupported
      const gradLabel = gradCheck.parentElement;
      const hessLabel = hessCheck.parentElement;
      if (gradLabel) gradLabel.classList.toggle('opt-disabled', !gradientSupported);
      if (hessLabel) hessLabel.classList.toggle('opt-disabled', !hessianSupported);

      // Build enriched molecule info table
      const nAtoms = atoms.length;
      const nElec = mol.numElectrons;
      const nAlpha = mol.numAlphaSpins;
      const nBeta = mol.numBetaSpins;
      const shellDetail = shellBreakdownDetail(mol);
      const nPrim = mol.primitiveShells.length;
      const nocc = Math.floor(nElec / 2);
      const nvirt = n - nocc;

      // Shell breakdown: "3s(3) + 2p(6) = 9"
      const shellRows = shellDetail
        .map(s => `<span class="mol-row-shell"><b>${s.label}</b><span class="mol-row-dim">(${s.nBasis})</span></span>`)
        .join('<span class="mol-row-dim"> + </span>')
        + `<span class="mol-row-dim"> = ${n}</span>`;

      const elecStr = mult > 1 ? `${nElec} <span class="mol-row-dim">(\u03B1${nAlpha} \u03B2${nBeta})</span>` : `${nElec}`;

      let html = '<div class="mol-row">';
      html += `<span class="mol-row-label">${t('info.atoms')}</span><span class="mol-row-val">${nAtoms}</span>`;
      html += `<span class="mol-row-label">${t('info.electrons')}</span><span class="mol-row-val">${elecStr}</span>`;
      html += '</div><div class="mol-row">';
      html += `<span class="mol-row-label">${t('info.basis')}</span><span class="mol-row-val">${n} <span class="mol-row-dim">(${nPrim} ${t('info.primitives')})</span></span>`;
      html += '</div><div class="mol-row">';
      html += `<span class="mol-row-label">${t('info.orbitals')}</span><span class="mol-row-val">${nocc} ${t('info.occVirt').split(' / ')[0]} / ${nvirt} ${t('info.occVirt').split(' / ')[1]}</span>`;
      html += '</div><div class="mol-row">';
      html += `<span class="mol-row-label">${t('info.shells')}</span><span class="mol-row-val mol-row-shells">${shellRows}</span>`;
      html += '</div>';

      molInfo.innerHTML = html;
    } catch {
      // Parse error — leave buttons as-is
    }
  }

  const runBtn = root.querySelector<HTMLButtonElement>('#run-btn')!;
  const output = root.querySelector<HTMLDivElement>('#output')!;
  const eriProgress = root.querySelector<HTMLProgressElement>('#eri-progress')!;
  const resultsPanel = root.querySelector<HTMLDivElement>('#results-panel')!;

  // Show calibration status
  {
    const calInfo = root.querySelector<HTMLSpanElement>('#cal-info')!;
    const cal = getCalibration();
    if (cal) {
      calInfo.textContent = `WASM PostHF>=${cal.wasmPostHFThreshold ?? 5} basis (${cal.timestamp})`;
    } else {
      calInfo.textContent = '';
    }
  }

  // Theme toggle
  const themeBtn = root.querySelector<HTMLButtonElement>('#nav-theme')!;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  let lastResultData: any = null;
  themeBtn.addEventListener('click', () => {
    const next = toggleTheme();
    themeBtn.textContent = next === 'dark' ? '\u2600' : '\u263E';
    if (lastResultData && resultsPanel.style.display !== 'none') {
      showResultsFromData(lastResultData.result, lastResultData.atoms, lastResultData.elapsed, lastResultData.iterations);
    }
  });

  // Language toggle — rebuild entire UI
  const langBtn = root.querySelector<HTMLButtonElement>('#nav-lang')!;
  langBtn.addEventListener('click', () => {
    toggleLang();
    initApp(root);
  });

  // Live-update preview as user edits XYZ (debounced)
  let xyzDebounce: number | null = null;
  xyzInput.addEventListener('input', () => {
    if (xyzDebounce !== null) clearTimeout(xyzDebounce);
    xyzDebounce = window.setTimeout(() => { updateMolPreview(); xyzDebounce = null; }, 300);
  });

  // Drag & drop .xyz files onto textarea
  xyzInput.addEventListener('dragover', (e) => {
    e.preventDefault();
    xyzInput.classList.add('dragover');
  });
  xyzInput.addEventListener('dragleave', () => {
    xyzInput.classList.remove('dragover');
  });
  xyzInput.addEventListener('drop', (e) => {
    e.preventDefault();
    xyzInput.classList.remove('dragover');
    const file = e.dataTransfer?.files[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = () => {
        xyzInput.value = reader.result as string;
        sampleSelect.value = '';
        updatePostHFButtons();
        updateMolPreview();
      };
      reader.readAsText(file);
    }
  });

  // Pre-fill Water (index 10)
  xyzInput.value = SAMPLE_MOLECULES[10].xyz;
  updatePostHFButtons();
  updateMolPreview();

  runBtn.addEventListener('click', async () => {
    runBtn.disabled = true;
    output.textContent = '';
    resultsPanel.style.display = 'none';
    eriProgress.style.display = 'none';

    const log = (msg: string, cls = '') => {
      const span = document.createElement('span');
      span.className = cls;
      span.textContent = msg + '\n';
      output.appendChild(span);
      output.scrollTop = output.scrollHeight;
    };

    let wasmBackend: 'wasm-simd' | 'wasm' | 'js' = 'js';

    const handleProgress = (msg: string) => {
      const eriMatch = msg.match(/Computing ERIs\.\.\. (\d+)%/);
      if (eriMatch) {
        eriProgress.style.display = 'block';
        eriProgress.value = parseInt(eriMatch[1], 10);
      } else if (eriProgress.style.display !== 'none') {
        eriProgress.style.display = 'none';
      }
      // Track backend tier from Config message
      if (msg.startsWith('Config: WASM SIMD')) wasmBackend = 'wasm-simd';
      else if (msg.startsWith('Config: WASM')) wasmBackend = 'wasm';
      // Detect backend from progress message markers [WASM] / [JS]
      if (tracker) {
        if (msg.includes('[WASM]')) tracker.setBackend(wasmBackend);
        else if (msg.includes('[JS]')) tracker.setBackend('js');
      }
      log(msg, 'progress');
    };

    const iterHistory: Array<{ iter: number; deltaE: number }> = [];
    let tracker: ProgressTracker | null = null;

    try {
      const xyzText = xyzInput.value.trim();
      if (!xyzText) { log(t('log.noXYZ'), 'error'); return; }

      const basisName = basisGroup.value;
      const method = methodGroup.value as HFMethod;
      const charge = parseInt(chargeGroup.value, 10);
      const multiplicity = parseInt(multGroup.value, 10);
      const postHF = postHFGroup.value;
      const runMP2 = postHF === 'mp2';
      const runMP3 = postHF === 'mp3';
      const runCCSD = postHF === 'ccsd';
      const runCIS = postHF === 'cis';
      const runADC2 = postHF === 'adc2';
      const runGradient = root.querySelector<HTMLInputElement>('#gradient-check')!.checked;
      const runHessian = root.querySelector<HTMLInputElement>('#hessian-check')!.checked;
      const dispersion = (dispersionGroup.value as 'none' | 'd2' | 'd3bj');
      const runD2 = dispersion === 'd2';
      const excitedTriplet = (root.querySelector<HTMLSelectElement>('#cis-spin-type')!.value === 'triplet');
      const cisNStates = +(root.querySelector<HTMLInputElement>('#cis-nstates')!.value) || 5;
      const fullCasida = fullCasidaCheck.checked;
      const initialGuess = guessGroup.value as 'core' | 'sad' | 'gwh';
      const eriBackend = backendGroup.value as EriBackend;
      const eriMethod = eriMethodGroup.value; // 'auto' | 'stored' | 'ri'
      const useRI = eriMethod === 'ri';
      const auxBasisMode = auxBasisValue;
      const scfAccelMethod = accelGroup.value as SCFAccelMethod;
      const scfAccelParams: SCFAccelParams = {
        diisHistory: +(root.querySelector<HTMLInputElement>('#p-diis-hist')!.value),
        dampingAlpha: +(root.querySelector<HTMLInputElement>('#p-damp-alpha')!.value),
        levelShift: +(root.querySelector<HTMLInputElement>('#p-ls-shift')!.value),
        ediisHistory: +(root.querySelector<HTMLInputElement>('#p-ediis-hist')!.value),
        adiisHistory: +(root.querySelector<HTMLInputElement>('#p-adiis-hist')!.value),
        adiisThreshold: +(root.querySelector<HTMLInputElement>('#p-adiis-thresh')!.value),
      };
      const isDFT = theoryGroup.value === 'dft';
      const dftConfig: DFTConfig | undefined = isDFT ? {
        functional: functionalGroup.value as FunctionalName,
        gridLevel: gridGroup.value as GridLevel,
      } : undefined;

      const theoryLabel = isDFT ? `${functionalGroup.value}/${method === 'RHF' ? 'RKS' : 'UKS'}` : method;
      const eriLabel = eriMethod === 'ri' ? 'RI' : eriMethod === 'auto' ? 'Auto' : eriBackend;
      log(`Theory: ${theoryLabel}  Basis: ${basisName}  Charge: ${charge}  Mult: ${multiplicity}  Guess: ${initialGuess.toUpperCase()}  ERI: ${eriLabel}  SCF: ${scfAccelMethod.toUpperCase()}${isDFT ? `  Grid: ${gridGroup.value}` : ''}${runMP2 ? '  +MP2' : ''}${runMP3 ? '  +MP3' : ''}${runCCSD ? '  +CCSD' : ''}${runCIS ? `  +${isDFT ? 'TDDFT' : 'CIS'}(${excitedTriplet ? 'T' : 'S'})` : ''}${runADC2 ? '  +ADC(2)' : ''}${runGradient ? '  +Grad' : ''}${runHessian ? '  +Hessian' : ''}`, 'progress');

      // Load basis set GBS text (needed for both worker and main-thread)
      log(`Loading basis set: ${basisName}...`, 'progress');
      const gbsText = await loadGBS(basisName);

      // Parse atoms for showResultsFromData (element names)
      const atoms = parseXYZ(xyzText);

      // Estimate calculation time and confirm if > 60s
      {
        const basisSet = BasisSet.fromGBS(gbsText);
        const betaToAlpha = Math.floor((multiplicity - 1) / 2);
        const mol = new Molecular(atoms, basisSet, charge, betaToAlpha);
        const N = mol.numBasis;
        const nocc = Math.ceil(mol.numElectrons / 2);
        // Empirical scaling estimates (seconds)
        let est = 2e-7 * N * N * N * 20; // SCF ~20 iterations
        if (runMP2) est += 5e-9 * nocc * N * N * N * N;
        if (runMP3) est += 1e-8 * nocc * nocc * N * N * N * N;
        if (runCCSD) est += 1e-8 * nocc * nocc * N * N * N * N;
        if (runGradient) est += 2e-8 * N * N * N * N;
        if (runHessian) est += 6 * atoms.length * 2e-8 * N * N * N * N; // 6N gradient evals
        if (est > 60) {
          const m = Math.floor(est / 60);
          const s = Math.round(est % 60);
          const timeStr = m > 0 ? `${m}m ${s}s` : `${s}s`;
          if (!confirm(t('confirm.time').replace('{time}', timeStr))) {
            runBtn.disabled = false;
            return;
          }
        }
      }

      // Build step definitions based on config
      const stepDefs: StepDef[] = [
        { id: 'setup', label: t('progress.setup') },
      ];
      if (useRI) {
        stepDefs.push({ id: 'ri-2c', label: t('progress.ri2c') });
        stepDefs.push({ id: 'ri-3c', label: t('progress.ri3c') });
        stepDefs.push({ id: 'ri-b', label: t('progress.riB') });
      } else {
        stepDefs.push({ id: 'integrals', label: t('progress.integrals') });
      }
      stepDefs.push({ id: 'scf', label: t('progress.scf') });
      if (runMP2) stepDefs.push({ id: 'mp2', label: t('progress.mp2') });
      if (runMP3) stepDefs.push({ id: 'mp3', label: t('progress.mp3') });
      if (runCCSD) stepDefs.push({ id: 'ccsd', label: t('progress.ccsd') });
      if (runCIS) stepDefs.push({ id: 'cis', label: isDFT ? 'TDDFT-TDA' : t('progress.cis') });
      if (runADC2) stepDefs.push({ id: 'adc2', label: t('progress.adc2') });
      if (runGradient) stepDefs.push({ id: 'gradient', label: t('progress.gradient') });
      if (runHessian) stepDefs.push({ id: 'hessian', label: t('progress.hessian') });
      stepDefs.push({ id: 'properties', label: t('progress.props') });

      let cancelWorker: (() => void) | null = null;

      tracker = new ProgressTracker(stepDefs, () => {
        cancelWorker?.();
      });

      const handleStep = (id: string, status: 'start' | 'done' | 'update', detail?: string) => {
        if (!tracker) return;
        if (status === 'start') {
          tracker.startStep(id, detail);
          // Set backend badge based on step type
          if (id === 'setup' || id === 'properties') {
            tracker.setBackend('js');
          } else if (id === 'integrals' || id === 'scf' || id === 'ri-2c' || id === 'ri-3c' || id === 'ri-b') {
            tracker.setBackend(wasmBackend);
          } else if (id === 'mp2' || id === 'mp3' || id === 'ccsd') {
            tracker.setBackend(wasmBackend);
          } else if (id === 'cis') {
            tracker.setBackend('js');
          }
        } else if (status === 'update') {
          tracker.updateStep(id, detail ?? '');
        } else {
          tracker.completeStep(id, detail);
        }
      };

      const handleIteration = (iter: number, deltaE: number) => {
        tracker?.updateStep('scf', `Iter ${iter + 1}  dE=${deltaE.toExponential(2)}`);
      };

      const t0 = performance.now();

      // Try Web Worker first, fall back to main thread
      // RI runs on main thread only (needs async auxiliary basis loading)
      const cisOpts = runCIS ? { runCIS: true, excitedTriplet, cisNStates, fullCasida } : undefined;
      let workerResult: WorkerResult;
      if (useRI) {
        workerResult = await runOnMainThread(xyzText, gbsText, method, charge, multiplicity, runMP2, runMP3, runCCSD, initialGuess, eriBackend, useRI, auxBasisMode, scfAccelMethod, scfAccelParams, atoms, log, handleProgress, iterHistory, handleStep, handleIteration, dftConfig, cisOpts);
      } else {
        try {
          const run = runInWorker(xyzText, gbsText, method, charge, multiplicity, runMP2, runMP3, runCCSD, initialGuess, eriBackend, scfAccelMethod, scfAccelParams, log, handleProgress, iterHistory, handleStep, handleIteration, dftConfig, cisOpts, runGradient, runHessian, runADC2, runD2, dispersion);
          cancelWorker = run.cancel;
          workerResult = await run.promise;
        } catch (workerErr) {
          // If cancelled, handle separately
          if (workerErr instanceof Error && workerErr.message === 'Cancelled') throw workerErr;
          log(t('log.workerFallback'), 'progress');
          cancelWorker = null;
          workerResult = await runOnMainThread(xyzText, gbsText, method, charge, multiplicity, runMP2, runMP3, runCCSD, initialGuess, eriBackend, useRI, auxBasisMode, scfAccelMethod, scfAccelParams, atoms, log, handleProgress, iterHistory, handleStep, handleIteration, undefined, cisOpts);
        }
      }

      const elapsed = ((performance.now() - t0) / 1000).toFixed(2);
      log('', '');
      log(`Total Energy: ${workerResult.totalEnergy.toFixed(10)} Hartree`, 'energy');
      if (workerResult.dispersionEnergy != null && workerResult.dispersionMethod && workerResult.dispersionMethod !== 'none') {
        const tag = workerResult.dispersionMethod === 'd2' ? 'D2'
                  : workerResult.dispersionMethod === 'd3bj' ? 'D3(BJ)' : workerResult.dispersionMethod;
        log(`${tag} dispersion: ${workerResult.dispersionEnergy.toFixed(10)} Hartree`, 'energy');
        log(`Total + ${tag}:   ${(workerResult.totalEnergy + workerResult.dispersionEnergy).toFixed(10)} Hartree`, 'energy');
      } else if (workerResult.d2Energy != null) {
        // Backward compat path
        log(`D2 dispersion: ${workerResult.d2Energy.toFixed(10)} Hartree`, 'energy');
        log(`Total + D2:    ${(workerResult.totalEnergy + workerResult.d2Energy).toFixed(10)} Hartree`, 'energy');
      }
      log(`Time: ${elapsed} s`, 'converged');

      lastResultData = { result: workerResult, atoms, elapsed, iterations: iterHistory };
      showResultsFromData(workerResult, atoms, elapsed, iterHistory);

      // Close modal after brief delay so user sees final state
      setTimeout(() => tracker?.close(), 600);

    } catch (e) {
      const isCancelled = e instanceof Error && e.message === 'Cancelled';
      if (isCancelled) {
        tracker?.cancel();
        log(t('log.cancelled'), 'progress');
      } else {
        tracker?.failRemaining();
        log(`Error: ${e instanceof Error ? e.message : String(e)}`, 'error');
        console.error(e);
        setTimeout(() => tracker?.close(), 2000);
      }
    } finally {
      runBtn.disabled = false;
    }
  });

  function runInWorker(
    xyzText: string, basisGBS: string, method: HFMethod, charge: number, multiplicity: number, runMP2: boolean, runMP3: boolean, runCCSD: boolean,
    initialGuess: 'core' | 'sad' | 'gwh', eriBackend: EriBackend, scfAccelMethod: SCFAccelMethod, scfAccelParams: SCFAccelParams,
    log: (msg: string, cls: string) => void,
    handleProgress: (msg: string) => void,
    iterHistory: Array<{ iter: number; deltaE: number }>,
    handleStep: (id: string, status: 'start' | 'done' | 'update', detail?: string) => void,
    handleIteration: (iter: number, deltaE: number) => void,
    dftConfig?: DFTConfig,
    cisOpts?: { runCIS: boolean; excitedTriplet: boolean; cisNStates: number; fullCasida?: boolean },
    runGradient?: boolean,
    runHessian?: boolean,
    runADC2?: boolean,
    runD2?: boolean,
    dispersion?: 'none' | 'd2' | 'd3bj',
  ): { promise: Promise<WorkerResult>; cancel: () => void } {
    let cancelFn: () => void = () => {};
    const promise = new Promise<WorkerResult>((resolve, reject) => {
      const worker = new Worker(
        new URL('../core/scfWorker.ts', import.meta.url),
        { type: 'module' },
      );

      cancelFn = () => {
        worker.terminate();
        reject(new Error('Cancelled'));
      };

      worker.onmessage = (e: MessageEvent<WorkerResponse>) => {
        const msg = e.data;
        switch (msg.type) {
          case 'progress':
            handleProgress(msg.message);
            break;
          case 'iteration': {
            const line = `Iter ${String(msg.iter).padStart(3)}: E = ${msg.energy.toFixed(10)}  dE = ${msg.deltaE.toExponential(4)}`;
            log(line, 'iteration');
            iterHistory.push({ iter: msg.iter, deltaE: msg.deltaE });
            handleIteration(msg.iter, msg.deltaE);
            break;
          }
          case 'step':
            handleStep(msg.id, msg.status, msg.detail);
            break;
          case 'result':
            worker.terminate();
            resolve(msg);
            break;
          case 'error':
            worker.terminate();
            reject(new Error(msg.message));
            break;
        }
      };

      worker.onerror = (e) => {
        worker.terminate();
        reject(new Error(e.message || 'Worker failed'));
      };

      worker.postMessage({
        type: 'run-scf', xyzText, basisGBS, method, charge, multiplicity, runMP2, runMP3, runCCSD, runGradient, runHessian, runD2, dispersion, initialGuess, eriBackend, scfAccelMethod, scfAccelParams,
        calibration: getCalibration(),
        baseUrl: import.meta.env.BASE_URL,
        dftConfig,
        ...(cisOpts?.runCIS ? { runCIS: true, excitedTriplet: cisOpts.excitedTriplet, cisNStates: cisOpts.cisNStates, fullCasida: cisOpts.fullCasida } : {}),
        ...(runADC2 ? { runADC2: true, excitedTriplet: root.querySelector<HTMLSelectElement>('#cis-spin-type')!.value === 'triplet', cisNStates: +(root.querySelector<HTMLInputElement>('#cis-nstates')!.value) || 5 } : {}),
      } satisfies WorkerRequest);
    });
    return { promise, cancel: () => cancelFn() };
  }

  async function runOnMainThread(
    xyzText: string, gbsText: string, method: HFMethod, charge: number, multiplicity: number, runMP2: boolean, runMP3: boolean, runCCSD: boolean,
    initialGuess: 'core' | 'sad' | 'gwh', eriBackend: EriBackend, useRI: boolean, auxBasisMode: string, scfAccelMethod: SCFAccelMethod, scfAccelParams: SCFAccelParams,
    atoms: ReturnType<typeof parseXYZ>,
    log: (msg: string, cls: string) => void,
    handleProgress: (msg: string) => void,
    iterHistory: Array<{ iter: number; deltaE: number }>,
    handleStep: (id: string, status: 'start' | 'done' | 'update', detail?: string) => void,
    handleIteration: (iter: number, deltaE: number) => void,
    dftConfig?: DFTConfig,
    cisOpts?: { runCIS: boolean; excitedTriplet: boolean; cisNStates: number; fullCasida?: boolean },
  ): Promise<WorkerResult> {
    // CCSD progress: extract short detail for the progress tracker popup
    const ccsdDetail = (msg: string): string => {
      const m = msg.match(/iter\s+(\d+).*ΔE\s*=\s*([^\s]+)/);
      if (m) return `Iter ${m[1]}  dE=${m[2]}`;
      if (msg.includes('MP2 initial guess')) return 'MP2 initial guess';
      const cm = msg.match(/converged after (\d+)/);
      if (cm) return `Converged (${cm[1]} iter)`;
      if (msg.includes('MO integral')) return 'MO transform';
      if (msg.includes('semicanonical')) return 'Semicanonical orbitals';
      return msg.replace(/^(U?R?O?CCSD:\s*)/i, '');
    };
    const ccsdProgress = (msg: string) => {
      handleProgress(msg);
      handleStep('ccsd', 'update', ccsdDetail(msg));
    };

    handleStep('setup', 'start');
    const betaToAlpha = Math.floor((multiplicity - 1) / 2);
    const basisSet = BasisSet.fromGBS(gbsText);
    const mol = new Molecular(atoms, basisSet, charge, betaToAlpha);
    log(`Atoms: ${atoms.length}, Basis: ${mol.numBasis}, Electrons: ${mol.numElectrons}`, 'progress');

    const ab = getActiveBackend();
    const mainBackend = ab === 'wasm-simd' ? 'WASM SIMD' : ab === 'wasm' ? 'WASM' : 'JS';
    log(`Config: ${mainBackend}, convergence 1e-8`, 'progress');

    const hf = buildHF(mol, method, dftConfig);
    if (dftConfig) {
      log(`DFT: ${dftConfig.functional}, grid=${dftConfig.gridLevel ?? 'medium'}, ${hf.grid?.length ?? 0} grid points`, 'progress');
    }

    // RI setup: load or generate auxiliary basis and precompute fitting coefficients
    handleStep('setup', 'done', `${atoms.length} atoms, ${mol.numBasis} basis`);

    let riData: RIData | undefined;
    if (useRI && (hf instanceof RHF || hf instanceof UHF || hf instanceof ROHF)) {
      let auxBasis: BasisSet | undefined;
      if (auxBasisMode === 'optimized') {
        const auxGBSText = await loadAuxGBS(basisGroup.value);
        if (auxGBSText) {
          auxBasis = BasisSet.fromGBS(auxGBSText);
          log(`RI: using optimized auxiliary basis (${AUX_BASIS_MAP[basisGroup.value]})`, 'progress');
        } else {
          log('Optimized auxiliary basis not found. Falling back to auto-generated.', 'progress');
        }
      }
      if (!auxBasis) {
        // Auto-generate auxiliary basis from product of primary basis exponents
        const elementNames = mol.atoms.map(a => atomicNumberToElementName(a.atomicNumber));
        auxBasis = generateAutoAuxBasis(basisSet, elementNames);
        log('RI: using auto-generated auxiliary basis (product basis)', 'progress');
      }
      handleStep('ri-2c', 'start');
      let riPhase: 'ri-2c' | 'ri-3c' | 'ri-b' = 'ri-2c';
      riData = await RIData.build(
        mol.primitiveShells, mol.cgtoNormalizationFactors, mol.numBasis,
        mol.atoms, auxBasis,
        (msg) => {
          handleProgress(msg);
          if (riPhase === 'ri-2c' && msg.includes('3-center')) {
            handleStep('ri-2c', 'done');
            handleStep('ri-3c', 'start');
            riPhase = 'ri-3c';
          } else if (riPhase === 'ri-3c' && msg.includes('B matrix')) {
            handleStep('ri-3c', 'done');
            handleStep('ri-b', 'start');
            riPhase = 'ri-b';
          }
          // Update current step detail
          if (riPhase === 'ri-2c') {
            if (msg.includes('2-center')) handleStep('ri-2c', 'update', '(P|Q)');
            else if (msg.includes('Cholesky')) handleStep('ri-2c', 'update', 'Cholesky');
          } else if (riPhase === 'ri-3c') {
            const m = msg.match(/(\d+)%/);
            if (m) handleStep('ri-3c', 'update', `${m[1]}%`);
          }
        },
      );
      handleStep(riPhase, 'done');
      hf.setRIData(riData);
      log(`RI: ${riData.naux} auxiliary basis functions`, 'progress');
    }

    if (initialGuess === 'sad') {
      handleProgress(t('log.sadInit'));
      const sadDensity = computeSOADDensity(mol.atoms, basisSet, mol.numBasis, mol.atomToBasisRange);
      hf.setInitialDensityGuess(sadDensity);
    } else if (initialGuess === 'gwh') {
      hf.setInitialGuessType('gwh');
    }

    const totalEnergy = await hf.solve({
      onIteration: (iter, energy, deltaE) => {
        const line = `Iter ${String(iter).padStart(3)}: E = ${energy.toFixed(10)}  dE = ${deltaE.toExponential(4)}`;
        log(line, 'iteration');
        iterHistory.push({ iter, deltaE });
        handleIteration(iter, deltaE);
      },
      onProgress: (msg) => handleProgress(msg),
      onPhase: (phase, status) => handleStep(phase, status),
      eriBackend,
      scfAccelMethod,
      scfAccelParams,
    });

    let orbitalEnergies: number[] = [];
    let mullikenCharges: number[] = [];
    let dipole = { x: 0, y: 0, z: 0, total: 0, debye: 0 };
    let numOccupied = 0;

    let mp2Energy: number | undefined;
    let mp3Energy: number | undefined;
    let ccsdEnergy: number | undefined;
    let moldenText: string | undefined;
    let cisStates: CISExcitedState[] | undefined;
    let cisIsTriplet = false;
    if (hf instanceof RHF) {
      orbitalEnergies = Array.from(hf.orbitalEnergies);
      numOccupied = mol.numAlphaSpins;
      const charges = computeMullikenCharges(hf.density, hf.overlap, mol.atoms, mol.atomToBasisRange);
      mullikenCharges = Array.from(charges);
      dipole = computeDipoleMoment(hf.density, mol.atoms, mol.primitiveShells, mol.cgtoNormalizationFactors, mol.numBasis);
      const useWasmPH = isWasmAvailable() && mol.numBasis >= getWasmPostHFThreshold();
      if (runMP2) {
        handleStep('mp2', 'start');
        if (riData) {
          handleProgress('Computing RI-MP2 correlation energy...');
          mp2Energy = computeRIMP2Energy(hf.coefficients, hf.orbitalEnergies, riData, numOccupied, mol.numBasis);
        } else if (useWasmPH) {
          handleProgress('Computing MP2 correlation energy (WASM)...');
          mp2Energy = computeMP2EnergyWasm(hf.eriStore.data, hf.coefficients.data, hf.orbitalEnergies, numOccupied, mol.numBasis)
            ?? computeMP2Energy(hf.coefficients, hf.orbitalEnergies, hf.eriStore, numOccupied, mol.numBasis);
        } else {
          handleProgress('Computing MP2 correlation energy...');
          mp2Energy = computeMP2Energy(hf.coefficients, hf.orbitalEnergies, hf.eriStore, numOccupied, mol.numBasis);
        }
        handleStep('mp2', 'done', `${mp2Energy.toFixed(6)} Eh`);
      }
      if (runMP3) {
        if (riData) {
          log('MP3 is not available with RI approximation. Skipped.', 'progress');
        } else {
          handleStep('mp3', 'start');
          handleProgress(`Computing MP3 correlation energy${useWasmPH ? ' (WASM)' : ''}...`);
          if (useWasmPH) {
            const wasmResult = computeMP3EnergyWasm(hf.eriStore.data, hf.coefficients.data, hf.orbitalEnergies, numOccupied, mol.numBasis);
            if (wasmResult) { mp2Energy = wasmResult.mp2; mp3Energy = wasmResult.mp3; }
            else {
              const r = computeMP3Energy(hf.coefficients, hf.orbitalEnergies, hf.eriStore, numOccupied, mol.numBasis, handleProgress);
              mp2Energy = r.mp2; mp3Energy = r.mp3;
            }
          } else {
            const r = computeMP3Energy(hf.coefficients, hf.orbitalEnergies, hf.eriStore, numOccupied, mol.numBasis, handleProgress);
            mp2Energy = r.mp2; mp3Energy = r.mp3;
          }
          handleStep('mp3', 'done', `${mp3Energy?.toFixed(6)} Eh`);
        }
      }
      if (runCCSD) {
        if (riData) {
          log('CCSD is not available with RI approximation. Skipped.', 'progress');
        } else {
          handleStep('ccsd', 'start');
          handleProgress(`Computing CCSD correlation energy${useWasmPH ? ' (WASM)' : ''}...`);
          if (useWasmPH) {
            ccsdEnergy = computeCCSDEnergyWasm(hf.eriStore.data, hf.coefficients.data, hf.orbitalEnergies, numOccupied, mol.numBasis)
              ?? computeCCSDEnergy(hf.coefficients, hf.orbitalEnergies, hf.eriStore, numOccupied, mol.numBasis, ccsdProgress);
          } else {
            ccsdEnergy = computeCCSDEnergy(hf.coefficients, hf.orbitalEnergies, hf.eriStore, numOccupied, mol.numBasis, ccsdProgress);
          }
          handleStep('ccsd', 'done', `${ccsdEnergy.toFixed(6)} Eh`);
        }
      }
      if (cisOpts?.runCIS) {
        handleStep('cis', 'start');
        const nStates = cisOpts.cisNStates;
        const isTriplet = cisOpts.excitedTriplet;
        handleProgress(`Computing CIS ${isTriplet ? 'triplet' : 'singlet'} excited states (${nStates} roots)...`);
        const cisResult = computeCIS(
          hf.coefficients, hf.orbitalEnergies, hf.eriStore,
          numOccupied, mol.numBasis, nStates, isTriplet,
          mol.primitiveShells, mol.cgtoNormalizationFactors,
          (msg) => handleProgress(msg),
        );
        cisStates = cisResult.states;
        cisIsTriplet = isTriplet;
        handleStep('cis', 'done', `${cisStates.length} states`);
      }
      moldenText = writeMolden({
        atoms: mol.atoms, basisSet, coefficients: hf.coefficients,
        orbitalEnergies: hf.orbitalEnergies, numOccupied, numBasis: mol.numBasis,
      });
    } else if (hf instanceof ROHF) {
      orbitalEnergies = Array.from(hf.orbitalEnergies);
      numOccupied = mol.numAlphaSpins;
      const charges = computeMullikenCharges(hf.density, hf.overlap, mol.atoms, mol.atomToBasisRange);
      mullikenCharges = Array.from(charges);
      dipole = computeDipoleMoment(hf.density, mol.atoms, mol.primitiveShells, mol.cgtoNormalizationFactors, mol.numBasis);
      const useWasmPH = isWasmAvailable() && mol.numBasis >= getWasmPostHFThreshold();
      if (runMP2) {
        if (riData) {
          log('ROMP2 is not available with RI approximation. Skipped.', 'progress');
        } else {
          handleStep('mp2', 'start');
          handleProgress(`Computing ROMP2 correlation energy${useWasmPH ? ' (WASM)' : ''}...`);
          if (useWasmPH) {
            mp2Energy = computeROMP2EnergyWasm(hf.eriStore.data, hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
              mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis)
              ?? computeROMP2Energy(hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
                hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis);
          } else {
            mp2Energy = computeROMP2Energy(hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
              hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis);
          }
          handleStep('mp2', 'done', `${mp2Energy.toFixed(6)} Eh`);
        }
      }
      if (runMP3) {
        if (riData) {
          log('MP3 is not available with RI approximation. Skipped.', 'progress');
        } else {
          handleStep('mp3', 'start');
          handleProgress(`Computing ROMP3 correlation energy${useWasmPH ? ' (WASM)' : ''}...`);
          if (useWasmPH) {
            const wasmResult = computeROMP3EnergyWasm(hf.eriStore.data, hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
              mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis);
            if (wasmResult) {
              mp2Energy = wasmResult.mp2;
              mp3Energy = wasmResult.mp3;
            } else {
              const mp3Result = computeROMP3Energy(hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
                hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis, handleProgress);
              mp2Energy = mp3Result.mp2;
              mp3Energy = mp3Result.mp3;
            }
          } else {
            const mp3Result = computeROMP3Energy(hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
              hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis, handleProgress);
            mp2Energy = mp3Result.mp2;
            mp3Energy = mp3Result.mp3;
          }
          handleStep('mp3', 'done', `${mp3Energy?.toFixed(6)} Eh`);
        }
      }
      if (runCCSD) {
        if (riData) {
          log('CCSD is not available with RI approximation. Skipped.', 'progress');
        } else {
          handleStep('ccsd', 'start');
          handleProgress(`Computing ROCCSD correlation energy${useWasmPH ? ' (WASM)' : ''}...`);
          if (useWasmPH) {
            ccsdEnergy = computeROCCSDEnergyWasm(hf.eriStore.data, hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
              mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis)
              ?? computeROCCSDEnergy(hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
                hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis, ccsdProgress);
          } else {
            ccsdEnergy = computeROCCSDEnergy(hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
              hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis, ccsdProgress);
          }
          handleStep('ccsd', 'done', `${ccsdEnergy.toFixed(6)} Eh`);
        }
      }
      moldenText = writeMolden({
        atoms: mol.atoms, basisSet, coefficients: hf.coefficients,
        orbitalEnergies: hf.orbitalEnergies, numOccupied, numBasis: mol.numBasis,
        numOccupiedBeta: mol.numBetaSpins,
      });
    } else if (hf instanceof UHF) {
      orbitalEnergies = Array.from(hf.orbitalEnergiesAlpha);
      numOccupied = mol.numAlphaSpins;
      const charges = computeMullikenCharges(hf.density, hf.overlap, mol.atoms, mol.atomToBasisRange);
      mullikenCharges = Array.from(charges);
      dipole = computeDipoleMoment(hf.density, mol.atoms, mol.primitiveShells, mol.cgtoNormalizationFactors, mol.numBasis);
      const useWasmPH = isWasmAvailable() && mol.numBasis >= getWasmPostHFThreshold();
      if (runMP2) {
        handleStep('mp2', 'start');
        if (riData) {
          handleProgress('Computing RI-UMP2 correlation energy...');
          mp2Energy = computeRIUMP2Energy(
            hf.coefficientsAlpha, hf.coefficientsBeta,
            hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
            riData, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis,
          );
        } else if (useWasmPH) {
          handleProgress('Computing UMP2 correlation energy (WASM)...');
          mp2Energy = computeUMP2EnergyWasm(hf.eriStore.data, hf.coefficientsAlpha.data, hf.coefficientsBeta.data,
            hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
            mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis)
            ?? computeUMP2Energy(hf.coefficientsAlpha, hf.coefficientsBeta,
              hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
              hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis);
        } else {
          handleProgress('Computing UMP2 correlation energy...');
          mp2Energy = computeUMP2Energy(hf.coefficientsAlpha, hf.coefficientsBeta,
            hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
            hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis);
        }
        handleStep('mp2', 'done', `${mp2Energy.toFixed(6)} Eh`);
      }
      if (runMP3) {
        if (riData) {
          log('MP3 is not available with RI approximation. Skipped.', 'progress');
        } else {
          handleStep('mp3', 'start');
          handleProgress(`Computing UMP3 correlation energy${useWasmPH ? ' (WASM)' : ''}...`);
          if (useWasmPH) {
            const wasmResult = computeUMP3EnergyWasm(hf.eriStore.data, hf.coefficientsAlpha.data, hf.coefficientsBeta.data,
              hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
              mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis);
            if (wasmResult) {
              mp2Energy = wasmResult.mp2;
              mp3Energy = wasmResult.mp3;
            } else {
              const mp3Result = computeUMP3Energy(hf.coefficientsAlpha, hf.coefficientsBeta,
                hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
                hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis, handleProgress);
              mp2Energy = mp3Result.mp2;
              mp3Energy = mp3Result.mp3;
            }
          } else {
            const mp3Result = computeUMP3Energy(hf.coefficientsAlpha, hf.coefficientsBeta,
              hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
              hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis, handleProgress);
            mp2Energy = mp3Result.mp2;
            mp3Energy = mp3Result.mp3;
          }
          handleStep('mp3', 'done', `${mp3Energy?.toFixed(6)} Eh`);
        }
      }
      if (runCCSD) {
        if (riData) {
          log('CCSD is not available with RI approximation. Skipped.', 'progress');
        } else {
          handleStep('ccsd', 'start');
          handleProgress(`Computing UCCSD correlation energy${useWasmPH ? ' (WASM)' : ''}...`);
          if (useWasmPH) {
            ccsdEnergy = computeUCCSDEnergyWasm(hf.eriStore.data, hf.coefficientsAlpha.data, hf.coefficientsBeta.data,
              hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
              mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis)
              ?? computeUCCSDEnergy(hf.coefficientsAlpha, hf.coefficientsBeta,
                hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
                hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis, ccsdProgress);
          } else {
            ccsdEnergy = computeUCCSDEnergy(hf.coefficientsAlpha, hf.coefficientsBeta,
              hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
              hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis, ccsdProgress);
          }
          handleStep('ccsd', 'done', `${ccsdEnergy.toFixed(6)} Eh`);
        }
      }
      moldenText = writeMolden({
        atoms: mol.atoms, basisSet, coefficients: hf.coefficientsAlpha,
        orbitalEnergies: hf.orbitalEnergiesAlpha, numOccupied, numBasis: mol.numBasis,
        coefficientsBeta: hf.coefficientsBeta,
        orbitalEnergiesBeta: hf.orbitalEnergiesBeta,
        numOccupiedBeta: mol.numBetaSpins,
      });
    }

    // Density matrix
    let densityMatrix: number[] | undefined;
    let densityMatrixAlpha: number[] | undefined;
    let densityMatrixBeta: number[] | undefined;
    if (hf instanceof RHF) {
      densityMatrix = Array.from(hf.density.data);
    } else if (hf instanceof ROHF) {
      densityMatrix = Array.from(hf.density.data);
      densityMatrixAlpha = Array.from(hf.densityAlphaMatrix.data);
      densityMatrixBeta = Array.from(hf.densityBetaMatrix.data);
    } else if (hf instanceof UHF) {
      densityMatrix = Array.from(hf.density.data);
      densityMatrixAlpha = Array.from(hf.densityAlphaMatrix.data);
      densityMatrixBeta = Array.from(hf.densityBetaMatrix.data);
    }

    // New analysis properties
    handleStep('properties', 'start');
    let s2: { exact: number; computed: number } | undefined;
    let lowdinCharges: number[] = [];
    let bondOrders: { i: number; j: number; order: number }[] = [];
    let energyComponents: { oneElectron: number; twoElectron: number; kinetic: number; nuclearAttraction: number } | undefined;

    if (hf instanceof RHF) {
      lowdinCharges = Array.from(computeLowdinCharges(hf.density, hf.overlap, mol.atoms, mol.atomToBasisRange));
      const n = mol.numBasis;
      const Phalf = new Matrix(n, n);
      for (let k = 0; k < n * n; k++) Phalf.data[k] = hf.density.data[k] * 0.5;
      bondOrders = computeWibergBondOrder(Phalf, Phalf, hf.overlap, mol.atoms, mol.atomToBasisRange);
      energyComponents = computeEnergyComponents(
        hf.density, Phalf, Phalf,
        hf.coreHamiltonian, hf.kinetic, hf.fock, hf.fock,
      );
    } else if (hf instanceof ROHF) {
      lowdinCharges = Array.from(computeLowdinCharges(hf.density, hf.overlap, mol.atoms, mol.atomToBasisRange));
      bondOrders = computeWibergBondOrder(hf.densityAlphaMatrix, hf.densityBetaMatrix, hf.overlap, mol.atoms, mol.atomToBasisRange);
      energyComponents = computeEnergyComponents(
        hf.density, hf.densityAlphaMatrix, hf.densityBetaMatrix,
        hf.coreHamiltonian, hf.kinetic, hf.fockAlphaMatrix, hf.fockBetaMatrix,
      );
    } else if (hf instanceof UHF) {
      lowdinCharges = Array.from(computeLowdinCharges(hf.density, hf.overlap, mol.atoms, mol.atomToBasisRange));
      s2 = computeS2(hf.densityAlphaMatrix, hf.densityBetaMatrix, hf.overlap, mol.numAlphaSpins, mol.numBetaSpins);
      bondOrders = computeWibergBondOrder(hf.densityAlphaMatrix, hf.densityBetaMatrix, hf.overlap, mol.atoms, mol.atomToBasisRange);
      energyComponents = computeEnergyComponents(
        hf.density, hf.densityAlphaMatrix, hf.densityBetaMatrix,
        hf.coreHamiltonian, hf.kinetic, hf.fockAlphaMatrix, hf.fockBetaMatrix,
      );
    }

    // UHF/ROHF beta orbital data
    let orbitalEnergiesBeta: number[] | undefined;
    let numOccupiedBeta: number | undefined;
    if (hf instanceof UHF) {
      orbitalEnergiesBeta = Array.from(hf.orbitalEnergiesBeta);
      numOccupiedBeta = mol.numBetaSpins;
    } else if (hf instanceof ROHF) {
      numOccupiedBeta = mol.numBetaSpins;
    }

    handleStep('properties', 'done');

    return {
      type: 'result',
      method,
      totalEnergy,
      nuclearRepulsion: hf.nuclearRepulsion,
      orbitalEnergies,
      mullikenCharges,
      dipole,
      numOccupied,
      numBasis: mol.numBasis,
      densityMatrix,
      densityMatrixAlpha,
      densityMatrixBeta,
      mp2Energy,
      mp3Energy,
      ccsdEnergy,
      moldenText,
      orbitalEnergiesBeta,
      numOccupiedBeta,
      s2,
      lowdinCharges,
      bondOrders,
      energyComponents,
      cisStates,
      excitedTriplet: cisIsTriplet,
    };
  }

  function showResultsFromData(
    result: WorkerResult,
    atoms: ReturnType<typeof parseXYZ>,
    elapsed: string,
    iterations: Array<{ iter: number; deltaE: number }>,
  ) {
    resultsPanel.style.display = 'block';

    // Energy summary
    const energySummary = root.querySelector<HTMLDivElement>('#energy-summary')!;
    let energyRows = `
        <tr title="${t('res.totalHFTip')}"><td>${t('res.totalHF')}</td><td>${result.totalEnergy.toFixed(10)} Hartree</td></tr>
        <tr title="${t('res.nuclearRepulsionTip')}"><td>${t('res.nuclearRepulsion')}</td><td>${result.nuclearRepulsion.toFixed(10)} Hartree</td></tr>
        <tr title="${t('res.electronicEnergyTip')}"><td>${t('res.electronicEnergy')}</td><td>${(result.totalEnergy - result.nuclearRepulsion).toFixed(10)} Hartree</td></tr>`;
    if (result.energyComponents) {
      const ec = result.energyComponents;
      energyRows += `
        <tr title="${t('res.oneElectronTip')}"><td>${t('res.oneElectron')}</td><td>${ec.oneElectron.toFixed(10)} Hartree</td></tr>
        <tr title="${t('res.kineticTip')}"><td>&nbsp;&nbsp;${t('res.kinetic')}</td><td>${ec.kinetic.toFixed(10)} Hartree</td></tr>
        <tr title="${t('res.nuclearAttractionTip')}"><td>&nbsp;&nbsp;${t('res.nuclearAttraction')}</td><td>${ec.nuclearAttraction.toFixed(10)} Hartree</td></tr>
        <tr title="${t('res.twoElectronTip')}"><td>${t('res.twoElectron')}</td><td>${ec.twoElectron.toFixed(10)} Hartree</td></tr>`;
      // Virial ratio: -V/T (should be ≈2 at equilibrium)
      if (ec.kinetic !== 0) {
        const V = ec.nuclearAttraction + ec.twoElectron + result.nuclearRepulsion;
        const virial = -V / ec.kinetic;
        energyRows += `<tr title="${t('res.virialTip')}"><td>${t('res.virialRatio')}</td><td>${virial.toFixed(6)}</td></tr>`;
      }
    }
    // HOMO-LUMO gap & Koopmans IP/EA
    if (result.orbitalEnergies.length > 0) {
      const nocc = result.numOccupied;
      const isUHF = result.orbitalEnergiesBeta != null && result.orbitalEnergiesBeta.length > 0;
      const HA_TO_EV = 27.2114;
      if (isUHF) {
        const eA = result.orbitalEnergies;
        const eB = result.orbitalEnergiesBeta!;
        const noccB = result.numOccupiedBeta ?? 0;
        const homoA = nocc > 0 ? eA[nocc - 1] : NaN;
        const homoB = noccB > 0 ? eB[noccB - 1] : NaN;
        const lumoA = nocc < eA.length ? eA[nocc] : NaN;
        const lumoB = noccB < eB.length ? eB[noccB] : NaN;
        const homo = Math.max(homoA, homoB);
        const lumo = Math.min(lumoA, lumoB);
        if (!isNaN(homo) && !isNaN(lumo)) {
          const gap = (lumo - homo) * HA_TO_EV;
          energyRows += `<tr title="${t('res.homoLumoTip')}"><td>${t('res.homoLumoGap')}</td><td>${gap.toFixed(4)} eV</td></tr>`;
        }
        if (!isNaN(homo)) energyRows += `<tr title="${t('res.ipTip')}"><td>${t('res.ip')}</td><td>${(-homo * HA_TO_EV).toFixed(4)} eV</td></tr>`;
        if (!isNaN(lumo)) energyRows += `<tr title="${t('res.eaTip')}"><td>${t('res.ea')}</td><td>${(-lumo * HA_TO_EV).toFixed(4)} eV</td></tr>`;
      } else {
        const homo = nocc > 0 ? result.orbitalEnergies[nocc - 1] : NaN;
        const lumo = nocc < result.orbitalEnergies.length ? result.orbitalEnergies[nocc] : NaN;
        if (!isNaN(homo) && !isNaN(lumo)) {
          const gap = (lumo - homo) * HA_TO_EV;
          energyRows += `<tr title="${t('res.homoLumoTip')}"><td>${t('res.homoLumoGap')}</td><td>${gap.toFixed(4)} eV</td></tr>`;
        }
        if (!isNaN(homo)) energyRows += `<tr title="${t('res.ipTip')}"><td>${t('res.ip')}</td><td>${(-homo * HA_TO_EV).toFixed(4)} eV</td></tr>`;
        if (!isNaN(lumo)) energyRows += `<tr title="${t('res.eaTip')}"><td>${t('res.ea')}</td><td>${(-lumo * HA_TO_EV).toFixed(4)} eV</td></tr>`;
      }
    }
    if (result.s2) {
      energyRows += `
        <tr title="${t('res.s2Tip')}"><td>${t('res.s2')}</td><td>${result.s2.computed.toFixed(4)} (${t('res.exact')} ${result.s2.exact.toFixed(4)})</td></tr>`;
    }
    if (result.mp3Energy != null && result.mp2Energy != null) {
      const totalMP3 = result.totalEnergy + result.mp2Energy + result.mp3Energy;
      energyRows += `
        <tr title="${t('res.mp2CorrTip')}"><td>${t('res.mp2Corr')}</td><td>${result.mp2Energy.toFixed(10)} Hartree</td></tr>
        <tr title="${t('res.mp3CorrTip')}"><td>${t('res.mp3Corr')}</td><td>${result.mp3Energy.toFixed(10)} Hartree</td></tr>
        <tr title="${t('res.totalMP3Tip')}"><td>${t('res.totalMP3')}</td><td>${totalMP3.toFixed(10)} Hartree</td></tr>`;
    } else if (result.mp2Energy != null) {
      const totalMP2 = result.totalEnergy + result.mp2Energy;
      energyRows += `
        <tr title="${t('res.mp2CorrTip')}"><td>${t('res.mp2Corr')}</td><td>${result.mp2Energy.toFixed(10)} Hartree</td></tr>
        <tr title="${t('res.totalMP2Tip')}"><td>${t('res.totalMP2')}</td><td>${totalMP2.toFixed(10)} Hartree</td></tr>`;
    }
    if (result.ccsdEnergy != null) {
      const totalCCSD = result.totalEnergy + result.ccsdEnergy;
      energyRows += `
        <tr title="${t('res.ccsdCorrTip')}"><td>${t('res.ccsdCorr')}</td><td>${result.ccsdEnergy.toFixed(10)} Hartree</td></tr>
        <tr title="${t('res.totalCCSDTip')}"><td>${t('res.totalCCSD')}</td><td>${totalCCSD.toFixed(10)} Hartree</td></tr>`;
    }
    if (result.dispersionEnergy != null && result.dispersionMethod && result.dispersionMethod !== 'none') {
      const tag = result.dispersionMethod === 'd2' ? 'D2'
                : result.dispersionMethod === 'd3bj' ? 'D3(BJ)' : result.dispersionMethod;
      const totalWithDisp = result.totalEnergy + result.dispersionEnergy;
      const tip = result.dispersionMethod === 'd3bj'
        ? `D3(BJ)-lite: BJ rational damping + C8 term (no CN-dependent C6)`
        : `Grimme D2 dispersion correction (Eh)`;
      energyRows += `
        <tr title="${tip}"><td>${tag} dispersion</td><td>${result.dispersionEnergy.toFixed(10)} Hartree</td></tr>
        <tr title="Total energy + ${tag} dispersion (Eh)"><td>Total + ${tag}</td><td>${totalWithDisp.toFixed(10)} Hartree</td></tr>`;
    } else if (result.d2Energy != null) {
      const totalWithD2 = result.totalEnergy + result.d2Energy;
      energyRows += `
        <tr title="Grimme D2 dispersion correction (Eh)"><td>D2 dispersion</td><td>${result.d2Energy.toFixed(10)} Hartree</td></tr>
        <tr title="Total energy + D2 dispersion (Eh)"><td>Total + D2</td><td>${totalWithD2.toFixed(10)} Hartree</td></tr>`;
    }
    energyRows += `<tr><td>${t('res.time')}</td><td>${elapsed} s</td></tr>`;
    energySummary.innerHTML = `<h2 title="${t('res.energySummaryTip')}">${t('res.energySummary')}</h2><table class="result-table">${energyRows}</table>`;

    if (result.orbitalEnergies.length > 0) {
      const energies = new Float64Array(result.orbitalEnergies); // Always FP64 for display
      const nocc = result.numOccupied;
      const isUHF = result.orbitalEnergiesBeta != null && result.orbitalEnergiesBeta.length > 0;
      const isROHF = result.method === 'ROHF';
      const noccBeta = result.numOccupiedBeta ?? 0;

      const orbitalDiv = root.querySelector<HTMLDivElement>('#orbital-diagram')!;

      if (isUHF) {
        // UHF: side-by-side alpha/beta tables
        const betaEnergies = new Float64Array(result.orbitalEnergiesBeta!);

        const makeTable = (label: string, e: Float64Array, nOcc: number) => {
          let html = `<table class="result-table"><tr><th>${t('res.orbCol.index')}</th><th>${t('res.orbCol.hartree')}</th><th>${t('res.orbCol.ev')}</th><th></th></tr>`;
          for (let i = 0; i < e.length; i++) {
            const eV = e[i] * 27.2114;
            let lbl = '';
            if (i === nOcc - 1) lbl = 'HOMO';
            else if (i === nOcc) lbl = 'LUMO';
            else if (i < nOcc) lbl = t('res.occ');
            else lbl = t('res.vir');
            const cls = i < nOcc ? 'occ' : 'vir';
            html += `<tr class="${cls}"><td>${i + 1}</td><td>${e[i].toFixed(6)}</td><td>${eV.toFixed(3)}</td><td>${lbl}</td></tr>`;
          }
          html += '</table>';
          return `<div><h3>${label}</h3>${html}</div>`;
        };

        orbitalDiv.innerHTML = `<h2 title="${t('res.orbitalEnergiesTip')}">${t('res.orbitalEnergies')}</h2>
          <div style="display:flex;gap:12px;flex-wrap:wrap;">
            ${makeTable(t('res.alphaSpin'), energies, nocc)}
            ${makeTable(t('res.betaSpin'), betaEnergies, noccBeta)}
          </div>`;

        // SVG orbital diagram (UHF side-by-side)
        const svgDiv = document.createElement('div');
        orbitalDiv.appendChild(svgDiv);
        renderUHFOrbitalDiagram(svgDiv, energies, nocc, betaEnergies, noccBeta);
      } else if (isROHF) {
        // ROHF: single table with docc/SOMO/vir labels
        let orbTable = `<h2 title="${t('res.orbitalEnergiesTip')}">${t('res.orbitalEnergies')}</h2><table class="result-table">`;
        orbTable += `<tr><th>${t('res.orbCol.index')}</th><th>${t('res.orbCol.hartree')}</th><th>${t('res.orbCol.ev')}</th><th></th></tr>`;
        for (let i = 0; i < energies.length; i++) {
          const eV = energies[i] * 27.2114;
          let label = '';
          if (i < noccBeta) label = t('res.docc');
          else if (i < nocc) label = 'SOMO';
          else if (i === nocc) label = 'LUMO';
          else label = t('res.vir');
          if (i === nocc - 1 && label === t('res.docc')) label = 'HOMO';
          const cls = i < nocc ? 'occ' : 'vir';
          orbTable += `<tr class="${cls}"><td>${i + 1}</td><td>${energies[i].toFixed(6)}</td><td>${eV.toFixed(3)}</td><td>${label}</td></tr>`;
        }
        orbTable += '</table>';
        orbitalDiv.innerHTML = orbTable;

        // SVG orbital diagram (ROHF)
        const svgDiv = document.createElement('div');
        orbitalDiv.appendChild(svgDiv);
        renderROHFOrbitalDiagram(svgDiv, energies, nocc, noccBeta);
      } else {
        // RHF: single table
        let orbTable = `<h2 title="${t('res.orbitalEnergiesTip')}">${t('res.orbitalEnergies')}</h2><table class="result-table">`;
        orbTable += `<tr><th>${t('res.orbCol.index')}</th><th>${t('res.orbCol.hartree')}</th><th>${t('res.orbCol.ev')}</th><th></th></tr>`;
        for (let i = 0; i < energies.length; i++) {
          const eV = energies[i] * 27.2114;
          let label = '';
          if (i === nocc - 1) label = 'HOMO';
          else if (i === nocc) label = 'LUMO';
          else if (i < nocc) label = t('res.occ');
          else label = t('res.vir');
          const cls = i < nocc ? 'occ' : 'vir';
          orbTable += `<tr class="${cls}"><td>${i + 1}</td><td>${energies[i].toFixed(6)}</td><td>${eV.toFixed(3)}</td><td>${label}</td></tr>`;
        }
        orbTable += '</table>';
        orbitalDiv.innerHTML = orbTable;

        // SVG orbital diagram (RHF)
        const svgDiv = document.createElement('div');
        orbitalDiv.appendChild(svgDiv);
        renderOrbitalDiagram(svgDiv, energies, nocc);
      }

      // Population analysis (Mulliken + Löwdin)
      const mullikenDiv = root.querySelector<HTMLDivElement>('#mulliken-panel')!;
      const hasLowdin = result.lowdinCharges && result.lowdinCharges.length > 0;
      let mulTable = `<h2 title="${t('res.populationTip')}">${t('res.population')}</h2><table class="result-table">`;
      mulTable += `<tr><th>${t('res.colAtom')}</th><th title="${t('res.mullikenTip')}">${t('res.colMulliken')}</th>${hasLowdin ? `<th title="${t('res.lowdinTip')}">${t('res.colLowdin')}</th>` : ''}</tr>`;
      for (let i = 0; i < atoms.length; i++) {
        const elem = atomicNumberToElementName(atoms[i].atomicNumber);
        mulTable += `<tr><td>${elem}${i + 1}</td><td>${result.mullikenCharges[i].toFixed(4)}</td>`;
        if (hasLowdin) mulTable += `<td>${result.lowdinCharges![i].toFixed(4)}</td>`;
        mulTable += '</tr>';
      }
      mulTable += '</table>';
      mullikenDiv.innerHTML = mulTable;

      // Dipole moment
      const dipoleDiv = root.querySelector<HTMLDivElement>('#dipole-panel')!;
      const d = result.dipole;
      dipoleDiv.innerHTML = `<h2 title="${t('res.dipoleTip')}">${t('res.dipole')}</h2>
        <table class="result-table">
          <tr title="${t('res.dipole.xTip')}"><td>x</td><td>${d.x.toFixed(6)} a.u.</td></tr>
          <tr title="${t('res.dipole.yTip')}"><td>y</td><td>${d.y.toFixed(6)} a.u.</td></tr>
          <tr title="${t('res.dipole.zTip')}"><td>z</td><td>${d.z.toFixed(6)} a.u.</td></tr>
          <tr title="${t('res.dipole.auTip')}"><td>|&mu;|</td><td>${d.total.toFixed(6)} a.u.</td></tr>
          <tr title="${t('res.dipole.dbyTip')}"><td>|&mu;|</td><td>${d.debye.toFixed(4)} Debye</td></tr>
        </table>`;

      // Bond orders
      const bondDiv = root.querySelector<HTMLDivElement>('#bondorder-panel')!;
      if (result.bondOrders && result.bondOrders.length > 0) {
        let bondTable = `<h2 title="${t('res.bondOrdersTip')}">${t('res.bondOrders')}</h2><table class="result-table">`;
        bondTable += `<tr><th>${t('res.colBond')}</th><th title="${t('res.orderTip')}">${t('res.colOrder')}</th></tr>`;
        const sorted = [...result.bondOrders].sort((a, b) => b.order - a.order);
        for (const bo of sorted) {
          const elemI = atomicNumberToElementName(atoms[bo.i].atomicNumber);
          const elemJ = atomicNumberToElementName(atoms[bo.j].atomicNumber);
          bondTable += `<tr><td>${elemI}${bo.i + 1}\u2013${elemJ}${bo.j + 1}</td><td>${bo.order.toFixed(3)}</td></tr>`;
        }
        bondTable += '</table>';
        bondDiv.innerHTML = bondTable;
      } else {
        bondDiv.innerHTML = '';
      }
    }

    // Nuclear gradient panel
    const gradDiv = root.querySelector<HTMLDivElement>('#gradient-panel')!;
    if (result.gradient && result.gradient.length > 0) {
      const grad = result.gradient;
      const nAtom = atoms.length;
      let maxForce = 0, rmsSum = 0;
      for (let i = 0; i < grad.length; i++) {
        maxForce = Math.max(maxForce, Math.abs(grad[i]));
        rmsSum += grad[i] * grad[i];
      }
      const rmsForce = Math.sqrt(rmsSum / grad.length);

      let html = `<h2 title="${t('res.gradientTip')}">${t('res.gradient')}</h2><table class="result-table">`;
      html += `<tr><th>${t('res.colAtom')}</th><th>dE/dx</th><th>dE/dy</th><th>dE/dz</th></tr>`;
      for (let i = 0; i < nAtom; i++) {
        const elem = atomicNumberToElementName(atoms[i].atomicNumber);
        html += `<tr><td>${elem}${i + 1}</td>`;
        for (let d = 0; d < 3; d++) html += `<td>${grad[3 * i + d].toFixed(8)}</td>`;
        html += '</tr>';
      }
      html += '</table>';
      html += `<div style="margin-top:6px;font-size:0.85em;color:#888;">`;
      html += `${t('res.maxForce')}: ${maxForce.toFixed(6)} Eh/bohr &nbsp; ${t('res.rmsForce')}: ${rmsForce.toFixed(6)} Eh/bohr</div>`;
      gradDiv.innerHTML = html;
    } else {
      gradDiv.innerHTML = '';
    }

    // Vibrational frequency panel
    const freqDiv = root.querySelector<HTMLDivElement>('#frequency-panel')!;
    if (result.frequencies && result.frequencies.length > 0) {
      const freqs = result.frequencies;
      const irInt = result.irIntensities ?? [];
      const imagFreqs = freqs.filter(f => f < -50);
      const hasIR = irInt.length > 0;

      // ── Frequency table with IR intensities ──
      let html = `<h2 title="${t('res.frequenciesTip')}">${t('res.frequencies')}</h2><table class="result-table">`;
      html += `<tr><th>#</th><th>${t('res.colFreq')}</th>${hasIR ? '<th>IR (km/mol)</th>' : ''}<th>${t('res.colType')}</th></tr>`;
      let modeIdx = 1;
      for (let fi = 0; fi < freqs.length; fi++) {
        const f = freqs[fi];
        if (Math.abs(f) < 50) continue;
        const type = f < 0 ? 'imaginary' : 'real';
        const ir = hasIR ? irInt[fi] : 0;
        html += `<tr${f < 0 ? ' style="color:var(--color-error)"' : ''}>`;
        html += `<td>${modeIdx++}</td><td>${Math.abs(f).toFixed(1)}</td>`;
        if (hasIR) html += `<td>${ir.toFixed(1)}</td>`;
        html += `<td>${type}</td></tr>`;
      }
      html += '</table>';
      if (result.zpe != null) {
        html += `<div style="margin-top:6px;font-size:0.82em;">ZPE: ${result.zpe.toFixed(6)} Eh (${(result.zpe * 627.509).toFixed(2)} kcal/mol)</div>`;
      }
      if (imagFreqs.length > 0) {
        html += `<div style="margin-top:4px;font-size:0.78em;color:var(--color-error);">${imagFreqs.length} imaginary freq \u2192 not a minimum</div>`;
      }

      // ── IR Spectrum SVG (Lorentzian broadening) ──
      if (hasIR) {
        const vibFreqs: number[] = [], vibIR: number[] = [];
        for (let fi = 0; fi < freqs.length; fi++) {
          if (freqs[fi] > 50 && irInt[fi] > 0.01) { vibFreqs.push(freqs[fi]); vibIR.push(irInt[fi]); }
        }
        if (vibFreqs.length > 0) {
          const W = 460, H = 200, PAD = { l: 50, r: 15, t: 24, b: 36 };
          const pw = W - PAD.l - PAD.r, ph = H - PAD.t - PAD.b;
          const xMin = 400, xMax = 4800; // cm⁻¹ range
          const gamma = 30; // Lorentzian half-width cm⁻¹

          // Build spectrum curve (transmittance-like: peaks go DOWN)
          const nPts = 300;
          const spectrum = new Float64Array(nPts);
          for (let p = 0; p < nPts; p++) {
            const nu = xMin + (xMax - xMin) * p / (nPts - 1);
            let absorbance = 0;
            for (let k = 0; k < vibFreqs.length; k++) {
              const d = nu - vibFreqs[k];
              absorbance += vibIR[k] * gamma * gamma / (d * d + gamma * gamma);
            }
            spectrum[p] = absorbance;
          }
          const maxAbs = Math.max(...spectrum) || 1;

          const sx = (nu: number) => PAD.l + (nu - xMin) / (xMax - xMin) * pw;
          const sy = (a: number) => PAD.t + (a / maxAbs) * ph; // 0=top (100%T), maxAbs=bottom (0%T)

          let svg = `<svg viewBox="0 0 ${W} ${H}" style="width:100%;max-width:${W}px;background:var(--color-surface);border:1px solid var(--color-border);border-radius:6px;margin-top:8px;">`;
          svg += `<text x="${W / 2}" y="15" text-anchor="middle" font-size="11" fill="var(--color-text)" font-weight="600">IR Spectrum</text>`;
          // Axes
          svg += `<line x1="${PAD.l}" y1="${PAD.t}" x2="${PAD.l}" y2="${H - PAD.b}" stroke="var(--color-border)"/>`;
          svg += `<line x1="${PAD.l}" y1="${H - PAD.b}" x2="${W - PAD.r}" y2="${H - PAD.b}" stroke="var(--color-border)"/>`;
          // X axis label & ticks
          svg += `<text x="${PAD.l + pw / 2}" y="${H - 5}" text-anchor="middle" font-size="9" fill="var(--color-text-dim)">Wavenumber (cm\u207b\u00b9)</text>`;
          for (let tick = 500; tick <= 4500; tick += 500) {
            if (tick < xMin || tick > xMax) continue;
            const x = sx(tick);
            svg += `<line x1="${x}" y1="${H - PAD.b}" x2="${x}" y2="${H - PAD.b + 3}" stroke="var(--color-border)"/>`;
            svg += `<text x="${x}" y="${H - PAD.b + 13}" text-anchor="middle" font-size="7" fill="var(--color-text-dim)">${tick}</text>`;
          }
          // Y axis label
          svg += `<text x="8" y="${PAD.t + ph / 2}" text-anchor="middle" font-size="8" fill="var(--color-text-dim)" transform="rotate(-90,8,${PAD.t + ph / 2})">Absorbance</text>`;
          // 100% T line
          svg += `<line x1="${PAD.l}" y1="${PAD.t}" x2="${W - PAD.r}" y2="${PAD.t}" stroke="var(--color-border)" stroke-dasharray="2,3"/>`;

          // Spectrum curve
          const pts: string[] = [];
          for (let p = 0; p < nPts; p++) {
            const nu = xMin + (xMax - xMin) * p / (nPts - 1);
            pts.push(`${sx(nu).toFixed(1)},${sy(spectrum[p]).toFixed(1)}`);
          }
          svg += `<polyline points="${pts.join(' ')}" fill="none" stroke="var(--color-accent)" stroke-width="1.5"/>`;

          // Fill under curve
          svg += `<polygon points="${sx(xMin).toFixed(1)},${PAD.t} ${pts.join(' ')} ${sx(xMax).toFixed(1)},${PAD.t}" fill="var(--color-accent)" opacity="0.15"/>`;

          // Peak labels
          for (let k = 0; k < vibFreqs.length; k++) {
            if (vibIR[k] < maxAbs * 0.05) continue; // skip weak
            const x = sx(vibFreqs[k]), y = sy(vibIR[k] * gamma * gamma / (0 + gamma * gamma));
            svg += `<text x="${x}" y="${Math.max(y - 4, PAD.t + 10)}" text-anchor="middle" font-size="7" fill="var(--color-text)">${vibFreqs[k].toFixed(0)}</text>`;
          }

          svg += '</svg>';
          html += svg;
        } else {
          html += `<div style="font-size:0.78rem;color:var(--color-text-dim);margin:6px 0;padding:8px;background:var(--color-surface-alt);border-radius:6px;">No IR-active modes \u2014 homonuclear molecules have no dipole change during vibration.</div>`;
        }
      }

      // Thermodynamic quantities table
      if (result.thermo) {
        const th = result.thermo;
        const f6 = (v: number) => v.toFixed(6);
        const f2 = (v: number) => v.toFixed(2);
        html += `<h3 style="margin-top:12px;font-size:0.82rem;">${t('res.thermo')} (${th.temperature.toFixed(1)} K, ${th.pressure.toFixed(1)} atm)</h3>`;
        html += '<table class="result-table" style="font-size:0.78rem;">';
        html += `<tr><td>${t('res.thermoZPE')}</td><td>${f2(th.ZPE)} kcal/mol</td></tr>`;
        html += `<tr><td>${t('res.thermoThermal')}</td><td>${f6(th.thermalCorr)} Eh</td></tr>`;
        html += `<tr><td>${t('res.thermoEnthalpy')}</td><td>${f6(th.enthalpyCorr)} Eh</td></tr>`;
        html += `<tr><td>${t('res.thermoGibbs')}</td><td>${f6(th.gibbsCorr)} Eh</td></tr>`;
        html += '<tr><td colspan="2" style="border-top:1px solid var(--color-border);padding-top:4px;"></td></tr>';
        html += `<tr><td>E + ${t('res.thermoThermal')}</td><td>${f6(th.E_tot)} Eh</td></tr>`;
        html += `<tr><td>H (${t('res.thermoH')})</td><td>${f6(th.H)} Eh</td></tr>`;
        html += `<tr><td>G (${t('res.thermoG')})</td><td>${f6(th.G)} Eh</td></tr>`;
        html += '<tr><td colspan="2" style="border-top:1px solid var(--color-border);padding-top:4px;"></td></tr>';
        html += `<tr><td>S (${t('res.thermoS')})</td><td>${f2(th.S_tot)} cal/(mol\u00b7K)</td></tr>`;
        html += `<tr><td>&nbsp;&nbsp;${t('res.thermoSTrans')}</td><td>${f2(th.S_trans)}</td></tr>`;
        html += `<tr><td>&nbsp;&nbsp;${t('res.thermoSRot')}</td><td>${f2(th.S_rot)}</td></tr>`;
        html += `<tr><td>&nbsp;&nbsp;${t('res.thermoSVib')}</td><td>${f2(th.S_vib)}</td></tr>`;
        html += `<tr><td>C\u1d65 (${t('res.thermoCv')})</td><td>${f2(th.Cv_tot)} cal/(mol\u00b7K)</td></tr>`;
        html += '</table>';
      }

      freqDiv.innerHTML = html;
    } else {
      freqDiv.innerHTML = '';
    }

    // SCF convergence graph
    const convDiv = root.querySelector<HTMLDivElement>('#convergence-graph')!;
    if (iterations.length >= 2) {
      renderConvergenceGraph(convDiv, iterations, 1e-8);
    } else {
      convDiv.innerHTML = '';
    }

    // Density matrix heatmap (with alpha/beta tabs for UHF)
    const densDiv = root.querySelector<HTMLDivElement>('#density-panel')!;
    if (result.densityMatrix && result.numBasis > 0) {
      const isUHFDensity = result.densityMatrixAlpha != null && result.densityMatrixBeta != null;
      if (isUHFDensity) {
        // Tabbed view: Total / α / β
        densDiv.innerHTML = '';
        const tabBar = document.createElement('div');
        tabBar.className = 'density-tabs';
        const labels = ['Total (P\u03B1+P\u03B2)', '\u03B1 Spin', '\u03B2 Spin'];
        const matrices = [result.densityMatrix, result.densityMatrixAlpha!, result.densityMatrixBeta!];
        const heatmapContainer = document.createElement('div');

        labels.forEach((label, idx) => {
          const btn = document.createElement('button');
          btn.className = 'density-tab' + (idx === 0 ? ' active' : '');
          btn.textContent = label;
          btn.addEventListener('click', () => {
            tabBar.querySelectorAll('.density-tab').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            renderDensityHeatmap(heatmapContainer, matrices[idx], result.numBasis);
          });
          tabBar.appendChild(btn);
        });

        densDiv.appendChild(tabBar);
        densDiv.appendChild(heatmapContainer);
        renderDensityHeatmap(heatmapContainer, result.densityMatrix, result.numBasis);
      } else {
        renderDensityHeatmap(densDiv, result.densityMatrix, result.numBasis);
      }
    } else {
      densDiv.innerHTML = '';
    }

    // Spectrum panel (CIS excited states)
    const spectrumDiv = root.querySelector<HTMLDivElement>('#spectrum-panel')!;
    if (result.cisStates && result.cisStates.length > 0) {
      renderSpectrumChart(spectrumDiv, result.cisStates, result.excitedTriplet ?? false);
    } else {
      spectrumDiv.innerHTML = '';
    }

    // Molden export panel — MOrbVis featured prominently
    const moldenDiv = root.querySelector<HTMLDivElement>('#molden-panel')!;
    if (result.moldenText) {
      moldenDiv.innerHTML =
        `<h2 title="${t('morbvis.headingTip')}">${t('morbvis.heading')}</h2>` +
        '<div class="morbvis-feature">' +
          `<p class="morbvis-desc">${t('morbvis.desc')}</p>` +
          '<div class="morbvis-feature-actions"></div>' +
        '</div>' +
        '<div class="molden-secondary"></div>';

      // MOrbVis primary button
      const morbvisBtn = document.createElement('button');
      morbvisBtn.className = 'morbvis-primary-btn';
      morbvisBtn.textContent = t('morbvis.open');
      morbvisBtn.addEventListener('click', () => {
        try {
          sessionStorage.setItem('gansu-molden', result.moldenText!);
          const morbvisUrl = new URL('/morbvis/', window.location.origin).href;
          window.open(morbvisUrl, '_blank');
        } catch {
          const blob = new Blob([result.moldenText!], { type: 'text/plain' });
          const url = URL.createObjectURL(blob);
          const a = document.createElement('a');
          a.href = url;
          a.download = 'gansu-web.molden';
          a.click();
          URL.revokeObjectURL(url);
        }
      });
      moldenDiv.querySelector('.morbvis-feature-actions')!.appendChild(morbvisBtn);

      // Molden download — secondary
      const moldenBtn = document.createElement('button');
      moldenBtn.className = 'molden-btn';
      moldenBtn.textContent = t('morbvis.download');
      moldenBtn.addEventListener('click', () => {
        const blob = new Blob([result.moldenText!], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'gansu-web.molden';
        a.click();
        URL.revokeObjectURL(url);
      });
      const secondary = moldenDiv.querySelector('.molden-secondary')!;
      secondary.appendChild(moldenBtn);
    } else {
      moldenDiv.innerHTML = '';
    }
  }
}
