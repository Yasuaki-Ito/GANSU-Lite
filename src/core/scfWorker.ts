/** Web Worker for SCF computation — runs HF off the main thread */

import { parseXYZ } from './parseXYZ';
import { BasisSet } from './basisSet';
import { Molecular } from './molecular';
import { buildHF } from './builder';
import type { HFMethod, DFTConfig } from './builder';
import { RHF } from './rhf';
import { UHF } from './uhf';
import { ROHF } from './rohf';
import { computeMullikenCharges, computeDipoleMoment, computeS2, computeLowdinCharges, computeWibergBondOrder, computeEnergyComponents } from './properties';
import { computeMP2Energy, computeUMP2Energy, computeROMP2Energy } from './mp2';
import { computeDispersion, type DispersionMethod } from './dispersion';
import { computeCCSDEnergy, computeUCCSDEnergy, computeROCCSDEnergy } from './ccsd';
import { computeMP3Energy, computeUMP3Energy, computeROMP3Energy } from './mp3';
import { computeCIS, type CISExcitedState } from './cis';
import { computeTDDFT } from './tddft';
import { computeRHFGradient } from './gradient';
import { computeADC2 } from './adc2';
import { computeHessianAuto, type HessianResult } from './hessian';
import { computeSOADDensity } from './initialGuess';
import { writeMolden } from './moldenWriter';
import { Matrix } from '../linalg/matrix';
import { RIData, generateAutoAuxBasis } from './ri';
import { atomicNumberToElementName } from './constants';
import { initWasm, isWasmAvailable, getActiveBackend, computeMP2EnergyWasm, computeMP3EnergyWasm, computeCCSDEnergyWasm, computeUMP2EnergyWasm, computeUMP3EnergyWasm, computeUCCSDEnergyWasm, computeROMP2EnergyWasm, computeROMP3EnergyWasm, computeROCCSDEnergyWasm } from './eriWasm';
import type { EriBackend } from './hf';
import type { SCFAccelMethod, SCFAccelParams } from './scfAccelerator';
import { setCalibration, getWasmPostHFThreshold, type CalibrationData } from './calibration';

export type { EriBackend, SCFAccelMethod, SCFAccelParams };

export interface WorkerRequest {
  type: 'run-scf';
  xyzText: string;
  basisGBS: string;
  method: HFMethod;
  charge: number;
  multiplicity: number;
  runMP2?: boolean;
  runMP3?: boolean;
  runCCSD?: boolean;
  runCIS?: boolean;
  runADC2?: boolean;
  excitedTriplet?: boolean;
  cisNStates?: number;
  fullCasida?: boolean;
  runGradient?: boolean;
  runHessian?: boolean;
  runD2?: boolean;
  dispersion?: DispersionMethod;
  initialGuess?: 'core' | 'sad' | 'gwh';
  eriBackend?: EriBackend;
  calibration?: CalibrationData | null;
  baseUrl?: string;
  scfAccelMethod?: SCFAccelMethod;
  scfAccelParams?: SCFAccelParams;
  dftConfig?: DFTConfig;
}

export interface WorkerProgress {
  type: 'progress';
  message: string;
}

export interface WorkerIteration {
  type: 'iteration';
  iter: number;
  energy: number;
  deltaE: number;
}

export interface WorkerResult {
  type: 'result';
  method?: HFMethod;
  totalEnergy: number;
  nuclearRepulsion: number;
  orbitalEnergies: number[];
  mullikenCharges: number[];
  dipole: { x: number; y: number; z: number; total: number; debye: number };
  numOccupied: number;
  numBasis: number;
  densityMatrix?: number[];
  densityMatrixAlpha?: number[];
  densityMatrixBeta?: number[];
  mp2Energy?: number;
  mp3Energy?: number;
  ccsdEnergy?: number;
  d2Energy?: number;
  dispersionEnergy?: number;
  dispersionMethod?: DispersionMethod;
  moldenText?: string;
  // UHF-specific: beta orbital energies and occupation
  orbitalEnergiesBeta?: number[];
  numOccupiedBeta?: number;
  // New analysis properties
  s2?: { exact: number; computed: number };
  lowdinCharges?: number[];
  bondOrders?: { i: number; j: number; order: number }[];
  energyComponents?: { oneElectron: number; twoElectron: number; kinetic: number; nuclearAttraction: number };
  cisStates?: CISExcitedState[];
  excitedTriplet?: boolean;
  gradient?: number[];
  frequencies?: number[];    // cm⁻¹
  irIntensities?: number[];  // km/mol
  zpe?: number;              // Hartree
  thermo?: import('./hessian').ThermoData;
}

export interface WorkerStep {
  type: 'step';
  id: string;
  status: 'start' | 'done' | 'update';
  detail?: string;
}

export interface WorkerError {
  type: 'error';
  message: string;
}

export type WorkerResponse = WorkerProgress | WorkerIteration | WorkerStep | WorkerResult | WorkerError;

/** Extract short detail text from CCSD onProgress messages for the progress tracker. */
function ccsdDetailFromMsg(msg: string): string {
  // "CCSD iter  5: E = -0.012345678901, ΔE = 1.2345e-06" → "Iter 5  dE=1.23e-06"
  const m = msg.match(/iter\s+(\d+).*ΔE\s*=\s*([^\s]+)/);
  if (m) return `Iter ${m[1]}  dE=${m[2]}`;
  // "CCSD iter  0: E = ... (MP2 initial guess)" → "MP2 initial guess"
  if (msg.includes('MP2 initial guess')) return 'MP2 initial guess';
  // "CCSD converged after N iterations"
  const cm = msg.match(/converged after (\d+)/);
  if (cm) return `Converged (${cm[1]} iter)`;
  // "CCSD: Full MO integral transformation..." → "MO transform"
  if (msg.includes('MO integral')) return 'MO transform';
  // "semicanonical" → pass through
  if (msg.includes('semicanonical')) return 'Semicanonical orbitals';
  return msg.replace(/^(U?R?O?CCSD:\s*)/i, '');
}

self.onmessage = async (e: MessageEvent<WorkerRequest>) => {
  const { xyzText, basisGBS, method, charge, multiplicity, runMP2, runMP3, runCCSD, runCIS, runADC2, excitedTriplet, cisNStates, fullCasida, runGradient, runHessian, runD2, dispersion, initialGuess, eriBackend, calibration, baseUrl, scfAccelMethod, scfAccelParams, dftConfig } = e.data;

  // Inject calibration data (Worker has no localStorage)
  if (calibration) setCalibration(calibration);

  try {
    const post = (msg: WorkerResponse) => self.postMessage(msg);
    const backend = eriBackend || 'auto';

    post({ type: 'step', id: 'setup', status: 'start' });
    post({ type: 'progress', message: 'Parsing molecule...' });
    const atoms = parseXYZ(xyzText);
    const basisSet = BasisSet.fromGBS(basisGBS);
    const betaToAlpha = Math.floor((multiplicity - 1) / 2);
    const mol = new Molecular(atoms, basisSet, charge, betaToAlpha);
    post({ type: 'progress', message: `Atoms: ${atoms.length}, Basis: ${mol.numBasis}, Electrons: ${mol.numElectrons}` });

    // Initialize backends
    const backends: string[] = [];
    if (backend === 'auto' || backend === 'wasm') {
      const wasmOk = await initWasm(baseUrl);
      if (wasmOk) backends.push('WASM');
    }
    backends.push('JS');

    const activeBackend = getActiveBackend();
    const fockBackend = activeBackend === 'wasm-simd' ? 'WASM SIMD' : activeBackend === 'wasm' ? 'WASM' : 'JS';
    post({ type: 'progress', message: `Config: ${fockBackend}, convergence 1e-8` });

    const hf = buildHF(mol, method, dftConfig);
    if (dftConfig) {
      post({ type: 'progress', message: `DFT: ${dftConfig.functional}, grid=${dftConfig.gridLevel ?? 'medium'}, ${hf.grid?.length ?? 0} grid points` });
    }

    // Auto RI-J for pure DFT (no HF exchange): skip expensive full ERI computation
    const isPureDFT = hf.xcFunctional && hf.xcFunctional.exactExchangeFraction === 0;
    if (isPureDFT) {
      post({ type: 'progress', message: 'Building RI-J auxiliary basis (auto)...' });
      const elementNames = [...new Set(mol.atoms.map(a => atomicNumberToElementName(a.atomicNumber)))];
      const auxBasis = generateAutoAuxBasis(basisSet, elementNames);
      const riData = await RIData.build(
        mol.primitiveShells, mol.cgtoNormalizationFactors, mol.numBasis,
        mol.atoms, auxBasis,
        (msg) => post({ type: 'progress', message: `RI-J: ${msg}` }),
      );
      (hf as unknown as { setRIData: (ri: RIData) => void }).setRIData(riData);
      post({ type: 'progress', message: `RI-J ready: ${riData.naux} aux functions` });
    }

    post({ type: 'step', id: 'setup', status: 'done', detail: `${atoms.length} atoms, ${mol.numBasis} basis` });

    if (initialGuess === 'sad') {
      post({ type: 'progress', message: 'Computing SAD initial density...' });
      const sadDensity = computeSOADDensity(mol.atoms, basisSet, mol.numBasis, mol.atomToBasisRange);
      hf.setInitialDensityGuess(sadDensity);
    } else if (initialGuess === 'gwh') {
      hf.setInitialGuessType('gwh');
    }

    const totalEnergy = await hf.solve({
      onIteration: (iter, energy, deltaE) => {
        post({ type: 'iteration', iter, energy, deltaE });
      },
      onProgress: (msg) => {
        post({ type: 'progress', message: msg });
      },
      onPhase: (phase, status) => {
        post({ type: 'step', id: phase, status });
      },
      eriBackend: backend,
      scfAccelMethod,
      scfAccelParams,
    });

    // Post-SCF properties
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
    let gradient: number[] | undefined;
    let frequencies: number[] | undefined;
    let irIntensities: number[] | undefined;
    let zpe: number | undefined;
    let thermo: import('./hessian').ThermoData | undefined;

    if (hf instanceof RHF) {
      orbitalEnergies = Array.from(hf.orbitalEnergies);
      numOccupied = mol.numAlphaSpins;
      const charges = computeMullikenCharges(hf.density, hf.overlap, mol.atoms, mol.atomToBasisRange);
      mullikenCharges = Array.from(charges);
      dipole = computeDipoleMoment(hf.density, mol.atoms, mol.primitiveShells, mol.cgtoNormalizationFactors, mol.numBasis);

      const useWasmPH = isWasmAvailable() && mol.numBasis >= getWasmPostHFThreshold();

      if (runMP2) {
        post({ type: 'step', id: 'mp2', status: 'start' });
        post({ type: 'progress', message: `Computing MP2 correlation energy${useWasmPH ? ' (WASM)' : ''}...` });
        if (useWasmPH) {
          mp2Energy = computeMP2EnergyWasm(hf.eriStore.data, hf.coefficients.data, hf.orbitalEnergies, numOccupied, mol.numBasis)
            ?? computeMP2Energy(hf.coefficients, hf.orbitalEnergies, hf.eriStore, numOccupied, mol.numBasis);
        } else {
          mp2Energy = computeMP2Energy(hf.coefficients, hf.orbitalEnergies, hf.eriStore, numOccupied, mol.numBasis);
        }
        post({ type: 'progress', message: `MP2 correlation energy: ${mp2Energy.toFixed(10)} Hartree` });
        post({ type: 'step', id: 'mp2', status: 'done', detail: `${mp2Energy.toFixed(6)} Eh` });
      }

      if (runMP3) {
        post({ type: 'step', id: 'mp3', status: 'start' });
        post({ type: 'progress', message: `Computing MP3 correlation energy${useWasmPH ? ' (WASM)' : ''}...` });
        if (useWasmPH) {
          const wasmResult = computeMP3EnergyWasm(hf.eriStore.data, hf.coefficients.data, hf.orbitalEnergies, numOccupied, mol.numBasis);
          if (wasmResult) {
            mp2Energy = wasmResult.mp2;
            mp3Energy = wasmResult.mp3;
          } else {
            const mp3Result = computeMP3Energy(hf.coefficients, hf.orbitalEnergies, hf.eriStore, numOccupied, mol.numBasis,
              (msg) => post({ type: 'progress', message: msg }));
            mp2Energy = mp3Result.mp2;
            mp3Energy = mp3Result.mp3;
          }
        } else {
          const mp3Result = computeMP3Energy(hf.coefficients, hf.orbitalEnergies, hf.eriStore, numOccupied, mol.numBasis,
            (msg) => post({ type: 'progress', message: msg }));
          mp2Energy = mp3Result.mp2;
          mp3Energy = mp3Result.mp3;
        }
        post({ type: 'step', id: 'mp3', status: 'done', detail: `${mp3Energy?.toFixed(6)} Eh` });
      }

      if (runCCSD) {
        post({ type: 'step', id: 'ccsd', status: 'start' });
        post({ type: 'progress', message: `Computing CCSD correlation energy${useWasmPH ? ' (WASM)' : ''}...` });
        const ccsdProgress = (msg: string) => {
          post({ type: 'progress', message: msg });
          post({ type: 'step', id: 'ccsd', status: 'update', detail: ccsdDetailFromMsg(msg) });
        };
        if (useWasmPH) {
          ccsdEnergy = computeCCSDEnergyWasm(hf.eriStore.data, hf.coefficients.data, hf.orbitalEnergies, numOccupied, mol.numBasis)
            ?? computeCCSDEnergy(hf.coefficients, hf.orbitalEnergies, hf.eriStore, numOccupied, mol.numBasis, ccsdProgress);
        } else {
          ccsdEnergy = computeCCSDEnergy(hf.coefficients, hf.orbitalEnergies, hf.eriStore, numOccupied, mol.numBasis, ccsdProgress);
        }
        post({ type: 'progress', message: `CCSD correlation energy: ${ccsdEnergy.toFixed(10)} Hartree` });
        post({ type: 'step', id: 'ccsd', status: 'done', detail: `${ccsdEnergy.toFixed(6)} Eh` });
      }

      if (runCIS) {
        post({ type: 'step', id: 'cis', status: 'start' });
        const nStates = cisNStates ?? 5;
        cisIsTriplet = excitedTriplet ?? false;
        if (dftConfig && hf instanceof RHF && hf.xcFunctional && hf.grid) {
          // TDDFT-TDA path (KS reference + XC kernel + RSH-aware exchange)
          post({ type: 'progress', message: `Computing TDDFT-TDA ${cisIsTriplet ? 'triplet' : 'singlet'} excited states (${nStates} roots)...` });
          const xcFn = hf.xcFunctional;
          const hfFrac = xcFn.exactExchangeFraction;
          const rsBeta = xcFn.rangeSeparation?.beta ?? 0;
          const tddftResult = computeTDDFT(
            hf.coefficients, hf.orbitalEnergies, hf.eriStore,
            (hf as unknown as { _eriLR?: typeof hf.eriStore })._eriLR ?? null,
            xcFn, hf.grid,
            mol.primitiveShells, mol.cgtoNormalizationFactors,
            numOccupied, mol.numBasis, nStates, cisIsTriplet,
            hfFrac, rsBeta,
            (msg) => {
              post({ type: 'progress', message: msg });
              // Strip leading prefix for compact display
              const detail = msg.replace(/^TDDFT-(TDA|Casida)( singlet| triplet)?:\s*/i, '');
              post({ type: 'step', id: 'cis', status: 'update', detail });
            },
            fullCasida ?? false,
          );
          cisStates = tddftResult.states;
          const tag = fullCasida ? 'TDDFT-Casida' : 'TDDFT-TDA';
          post({ type: 'progress', message: `${tag}: ${cisStates.length} excited states computed` });
          post({ type: 'step', id: 'cis', status: 'done', detail: `${cisStates.length} states (${tag})` });
        } else {
          // Standard CIS (HF reference)
          post({ type: 'progress', message: `Computing CIS ${cisIsTriplet ? 'triplet' : 'singlet'} excited states (${nStates} roots)...` });
          const cisResult = computeCIS(
            hf.coefficients, hf.orbitalEnergies, hf.eriStore,
            numOccupied, mol.numBasis, nStates, cisIsTriplet,
            mol.primitiveShells, mol.cgtoNormalizationFactors,
            (msg) => post({ type: 'progress', message: msg }),
          );
          cisStates = cisResult.states;
          post({ type: 'progress', message: `CIS: ${cisStates.length} excited states computed` });
          post({ type: 'step', id: 'cis', status: 'done', detail: `${cisStates.length} states` });
        }
      }

      if (runADC2) {
        post({ type: 'step', id: 'adc2', status: 'start' });
        const nStates = cisNStates ?? 5;
        const adc2Triplet = excitedTriplet ?? false;
        post({ type: 'progress', message: `Computing ADC(2) ${adc2Triplet ? 'triplet' : 'singlet'} excited states (${nStates} roots)...` });
        const adc2Result = computeADC2(
          hf.coefficients, hf.orbitalEnergies, hf.eriStore,
          numOccupied, mol.numBasis, nStates,
          mol.primitiveShells, mol.cgtoNormalizationFactors,
          (msg) => post({ type: 'progress', message: msg }),
          adc2Triplet,
        );
        cisStates = adc2Result.states;
        cisIsTriplet = adc2Triplet;
        // Use ADC(2)'s MP2 energy if MP2 wasn't requested separately
        if (!mp2Energy) mp2Energy = adc2Result.mp2Energy;
        post({ type: 'progress', message: `ADC(2): ${cisStates.length} excited states computed` });
        post({ type: 'step', id: 'adc2', status: 'done', detail: `${cisStates.length} states` });
      }

      if (runGradient) {
        post({ type: 'step', id: 'gradient', status: 'start' });
        const isDFTGrad = !!(hf.xcFunctional && hf.grid);
        const isRSHGrad = !!hf.xcFunctional?.rangeSeparation;
        post({ type: 'progress', message: `Computing nuclear gradient${isDFTGrad ? ' (KS analytic)' : ''}...` });
        if (isRSHGrad) {
          post({ type: 'progress', message: 'Note: range-separated long-range exchange gradient is not yet implemented; force may have small inaccuracy.' });
        }
        const gradResult = computeRHFGradient(
          mol.primitiveShells, mol.atoms, mol.cgtoNormalizationFactors,
          mol.numBasis, numOccupied, hf.density, hf.coefficients, hf.orbitalEnergies,
          undefined,
          isDFTGrad ? { functional: hf.xcFunctional!, grid: hf.grid! } : undefined,
        );
        gradient = Array.from(gradResult.total);
        const maxF = Math.max(...gradient.map(Math.abs));
        post({ type: 'progress', message: `Gradient computed: max |force| = ${maxF.toExponential(4)} Eh/bohr` });
        post({ type: 'step', id: 'gradient', status: 'done', detail: `max ${maxF.toExponential(2)}` });
      }

      if (runHessian) {
        post({ type: 'step', id: 'hessian', status: 'start' });
        const refCoords = new Float64Array(3 * mol.atoms.length);
        for (let i = 0; i < mol.atoms.length; i++) {
          refCoords[3*i] = mol.atoms[i].coordinate.x;
          refCoords[3*i+1] = mol.atoms[i].coordinate.y;
          refCoords[3*i+2] = mol.atoms[i].coordinate.z;
        }
        const hessResult = await computeHessianAuto(
          mol.atoms.map(a => a.atomicNumber), refCoords, basisSet, charge, 5e-4,
          (msg) => post({ type: 'progress', message: msg }),
          dftConfig,
        );
        frequencies = hessResult.frequencies;
        irIntensities = hessResult.intensities;
        zpe = hessResult.zpe;
        // Set E_elec and recompute thermo totals
        if (hessResult.thermo) {
          const th = hessResult.thermo;
          th.E_elec = totalEnergy;
          th.E_tot = totalEnergy + th.thermalCorr;
          th.H = totalEnergy + th.enthalpyCorr;
          th.G = totalEnergy + th.gibbsCorr;
          thermo = th;
        }
        const nReal = frequencies.filter(f => f > 50).length;
        post({ type: 'progress', message: `Hessian: ${nReal} vibrational modes computed` });
        post({ type: 'step', id: 'hessian', status: 'done', detail: `${nReal} modes` });
      }

      // Generate Molden file
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
        post({ type: 'step', id: 'mp2', status: 'start' });
        post({ type: 'progress', message: `Computing ROMP2 correlation energy${useWasmPH ? ' (WASM)' : ''}...` });
        if (useWasmPH) {
          mp2Energy = computeROMP2EnergyWasm(hf.eriStore.data, hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
            mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis)
            ?? computeROMP2Energy(hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
              hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis);
        } else {
          mp2Energy = computeROMP2Energy(hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
            hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis);
        }
        post({ type: 'progress', message: `MP2 correlation energy: ${mp2Energy.toFixed(10)} Hartree` });
        post({ type: 'step', id: 'mp2', status: 'done', detail: `${mp2Energy.toFixed(6)} Eh` });
      }

      if (runMP3) {
        post({ type: 'step', id: 'mp3', status: 'start' });
        post({ type: 'progress', message: `Computing ROMP3 correlation energy${useWasmPH ? ' (WASM)' : ''}...` });
        if (useWasmPH) {
          const wasmResult = computeROMP3EnergyWasm(hf.eriStore.data, hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
            mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis);
          if (wasmResult) {
            mp2Energy = wasmResult.mp2;
            mp3Energy = wasmResult.mp3;
          } else {
            const mp3Result = computeROMP3Energy(hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
              hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis,
              (msg) => post({ type: 'progress', message: msg }));
            mp2Energy = mp3Result.mp2;
            mp3Energy = mp3Result.mp3;
          }
        } else {
          const mp3Result = computeROMP3Energy(hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
            hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis,
            (msg) => post({ type: 'progress', message: msg }));
          mp2Energy = mp3Result.mp2;
          mp3Energy = mp3Result.mp3;
        }
        post({ type: 'step', id: 'mp3', status: 'done', detail: `${mp3Energy?.toFixed(6)} Eh` });
      }

      if (runCCSD) {
        post({ type: 'step', id: 'ccsd', status: 'start' });
        post({ type: 'progress', message: 'Computing ROCCSD correlation energy...' });
        const ccsdProgress = (msg: string) => {
          post({ type: 'progress', message: msg });
          post({ type: 'step', id: 'ccsd', status: 'update', detail: ccsdDetailFromMsg(msg) });
        };
        ccsdEnergy = computeROCCSDEnergy(hf.coefficients, hf.fockAlphaMatrix, hf.fockBetaMatrix,
          hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis, ccsdProgress);
        post({ type: 'progress', message: `CCSD correlation energy: ${ccsdEnergy.toFixed(10)} Hartree` });
        post({ type: 'step', id: 'ccsd', status: 'done', detail: `${ccsdEnergy.toFixed(6)} Eh` });
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
        post({ type: 'step', id: 'mp2', status: 'start' });
        post({ type: 'progress', message: `Computing UMP2 correlation energy${useWasmPH ? ' (WASM)' : ''}...` });
        if (useWasmPH) {
          mp2Energy = computeUMP2EnergyWasm(hf.eriStore.data, hf.coefficientsAlpha.data, hf.coefficientsBeta.data,
            hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
            mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis)
            ?? computeUMP2Energy(hf.coefficientsAlpha, hf.coefficientsBeta,
              hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
              hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis);
        } else {
          mp2Energy = computeUMP2Energy(hf.coefficientsAlpha, hf.coefficientsBeta,
            hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
            hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis);
        }
        post({ type: 'progress', message: `MP2 correlation energy: ${mp2Energy.toFixed(10)} Hartree` });
        post({ type: 'step', id: 'mp2', status: 'done', detail: `${mp2Energy.toFixed(6)} Eh` });
      }

      if (runMP3) {
        post({ type: 'step', id: 'mp3', status: 'start' });
        post({ type: 'progress', message: `Computing UMP3 correlation energy${useWasmPH ? ' (WASM)' : ''}...` });
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
              hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis,
              (msg) => post({ type: 'progress', message: msg }));
            mp2Energy = mp3Result.mp2;
            mp3Energy = mp3Result.mp3;
          }
        } else {
          const mp3Result = computeUMP3Energy(hf.coefficientsAlpha, hf.coefficientsBeta,
            hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
            hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis,
            (msg) => post({ type: 'progress', message: msg }));
          mp2Energy = mp3Result.mp2;
          mp3Energy = mp3Result.mp3;
        }
        post({ type: 'step', id: 'mp3', status: 'done', detail: `${mp3Energy?.toFixed(6)} Eh` });
      }

      if (runCCSD) {
        post({ type: 'step', id: 'ccsd', status: 'start' });
        post({ type: 'progress', message: 'Computing UCCSD correlation energy...' });
        const ccsdProgress = (msg: string) => {
          post({ type: 'progress', message: msg });
          post({ type: 'step', id: 'ccsd', status: 'update', detail: ccsdDetailFromMsg(msg) });
        };
        ccsdEnergy = computeUCCSDEnergy(hf.coefficientsAlpha, hf.coefficientsBeta,
          hf.orbitalEnergiesAlpha, hf.orbitalEnergiesBeta,
          hf.eriStore, mol.numAlphaSpins, mol.numBetaSpins, mol.numBasis, ccsdProgress);
        post({ type: 'progress', message: `CCSD correlation energy: ${ccsdEnergy.toFixed(10)} Hartree` });
        post({ type: 'step', id: 'ccsd', status: 'done', detail: `${ccsdEnergy.toFixed(6)} Eh` });
      }

      moldenText = writeMolden({
        atoms: mol.atoms, basisSet, coefficients: hf.coefficientsAlpha,
        orbitalEnergies: hf.orbitalEnergiesAlpha, numOccupied, numBasis: mol.numBasis,
        coefficientsBeta: hf.coefficientsBeta,
        orbitalEnergiesBeta: hf.orbitalEnergiesBeta,
        numOccupiedBeta: mol.numBetaSpins,
      });
    }

    // New analysis properties
    post({ type: 'step', id: 'properties', status: 'start' });
    let s2: { exact: number; computed: number } | undefined;
    let lowdinCharges: number[] = [];
    let bondOrders: { i: number; j: number; order: number }[] = [];
    let energyComponents: { oneElectron: number; twoElectron: number; kinetic: number; nuclearAttraction: number } | undefined;

    if (hf instanceof RHF) {
      lowdinCharges = Array.from(computeLowdinCharges(hf.density, hf.overlap, mol.atoms, mol.atomToBasisRange));
      // RHF: Pα = Pβ = P/2
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

    let densityMatrix: number[] | undefined;
    let densityMatrixAlpha: number[] | undefined;
    let densityMatrixBeta: number[] | undefined;
    if (hf instanceof RHF) {
      densityMatrix = Array.from(hf.density.data) as number[];
    } else if (hf instanceof ROHF) {
      densityMatrix = Array.from(hf.density.data) as number[];
      densityMatrixAlpha = Array.from(hf.densityAlphaMatrix.data) as number[];
      densityMatrixBeta = Array.from(hf.densityBetaMatrix.data) as number[];
    } else if (hf instanceof UHF) {
      densityMatrix = Array.from(hf.density.data) as number[];
      densityMatrixAlpha = Array.from(hf.densityAlphaMatrix.data) as number[];
      densityMatrixBeta = Array.from(hf.densityBetaMatrix.data) as number[];
    }

    // UHF beta orbital data
    let orbitalEnergiesBeta: number[] | undefined;
    let numOccupiedBeta: number | undefined;
    if (hf instanceof UHF) {
      orbitalEnergiesBeta = Array.from(hf.orbitalEnergiesBeta);
      numOccupiedBeta = mol.numBetaSpins;
    } else if (hf instanceof ROHF) {
      numOccupiedBeta = mol.numBetaSpins;
    }

    post({ type: 'step', id: 'properties', status: 'done' });

    // ── Dispersion correction (post-SCF) ──
    // Backwards compatible: legacy `runD2` boolean still accepted; new
    // `dispersion` field selects between 'd2' and 'd3bj'.
    let d2Energy: number | undefined;
    let dispersionEnergy: number | undefined;
    let dispersionMethodOut: DispersionMethod | undefined;
    const dispMethod: DispersionMethod = dispersion ?? (runD2 ? 'd2' : 'none');
    if (dispMethod !== 'none') {
      const functionalKey = dftConfig?.functional ?? method;
      const r = computeDispersion(mol.atoms, dispMethod, functionalKey);
      dispersionEnergy = r.energy;
      dispersionMethodOut = r.method;
      if (dispMethod === 'd2') d2Energy = r.energy;  // legacy field
      const params = dispMethod === 'd2'
        ? `s6=${r.s6}`
        : `s6=${r.s6}, s8=${r.s8}, a1=${r.a1}, a2=${r.a2}`;
      post({ type: 'progress', message: `${dispMethod.toUpperCase()} dispersion (${params}): ${r.energy.toFixed(8)} Eh` });
    }

    post({
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
      d2Energy,
      dispersionEnergy,
      dispersionMethod: dispersionMethodOut,
      moldenText,
      orbitalEnergiesBeta,
      numOccupiedBeta,
      s2,
      lowdinCharges,
      bondOrders,
      energyComponents,
      cisStates,
      excitedTriplet: cisIsTriplet,
      gradient,
      frequencies,
      irIntensities,
      zpe,
      thermo,
    });
  } catch (err) {
    self.postMessage({
      type: 'error',
      message: err instanceof Error ? err.message : String(err),
    } as WorkerError);
  }
};
