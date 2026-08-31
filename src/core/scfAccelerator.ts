/** SCF convergence accelerator — Strategy Pattern interface */

import type { FloatArray } from '../linalg/matrix';

export type SCFAccelMethod = 'diis' | 'damping' | 'level-shift' | 'ediis' | 'adiis-diis';

/** Per-method tuneable parameters */
export interface SCFAccelParams {
  diisHistory?: number;        // DIIS history size (default 8)
  dampingAlpha?: number;       // Damping mixing ratio (default 0.5)
  levelShift?: number;         // Level shift in Hartree (default 0.5)
  ediisHistory?: number;       // EDIIS history size (default 8)
  adiisHistory?: number;       // ADIIS history size (default 8)
  adiisThreshold?: number;     // ADIIS→DIIS switch threshold (default 0.1)
}

export interface SCFAccelInput {
  fock: FloatArray[];      // per-channel (RHF=[1], UHF=[alpha,beta])
  density: FloatArray[];   // per-channel
  error: FloatArray[];     // FPS-SPF commutator per-channel
  energy: number;
  overlap: FloatArray;
  numBasis: number;
  iteration: number;
}

export interface SCFAccelOutput {
  fock?: FloatArray[];     // updated Fock matrix
}

export interface SCFAccelerator {
  accelerate(input: SCFAccelInput): SCFAccelOutput;
  reset(): void;
  readonly name: string;
}
