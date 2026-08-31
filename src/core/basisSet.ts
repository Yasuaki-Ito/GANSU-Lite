/** Basis set classes — port of GANSU basis_set.hpp + basis_set.cpp */

import type { PrimitiveGauss } from './types';
import { shellNameToType, ANGULAR_MOMENTUMS } from './constants';

export class ContractedGauss {
  readonly type: string; // "S", "P", "D", ...
  readonly primitives: PrimitiveGauss[] = [];

  constructor(type: string) {
    this.type = type;
  }

  addPrimitive(exponent: number, coefficient: number) {
    this.primitives.push({ exponent, coefficient });
  }

  /** Compute CGTO normalization factors for each Cartesian component.
   *  Order matches AngularMomentums for this shell type. */
  getNormalizationFactors(): number[] {
    const shellType = shellNameToType(this.type);
    const angularMomentumList = ANGULAR_MOMENTUMS[shellType];
    const factors: number[] = [];

    for (const am of angularMomentumList) {
      const lx = am[0], ly = am[1], lz = am[2];
      const L = lx + ly + lz;
      let norm = 0.0;

      for (const a of this.primitives) {
        const nA = Math.pow(2.0 * Math.PI, -3.0 / 4.0) * Math.pow(4.0 * a.exponent, L / 2.0 + 0.75);
        for (const b of this.primitives) {
          const nB = Math.pow(2.0 * Math.PI, -3.0 / 4.0) * Math.pow(4.0 * b.exponent, L / 2.0 + 0.75);
          norm += a.coefficient * nA * b.coefficient * nB
            / Math.pow(2.0 * (a.exponent + b.exponent), L + 1.5);
        }
      }

      norm *= Math.pow(2.0 * Math.PI, 1.5);
      factors.push(Math.pow(norm, -0.5));
    }

    return factors;
  }
}

export class ElementBasisSet {
  elementName = '';
  readonly contractedGausses: ContractedGauss[] = [];

  addContractedGauss(cg: ContractedGauss) {
    this.contractedGausses.push(cg);
  }
}

export class BasisSet {
  private readonly elementBasisSets = new Map<string, ElementBasisSet>();

  add(ebs: ElementBasisSet) {
    this.elementBasisSets.set(ebs.elementName, ebs);
  }

  get(elementName: string): ElementBasisSet {
    const ebs = this.elementBasisSets.get(elementName);
    if (!ebs) throw new Error(`Basis set does not include element: ${elementName}`);
    return ebs;
  }

  /** Parse a Gaussian basis set (.gbs) file text */
  static fromGBS(text: string): BasisSet {
    const basisSet = new BasisSet();
    const lines = text.split(/\r?\n/);
    let i = 0;

    // Skip header lines (comments starting with '!' or non-alpha first char)
    while (i < lines.length) {
      const line = lines[i].trim();
      if (line.length > 0 && /^[A-Za-z]/.test(line)) break;
      i++;
    }

    while (i < lines.length) {
      const elementLine = lines[i]?.trim();
      if (!elementLine || elementLine === '****') { i++; continue; }

      const ebs = new ElementBasisSet();
      const parts = elementLine.split(/\s+/);
      ebs.elementName = parts[0];
      i++;

      // Read contracted gaussians until '****'
      while (i < lines.length) {
        const line = lines[i]?.trim();
        if (!line || line === '****') { i++; break; }

        const header = line.split(/\s+/);
        const type = header[0]; // "S", "P", "SP", "D", ...
        const numPrimitives = parseInt(header[1], 10);
        i++;

        if (type.length === 1) {
          // Single type (S, P, D, F, ...)
          const cg = new ContractedGauss(type);
          for (let j = 0; j < numPrimitives; j++) {
            const primLine = lines[i]?.replace(/D/g, 'E').trim();
            i++;
            const vals = primLine.split(/\s+/);
            cg.addPrimitive(parseFloat(vals[0]), parseFloat(vals[1]));
          }
          ebs.addContractedGauss(cg);
        } else if (type.length === 2) {
          // Dual type (SP, etc.)
          const cg0 = new ContractedGauss(type[0]);
          const cg1 = new ContractedGauss(type[1]);
          for (let j = 0; j < numPrimitives; j++) {
            const primLine = lines[i]?.replace(/D/g, 'E').trim();
            i++;
            const vals = primLine.split(/\s+/);
            cg0.addPrimitive(parseFloat(vals[0]), parseFloat(vals[1]));
            cg1.addPrimitive(parseFloat(vals[0]), parseFloat(vals[2]));
          }
          ebs.addContractedGauss(cg0);
          ebs.addContractedGauss(cg1);
        } else {
          throw new Error(`Invalid basis function type: ${type}`);
        }
      }

      basisSet.add(ebs);
    }

    return basisSet;
  }
}
