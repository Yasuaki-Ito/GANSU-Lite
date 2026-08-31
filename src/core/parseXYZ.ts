/** XYZ file parser — port of GANSU parseXYZ.cpp */

import type { Atom } from './types';
import { angstromToBohr, elementNameToAtomicNumber } from './constants';

export function parseXYZ(text: string): Atom[] {
  const lines = text.split(/\r?\n/);
  if (lines.length < 3) throw new Error('XYZ file too short');

  const atomCount = parseInt(lines[0].trim(), 10);
  if (isNaN(atomCount)) throw new Error('Invalid atom count');
  // Line 1 is comment (ignored)

  const atoms: Atom[] = [];
  const DUPLICATE_THRESHOLD = 1e-4; // in bohr

  for (let i = 2; i < 2 + atomCount; i++) {
    const line = lines[i]?.trim();
    if (!line) throw new Error(`Missing atom data at line ${i + 1}`);

    const parts = line.split(/\s+/);
    if (parts.length < 4) throw new Error(`Invalid line format at line ${i + 1}: ${line}`);

    const symbol = parts[0];
    const x = angstromToBohr(parseFloat(parts[1]));
    const y = angstromToBohr(parseFloat(parts[2]));
    const z = angstromToBohr(parseFloat(parts[3]));

    // Check for duplicate coordinates
    for (const existing of atoms) {
      const dx = x - existing.coordinate.x;
      const dy = y - existing.coordinate.y;
      const dz = z - existing.coordinate.z;
      if (dx * dx + dy * dy + dz * dz < DUPLICATE_THRESHOLD * DUPLICATE_THRESHOLD) {
        throw new Error(`Overlapping atoms at line ${i + 1}`);
      }
    }

    atoms.push({
      atomicNumber: elementNameToAtomicNumber(symbol),
      coordinate: { x, y, z },
      atomIndex: atoms.length,
    });
  }

  if (atoms.length !== atomCount) {
    throw new Error(`Atom count mismatch: expected ${atomCount}, got ${atoms.length}`);
  }

  return atoms;
}
