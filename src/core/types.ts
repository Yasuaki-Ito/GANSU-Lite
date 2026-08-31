/** Core data structures mirroring GANSU C++ types.hpp */

export interface Vec3 {
  x: number;
  y: number;
  z: number;
}

export interface Atom {
  atomicNumber: number;
  coordinate: Vec3;
  atomIndex: number;
}

export interface PrimitiveShell {
  exponent: number;
  coefficient: number;
  coordinate: Vec3;
  shellType: number; // 0=s, 1=p, 2=d, 3=f, ...
  basisIndex: number;
  atomIndex: number;
}

export interface ShellTypeInfo {
  count: number;
  startIndex: number;
}

export interface BasisRange {
  startIndex: number;
  endIndex: number;
}

export interface PrimitiveGauss {
  exponent: number;
  coefficient: number;
}
