/** Factory for SCF accelerators */

import type { SCFAccelMethod, SCFAccelParams, SCFAccelerator } from './scfAccelerator';
import { AccelDIIS } from './accelDIIS';
import { AccelDamping } from './accelDamping';
import { AccelLevelShift } from './accelLevelShift';
import { AccelEDIIS } from './accelEDIIS';
import { AccelADIIS } from './accelADIIS';

export function createAccelerator(method: SCFAccelMethod, params?: SCFAccelParams): SCFAccelerator {
  switch (method) {
    case 'diis': return new AccelDIIS(params?.diisHistory);
    case 'damping': return new AccelDamping(params?.dampingAlpha);
    case 'level-shift': return new AccelLevelShift(params?.levelShift);
    case 'ediis': return new AccelEDIIS(params?.ediisHistory);
    case 'adiis-diis': return new AccelADIIS(params?.adiisHistory, params?.adiisThreshold);
  }
}
