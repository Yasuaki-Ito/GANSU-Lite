/** Device fingerprint recorded alongside every benchmark and stress-test run. */
export interface DeviceInfo {
  ua: string;
  platform: string;
  deviceMemoryGB: number | null;
  hardwareConcurrency: number;
  screen: string;
  language: string;
  timestamp: string;
}

export function gatherDeviceInfo(): DeviceInfo {
  const nav = navigator as Navigator & { deviceMemory?: number };
  return {
    ua: nav.userAgent,
    platform: nav.platform ?? 'unknown',
    deviceMemoryGB: nav.deviceMemory ?? null,
    hardwareConcurrency: nav.hardwareConcurrency ?? 0,
    screen: `${screen.width}×${screen.height}`,
    language: nav.language || 'unknown',
    timestamp: new Date().toISOString(),
  };
}

/** Full page load time from the navigation timing entry, in ms. */
export function getInitialLoadMs(): number | null {
  try {
    const nav = performance.getEntriesByType('navigation')[0] as PerformanceNavigationTiming | undefined;
    if (!nav) return null;
    return Math.round(nav.loadEventEnd - nav.fetchStart);
  } catch {
    return null;
  }
}
