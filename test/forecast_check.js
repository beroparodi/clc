// Forecast test of the second ten comets with the tool's own model code: K and the ice from the
// first half (as in the presets), empirical law fitted to the first half; run: node test/forecast_check.js
const M = require('../src/model.js'), fs = require('fs');
const iso = (s) => Date.parse(s + 'T00:00:00Z');
const P = (q, e, i, Om, w, Tp, ice, cz, lk, cut) => ({ el: { q, e, i, Om, w, Tp, full: true }, ice, cz, lk, cut });
const C = {
  q2: P(1.290358238881392, 0.9977738520356274, 80.30130181438155, 94.97505909428631, 12.39526184803142, 2457052.569249606, 'H2O', 0.25, 5.061, '2015-04-18'),
  k2: P(1.797281961376049, 1.000616562373503, 87.55487509681637, 88.24109127986824, 236.1923489989567, 2459933.227066954, 'CO', 1, 5.963, '2021-10-24'),
  t2: P(1.615095955888523, 0.9996765980850288, 57.23142385337199, 64.37877215685057, 92.99331253710214, 2458974.449996825, 'CO2', 1, 4.770, '2019-10-30'),
  p1: P(1.550535584223489, 1.000997435576857, 106.1772916768743, 325.9972285528648, 90.74706541748093, 2455919.177555371, 'CO', 1, 5.548, '2011-10-20'),
  e3: P(1.11224437022534, 1.000301905819192, 109.169480756749, 302.5550197168474, 145.81492879, 2459957.285198830, 'H2O', 1, 4.070, '2022-10-22'),
  l3: P(3.554492702308672, 1.001405275505796, 48.36374698279935, 290.7888867057587, 171.6137051016775, 2459589.132623246, 'CO2', 0.25, 5.308, '2021-09-11'),
  l4: P(0.3015442238164768, 1.000032767108972, 84.20819702923397, 65.66588639150341, 333.6516376633611, 2456361.669960002, 'CO2', 0.25, 5.019, '2013-04-02'),
  a6: P(0.5299041302912412, 0.9955362293967411, 143.663707214442, 108.0979873301086, 132.9688404751593, 2460988.039387102, 'CO', 1, 4.252, '2026-02-06'),
  v2j: P(1.636950053071025, 1.001702614530096, 49.87481167153636, 69.8521003063327, 164.8980978916873, 2457916.847458616, 'CO', 1, 4.644, '2017-01-27'),
  v2z: P(2.227959873700067, 1.00089953037377, 131.6110418896687, 212.3712351547342, 162.4243172377264, 2460073.056589952, 'CO2', 0.25, 5.179, '2023-02-15'),
};
const data = {}; let cur = null;
for (const l of fs.readFileSync(__dirname + '/../data/cobs_lightcurves.txt', 'utf8').split('\n')) {
  if (l.startsWith('#')) { cur = l.slice(1).trim(); data[cur] = []; } else if (l.trim() && !l.startsWith('%')) { const [d, m, s] = l.split(/\s+/); data[cur].push({ ms: iso(d), m: +m, s: +s }); }
}
const S = { HN: 15, beta: 0.04, j: 0, ap: 20000, msun: -26.76 };
for (const [k, c] of Object.entries(C)) {
  const cut = iso(c.cut), pts = data[k].map((o) => { const g = M.ephem(c.el, o.ms); return { ...o, r: g.r, D: g.D, al: g.alpha }; });
  const fit = pts.filter((p) => p.ms < cut), test = pts.filter((p) => p.ms >= cut);
  let S0 = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
  for (const p of fit) { const w = 1 / p.s ** 2, x = 2.5 * Math.log10(p.r), y = p.m - 5 * Math.log10(p.D); S0 += w; Sx += w * x; Sy += w * y; Sxx += w * x * x; Sxy += w * x * y; }
  const n = (S0 * Sxy - Sx * Sy) / (S0 * Sxx - Sx * Sx), m0 = (Sy - n * Sx) / S0, mu = M.muFrom(10 ** c.lk, S);
  const rms = (f) => Math.sqrt(test.reduce((a, p) => a + (p.m - f(p)) ** 2, 0) / test.length);
  const rm = rms((p) => M.magnitude(p.r, p.D, p.al, M.eta(p.r, c.ice, c.cz), S, mu).m), re = rms((p) => m0 + 5 * Math.log10(p.D) + 2.5 * n * Math.log10(p.r));
  console.log(k, fit.length, test.length, 'model', rm.toFixed(2), 'empirical', re.toFixed(2), 'n', n.toFixed(1));
}
