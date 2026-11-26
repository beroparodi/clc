const M = require('../src/model.js');
const iso = (s) => Date.parse(s + 'T00:00:00Z');
const un = { q: 10.94525581, e: 0.99869477, i: 95.55335892, Om: 190.0846361196, w: 326.33161880945, Tp: 2462892.9068, full: true };
for (const d of ['2021-12-15', '2024-02-15', '2026-10-01', '2028-07-01', '2031-01-25']) {
  const g = M.ephem(un, iso(d)); console.log(d, g.r.toFixed(4), g.D.toFixed(4), g.alpha.toFixed(3));
}
const s = { HN: 8.40, beta: 0.04, j: 1, ap: 6, msun: -26.93 };
const K = 4110019.53, mu = M.muFrom(K, s); console.log('mu', mu.toFixed(2));
for (const d of ['2026-10-01', '2028-07-01', '2031-01-25', '2021-12-15', '2024-02-15']) {
  const g = M.ephem(un, iso(d)); const r = M.magnitude(g.r, g.D, g.alpha, M.eta(g.r, 'CO2', 1), s, mu); console.log(d, r.m.toFixed(3), r.mN.toFixed(3));
}
for (const r of [1, 3, 10.95, 20]) console.log(r, ['H2O', 'CO2', 'CO'].map(x => [x, M.balance(r, x, 1).T.toFixed(2), M.eta(r, x, 1).toExponential(4), M.eta(r, x, 0.25).toExponential(4)].join(' ')).join(' | '));
// T5 (hyperbolic) vs Python Barker approximation
const t5 = { q: 3.84251371, e: 1.00240368, i: 52.38681511, Om: 100.67985066075, w: 352.41723883275, Tp: M.jdOf(Date.parse('2027-05-06T02:38:00Z')), full: true };
for (const d of ['2025-11-02', '2026-11-10', '2027-05-06', '2027-11-01']) { const g = M.ephem(t5, iso(d)); console.log('T5', d, g.r.toFixed(3), g.D.toFixed(3), g.alpha.toFixed(1), g.elong.toFixed(0)); }
// Y1 parabolic r
const y1 = { q: 2.6, e: 1, Tp: M.jdOf(iso('2027-10-15')), full: false };
for (const d of ['2025-12-26', '2026-12-20', '2027-10-15']) console.log('Y1', d, M.ephem(y1, iso(d)).r.toFixed(3));
// growth factors (paper: Y1 from 4 au: 1.9, 2.0, 3.7, 1300; T5 from 4.14 au: 1.12 1.13 1.37 4.7)
const gr = (r0, r1, x, c) => M.eta(r1, x, c) / M.eta(r0, x, c) * Math.pow(r1 / r0, -1.5);
console.log('Y1', [['CO',1],['CO2',1],['H2O',1],['H2O',0.25]].map(([x,c]) => gr(4, 2.6, x, c).toPrecision(3)).join(' '));
const r0 = M.ephem(t5, iso('2026-11-10')).r;
console.log('T5 r0', r0.toFixed(3), [['CO',1],['CO2',1],['H2O',1],['H2O',0.25]].map(([x,c]) => gr(r0, 3.8425, x, c).toPrecision(3)).join(' '));
// n_eff at 6, 8, 10 au for CO2 isothermal (paper 4.7, 6.4, 10.7)
console.log('neff CO2 iso', [6, 8, 10, 12, 15].map(r => (2 + M.kSlope(r, 'CO2', 0.25)).toFixed(1)).join(' '));
console.log('neff CO2 sub', [6, 8, 10].map(r => (2 + M.kSlope(r, 'CO2', 1)).toFixed(2)).join(' '), 'CO sub', [6, 8, 10].map(r => (2 + M.kSlope(r, 'CO', 1)).toFixed(2)).join(' '));
// power budget
console.log('NEOWISE f(3.5km)', (5.27e30 / M.qMax(3.5, 0.337, 'H2O', 1)).toFixed(3), '3I f(1.3km)', (4e28 / M.qMax(1.3, 1.49, 'H2O', 1)).toFixed(2), (5.15e28 / M.qMax(1.3, 1.49, 'H2O', 1)).toFixed(2));
console.log('envelope q=1', M.envelope(1, 0.25).toFixed(2), M.envelope(1, 1).toFixed(2));
// a short-period ellipse sanity check: 67P-like orbit r at perihelion and aphelion
const p67 = { q: 1.243, e: 0.641, i: 7.04, Om: 50.1, w: 12.8, Tp: M.jdOf(iso('2015-08-13')), full: true };
const P = 2 * Math.PI / (0.01720209895 / Math.pow(1.243 / (1 - 0.641), 1.5));
console.log('67P r(T)', M.ephem(p67, iso('2015-08-13')).r.toFixed(3), 'r(T+P/2)', M.ephem(p67, iso('2015-08-13') + P / 2 * 86400000).r.toFixed(3), 'Q', (1.243 * (1 + 0.641) / (1 - 0.641)).toFixed(3));
// dust mantle (paper: 89 of 123 pass; thick mantle at 1 au Ts ~ 403 K, Ti ~ 174 K; envelope q=1 [-3.20, -1.01])
const bm = M.balance(1, 'H2O', 1, 0.03); console.log('mantle Ts Ti', bm.T.toFixed(0), bm.Ti.toFixed(0), 'env q=1', M.envelopeMantle(1).map(v => v.toFixed(2)).join(' '), 'q=0.1', M.envelopeMantle(0.1).map(v => v.toFixed(2)).join(' '));
// two ices on separate patches (Eq. 14): k = 1 gives the single ice; near full activity n is close to 7/2; water's dust share for C/2019 L3 (paper: 4 % at 7 au, 90 % at 4 au)
{
  const one = [{ x: 'CO', w: 1 }], mix = [{ x: 'CO', w: 0.04 }, { x: 'H2O', w: 0.96 }];
  const d1 = Math.abs(M.etaMix(5, one, 1) - M.eta(5, 'CO', 1)), d2 = Math.abs(M.kSlopeMix(5, [{ x: 'CO', w: 1 }, { x: 'H2O', w: 0 }], 1) - M.kSlope(5, 'CO', 1));
  const nIn = 2 + M.kSlopeMix(1, mix, 1), nOut = [7, 5.5, 4].map((r) => (2 + M.kSlopeMix(r, mix, 1)).toFixed(2));
  const w = (r) => (0.96 * M.eta(r, 'H2O', 1) / M.etaMix(r, mix, 1)).toFixed(2);
  console.log('mix: single-ice limits', d1 < 1e-12 && d2 < 1e-9 ? 'ok' : 'FAIL', '| n at 1 au', nIn.toFixed(3), '| n at 7, 5.5, 4 au', nOut.join(' '), '| water dust share at 7, 4 au', w(7), w(4));
  if (!(d1 < 1e-12 && d2 < 1e-9 && Math.abs(nIn - 3.5) < 0.15)) process.exitCode = 1;
}
