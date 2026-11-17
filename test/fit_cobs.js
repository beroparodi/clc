// Fits the coma scale K of the bundled COBS light curves (median visual magnitudes), as the tool's
// Fit button does, for the five ice scenarios; run: node test/fit_cobs.js
const M = require('../src/model.js'), fs = require('fs');
const iso = (s) => Date.parse(s.length > 10 ? s.replace(' ', 'T') + ':00Z' : s + 'T00:00:00Z');
const EL = {
  hyak: { q: 0.23022353, e: 0.99989165, i: 124.92204939, Om: 188.04513199216, w: 130.1751209781, Tp: 2450204.8941, HN: 15.71 },
  hb: { q: 0.9141, e: 0.9951, i: 89.43, Om: 282.47, w: 130.59, Tp: M.jdOf(iso('1997-04-01 03:19')), HN: 11.39 },
  neo: { q: 0.295, e: 0.99921, i: 128.94, Om: 61.010, w: 37.279, Tp: M.jdOf(iso('2020-07-03 12:00')), HN: 15 },
  p12: { q: 0.78086475, e: 0.95456023, i: 74.19083929, Om: 255.85505059528, w: 198.98800614657, Tp: 2460421.6316, HN: 13 },
  a3: { q: 0.39145906, e: 1.00012475, i: 139.11502467, Om: 21.558609111275, w: 308.4872848828, Tp: 2460581.2353, HN: 15 },
  i3: { q: 1.35638454, e: 6.13941774, i: 175.11310480, Om: 322.15684249105, w: 128.01051366675, Tp: 2460977.9814, HN: 17.04 },
  t5: { q: 3.84251371, e: 1.00240368, i: 52.38681511, Om: 100.67985066075, w: 352.41723883275, Tp: M.jdOf(iso('2027-05-06 02:38')), HN: 15 },
  y1: { q: 2.601924069, e: 0.999552254, i: 14.55698337, Om: 166.5518278289, w: 124.2978493120, Tp: 2461692.5495, HN: 15 },
  sw1: { q: 5.77512802, e: 0.04479814, i: 9.36477141, Om: 312.39091654503, w: 50.662700071797, Tp: 2458589.9259, HN: 10.22 },
};
const SCEN = [['CO', 1], ['CO2', 1], ['CO2', 0.25], ['H2O', 1], ['H2O', 0.25]];
const data = {}; let cur = null;
for (const l of fs.readFileSync(process.argv[2], 'utf8').split('\n')) {
  if (l.startsWith('#')) { cur = l.slice(1).trim(); data[cur] = []; } else if (l.trim() && !l.startsWith('%')) { const [d, m, s] = l.split(/\s+/); data[cur].push({ ms: iso(d), m: +m, s: +s }); }
}
function fit(key, x, c, HN) {
  const el = { ...EL[key], full: true }, P = { HN, beta: 0.04, j: 0, ap: 20000, msun: -26.76 };
  const pts = data[key].map((o) => { const g = M.ephem(el, o.ms); return { ...o, r: g.r, D: g.D, al: g.alpha, et: M.eta(g.r, x, c) }; });
  const chi = (lk) => { const mu = M.muFrom(10 ** lk, P); let s2 = 0; for (const p of pts) { const d = (p.m - M.magnitude(p.r, p.D, p.al, p.et, P, mu).m) / p.s; s2 += d * d; } return s2; };
  let bl = 0, bv = Infinity; for (let lk = -1; lk <= 10; lk += 0.05) { const v = chi(lk); if (v < bv) { bv = v; bl = lk; } }
  let a = bl - 0.05, b = bl + 0.05; const gr = (Math.sqrt(5) - 1) / 2;
  for (let k = 0; k < 40; k++) { const c1 = b - gr * (b - a), c2 = a + gr * (b - a); if (chi(c1) < chi(c2)) b = c2; else a = c1; }
  const lk = (a + b) / 2; return { lk, chi2nu: chi(lk) / (pts.length - 1) };
}
for (const key of Object.keys(data).filter((k) => EL[k])) {
  const line = SCEN.map(([x, c]) => { const r = fit(key, x, c, EL[key].HN); return `${x}${c === 1 ? '' : '/iso'} lk=${r.lk.toFixed(3)} chi2nu=${r.chi2nu.toFixed(1)}`; });
  console.log(key, data[key].length, '\n  ' + line.join('\n  '));
}
