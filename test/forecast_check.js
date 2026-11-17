// Forecast test of the comets that carry a cut (the fifteen added to the paper's test, Kohoutek and ISON) with the tool's own model code: the presets carry K, the ice
// and the cut of the paper (first half); the empirical law is fitted to the first half; run: node test/forecast_check.js
const M = require('../src/model.js'), fs = require('fs');
const html = fs.readFileSync(__dirname + '/../index.html', 'utf8');
const iso = (s) => Date.parse(s + 'T00:00:00Z');
const data = {}; let cur = null;
for (const l of fs.readFileSync(__dirname + '/../data/cobs_lightcurves.txt', 'utf8').split('\n')) {
  if (l.startsWith('#')) { cur = l.slice(1).trim(); data[cur] = []; } else if (l.trim() && !l.startsWith('%')) { const [d, m, s] = l.split(/\s+/); data[cur].push({ ms: iso(d), m: +m, s: +s }); }
}
const re = /  (\w+): MORE\(\{ q: "([\d.]+)", e: "([\d.]+)", i: "([\d.]+)", Om: "([\d.]+)", w: "([\d.]+)", T: fmtDT\(msOf\(([\d.]+)\)\) \},\n    \{ lk: ([\d.]+), ice: "(\w+)", cz: ([\d.]+), HN: 15, y0: \d+, y1: \d+, ref: "[\d-]+", cut: "([\d-]+)"/g;
const S = { HN: 15, beta: 0.04, j: 0, ap: 20000, msun: -26.76 };
for (const m of html.matchAll(re)) {
  const [, k, q, e, i, Om, w, Tp, lk, ice, cz, cutS] = m;
  const el = { q: +q, e: +e, i: +i, Om: +Om, w: +w, Tp: +Tp, full: true }, cut = iso(cutS);
  const pts = data[k].map((o) => { const g = M.ephem(el, o.ms); return { ...o, r: g.r, D: g.D, al: g.alpha }; });
  const fit = pts.filter((p) => p.ms < cut), test = pts.filter((p) => p.ms >= cut);
  let S0 = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
  for (const p of fit) { const wt = 1 / p.s ** 2, x = 2.5 * Math.log10(p.r), y = p.m - 5 * Math.log10(p.D); S0 += wt; Sx += wt * x; Sy += wt * y; Sxx += wt * x * x; Sxy += wt * x * y; }
  const n = (S0 * Sxy - Sx * Sy) / (S0 * Sxx - Sx * Sx), m0 = (Sy - n * Sx) / S0, mu = M.muFrom(10 ** +lk, S);
  const rms = (f) => Math.sqrt(test.reduce((a, p) => a + (p.m - f(p)) ** 2, 0) / test.length);
  console.log(k, fit.length, test.length, 'model', rms((p) => M.magnitude(p.r, p.D, p.al, M.eta(p.r, ice, +cz), S, mu).m).toFixed(2),
    'empirical', rms((p) => m0 + 5 * Math.log10(p.D) + 2.5 * n * Math.log10(p.r)).toFixed(2), 'n', n.toFixed(1));
}
